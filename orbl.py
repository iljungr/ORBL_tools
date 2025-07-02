#!/usr/bin/env python
# Copyright 2025 Irwin Jungreis
# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at http://www.apache.org/licenses/LICENSE-2.0
# Unless required by applicable law or agreed to in writing, software
# distributed under the License is distributed on an "AS IS" BASIS,
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
# See the License for the specific language governing permissions and
# limitations under the License.
"""
Utility to calculate ORF Relative Branch Length (ORBL) scores for measuring ORFness
conservation and constraint.
"""
from __future__ import division, print_function
import sys, os, itertools, subprocess, shlex
ThisDir = os.path.dirname(os.path.realpath(__file__))
sys.path.append(os.path.dirname(ThisDir)) # ORBL_tools's parent
from ORBL_tools.IntervalUtils import regionString_to_triples, get_intervals_length
from ORBL_tools.CmdLineUtils import check_arg, assert_num_args
from ORBL_tools.CalcORFrelBLs import calc_ORF_relBLs
from ORBL_tools.DownloadLocalAlignment import (
    download_local_alignment, download_alnset_tree)
from ORBL_tools.orblq import (
    ORBLqCalculator, BiotypesWithFS, get_untrans_ORBLv_dict_name, get_supported_alnsets)

VersionStr = 'Beta'

MaxMaxCodons = 40000

print_err = lambda *pArgs : print(*pArgs, file = sys.stderr)

UsageStr = ('(ALIGNMENT_SET (e.g., hg38_120mammals_primate) [--orblq] [--components] '
            '[FILE] |\n        -h,--help | -v,--version)')
AlnsetsURL = 'https://data.broadinstitute.org/compbio1/cav.php?Alnsets'

def main() :
    ## Parse arguments; print user guide if requested
    if (check_arg('-h', remove = True, silent = True) or
        check_arg('--help', remove = True, silent = True)) :
        assert_num_args(0, UsageStr, exact = True)
        # Display README.md, converting some of the markdown to terminal escape codes.
        with open(os.path.join(ThisDir, 'README.md'), 'rt') as readmeFile :
            outStr = poor_mans_markdown_interpreter(readmeFile.read())
            interactive = sys.stdout.isatty()
            # If interactive, use "less -R" (the -R is Needed to handle escape sequences)
            # otherwise use "cat" and write output to sys.stdout for use in unit test
            # (note that the subprocess ignores programmatic redirection of sys.stdout)
            pageProcess = subprocess.Popen(
                ['less', '-R'] if interactive else ['cat'],
                stdin = subprocess.PIPE,
                stdout = None if interactive else subprocess.PIPE)
            try :
                output, _ = pageProcess.communicate(
                    input = outStr.encode('utf-8')) # Convert str to bytes
            except KeyboardInterrupt : # pragma: no cover
                pass
            if not interactive :
                sys.stdout.write(output.decode('utf-8')) # stdout expects str not bytes
        return
    if (check_arg('-v', remove = True, silent = True) or
        check_arg('--version', remove = True, silent = True)) :
        assert_num_args(0, UsageStr, exact = True)
        print_err('Version: %s' % VersionStr)
        return
    calcOrblQ = check_arg('--orblq', remove = True, silent = True)
    outputComponents = check_arg('--components', remove = True, silent = True)
    assert_num_args(1, UsageStr, maxNumAllowed = 2)
    if any(arg.startswith('-') for arg in sys.argv[1:]) : # Non-existent option
        assert_num_args(1, UsageStr, maxNumAllowed = 2, forceFailure = True)
    alnset = sys.argv[1]
    inFileName = None if len(sys.argv) == 2 else sys.argv[2]

    # Suppress IDE warnings about possible uninitialized variables
    orblqCalc = biotypeWithFS = None

    ## Get ORBLq calculator
    if calcOrblQ :
        untransORBLvFileName = get_untrans_ORBLv_dict_name(alnset)
        if not os.path.exists(untransORBLvFileName) :
            print_err('ERROR: ORBLq is not implemented for alignment set "%s".\n' % alnset
                      + 'Supported alignment sets: ' + ', '.join(get_supported_alnsets()))
            raise SystemExit(1)
        orblqCalc = ORBLqCalculator(alnset)

    ## Download tree file; check for valid alnset
    try :
        tree = download_alnset_tree(alnset)
    except NotImplementedError :
        print_err('ERROR: Could not download tree file for "%s".\n' % alnset
                  + 'Choose alignment set from list at: %s' % AlnsetsURL)
        raise SystemExit(1)
    tree.set_descendants()  # Improves speed of repeated calls to subtree_branch_length.
    fullBL = tree.subtree_branch_length()

    ## For interactive input, print instruction message
    interactive = inFileName is None and sys.stdin.isatty()
    if interactive :
        if calcOrblQ :
            print('Enter: intervals strand biotype+frame_shift (tab-separated) '
                  '(e.g., chr10:1042727-1042762+chr10:1043301-1043321   -   intORF+1)',
                  file = sys.stderr)
        else :
            print('Enter: intervals strand (tab-separated) '
                  '(e.g., chr10:1042727-1042762+chr10:1043301-1043321   -)',
                  file = sys.stderr)

    ## Iterate through input lines and print output lines
    inFile = sys.stdin if inFileName is None else open(inFileName, 'rt')
    try :
        for lineNumber in itertools.count() : # "for line in ..." waits for EOF to start
            ## Parse input line
            inputLine = inFile.readline()
            if inputLine == '' : # End of file
                break
            inputLine = inputLine.strip('\n')
            if interactive and inputLine == '' : # Empty line terminates when interactive
                break

            def process_line(line) :
                """Process input line and return the corresponding output line."""

                if line.startswith('#') : # Comment line. Output as is.
                    return line

                # It's convenient for command line input to allow ' ' instead of tab
                #   between region, strand, and biotype (and convert to tab for downstream
                #   processing), but don't convert spaces farther down in the input line.
                outputLine = '\t'.join(line.split(None, 2))

                ## Parse and check region and strand
                words = outputLine.split('\t') # ''.split('\t') is [''], so length > 0
                numInputWords = len(words)
                intervalsStr = words[0]
                try :
                    # intervalTriples = [(chrom, start, end), ...]
                    intervalTriples = regionString_to_triples(intervalsStr)
                except ValueError :
                    return outputLine + '\t' + 'InvalidRegion'
                if not all(start <= end for chrom, start, end in intervalTriples) :
                    return outputLine + '\t' + 'BackwardsInterval'
                if numInputWords < 2 :
                    return outputLine + '\t' + 'MissingStrand'
                strand = words[1]
                if strand not in ['+', '-'] :
                    return outputLine + '\t' + 'InvalidStrand'

                ## Download alignment
                try :
                    aliSeg = download_local_alignment(intervalsStr, strand, alnset,
                                                      maxCodons = MaxMaxCodons)
                except NotImplementedError :
                    return outputLine + '\t' + 'UnableToDownloadAlignment'

                ## Calculate ORBLv and add it to the output line
                relBLs = calc_ORF_relBLs(aliSeg, tree, fullBL = fullBL)
                orblv = relBLs[3]
                outputLine += '\t%s' % orblv

                ## Calculate ORBLq and add it to the output line, if requested
                if calcOrblQ :
                    biotypeWithFS = words[2] if numInputWords >= 3 else None
                    if numInputWords < 3 :
                        orblq = 'MissingBiotypeWithFS'
                    elif biotypeWithFS not in BiotypesWithFS + ['mixed'] :
                        orblq = 'InvalidBiotypeWithFS'
                    elif biotypeWithFS == 'mixed' or orblv == 'NA' :
                        orblq = 'NA'
                    else :
                        numCodonsWithStop = get_intervals_length(intervalTriples) // 3
                        orblq = orblqCalc(biotypeWithFS, orblv, numCodonsWithStop)
                    outputLine += '\t%s' % orblq

                ## Add ORBLv components to the output line, if requested
                if outputComponents :
                    atgRelBL, stopRelBL, frameRelBL = relBLs[:3]
                    outputLine += '\t%s\t%s\t%s' % (atgRelBL, stopRelBL, frameRelBL)

                return outputLine

            ## Print output line to stdout
            print(process_line(inputLine))
            sys.stdout.flush()
    finally :
        if inFileName is not None :
            inFile.close()

"""
Potential future arguments: 
    [--nTermExt] [--secondORF] [--allowAltStart] [--allowAllNearCognate]
"""

def poor_mans_markdown_interpreter(inStr) :
    """
    Take a string in Markdown format and adjust it to something that can be printed on
        a unix terminal.
    Handles only an extremely limited subset of Markdown syntax:
    - Converts headings (anything starting with '# ', '## ', etc) to bold text. For
      level 1 heading, bold and underlined.
    - Indents code blocks three spaces (but only recognizes them if
      the starting and ending ``` are on lines with nothing else).
    This is a pure python implementation. If that isn't required, use "rich" or
        "markdown2" plus "html2text" for better Markdown interpretation.
    """
    escSeq = '\033['
    inLines = inStr.split('\n')
    outStr = ''
    inCode = False
    for line in inLines :
        if line.strip() == '```' :
            inCode = not inCode
        elif inCode :
            outStr += ' ' * 3 + line
        elif line.startswith('#') and line.lstrip('#').startswith(' ') :
            formatStr = '1;4m' if line.startswith('# ') else '1m'
            outStr += escSeq + formatStr + line.lstrip('#')[1:] + escSeq + '0m'
        else :
            outStr += line
        outStr += '\n'
    return outStr

if __name__ == '__main__' : # pragma: no cover
    main()
