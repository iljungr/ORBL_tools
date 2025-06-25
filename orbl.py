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
import sys, os, itertools, subprocess
ThisDir = os.path.dirname(os.path.realpath(__file__))
sys.path.append(os.path.dirname(ThisDir)) # ORBL_tools's parent
from ORBL_tools.IntervalUtils import regionString_to_triples
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
        # Use "less" to display README.md, converting some of the markdown to terminal
        #   escape codes.
        with open(os.path.join(ThisDir, 'README.md'), 'rt') as readmeFile :
            outStr = poor_mans_markdown_interpreter(readmeFile.read())
            less = subprocess.Popen(['less', '-R'], stdin = subprocess.PIPE)
            try :
                less.communicate(input = outStr.encode('utf-8'))
            except KeyboardInterrupt :
                pass
        return
    if (check_arg('-v', remove = True, silent = True) or
        check_arg('--version', remove = True, silent = True)) :
        assert_num_args(0, UsageStr, exact = True)
        print_err('Version: %s' % VersionStr)
        return
    calcOrblQ = check_arg('--orblq', remove = True, silent = True)
    outputComponents = check_arg('--components', remove = True, silent = True)
    assert_num_args(1, UsageStr, maxNumAllowed = 2)
    alnset = sys.argv[1]
    if alnset.startswith('-') :
        print_err(UsageStr)
        raise SystemExit(1)
    inFileName = None if len(sys.argv) == 2 else sys.argv[2]

    # Suppress IDE warnings about possible uninitialized variables
    orblqCalc = biotypeWithFS = None

    ## Get ORBLq calculator
    if calcOrblQ :
        untransORBLvFileName = get_untrans_ORBLv_dict_name(alnset)
        if not os.path.exists(untransORBLvFileName) :
            print_err('\nERROR: ORBLq is not implemented for alignment set "%s".\n' % alnset
                      + 'Supported alignment sets: ' + ', '.join(get_supported_alnsets()))
            raise SystemExit(1)
        orblqCalc = ORBLqCalculator(alnset)

    ## Download tree file; check for valid alnset
    try :
        tree = download_alnset_tree(alnset)
    except NotImplementedError :
        print_err('\nERROR: Invalid alignment set or missing tree file for "%s".\n' % alnset
                  + 'Choose from list at: %s' % AlnsetsURL)
        raise SystemExit(1)
    tree.set_descendants()  # Improves speed of repeated calls to subtree_branch_length.
    fullBL = tree.subtree_branch_length()

    ## For interactive input, print instruction message
    if inFileName is None and sys.stdin.isatty() :
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
            inputLine = inputLine.strip('\n')
            if inputLine == '' :
                break
            words = inputLine.split('\t')
            # It's convenient for command line input to allow ' ' instead of tab
            #     between region, strand, and biotype (but convert to tab for downstream).
            if len(words[0].split()) == 2 or (calcOrblQ and len(words[0].split()) == 3) :
                if len(words[0].split()) == 2 :
                    words = words[0].split(None, 1) + words[1:]
                else :
                    words = words[0].split(None, 2) + words[1:]
            elif calcOrblQ and len(words) > 1 and len(words[1].split()) == 2 :
                words = words[:1] + words[1].split() + words[2:]
            inputLine = '\t'.join(words)
            intervalsStr = words[0]
            if len(words) < 2 :
                print_err('\nERROR: Line %d missing region or strand: "%s"' % (
                    lineNumber + 1, inputLine))
                raise SystemExit(1)
            if calcOrblQ and len(words) < 3 :
                print_err('\nERROR: Line %d missing region, strand, or biotype: "%s"' % (
                    lineNumber + 1, inputLine))
                raise SystemExit(1)
            strand = words[1]
            if strand not in ['+', '-'] :
                print_err('\nERROR: Line %d invalid strand: "%s".' % (lineNumber + 1,
                                                                      strand))
                raise SystemExit(1)
            if calcOrblQ :
                biotypeWithFS = words[2]
                if biotypeWithFS not in BiotypesWithFS + ['mixed'] :
                    print_err('\nERROR: Line %d invalid biotype: "%s".' % (
                        lineNumber + 1, biotypeWithFS))
                    print_err('Valid biotypes: %s' % ','.join(BiotypesWithFS + ['mixed']))
                    raise SystemExit(1)

            # Download alignment
            try :
                aliSeg = download_local_alignment(intervalsStr, strand, alnset,
                                                  maxCodons = MaxMaxCodons)
            except NotImplementedError as ex :
                msg = ex.args[0]
                print_err('\nERROR: Unable to download alignment for line %d: "%s"\n%s' %
                        (lineNumber + 1, inputLine, msg))
                raise SystemExit(1)

            outputLine = inputLine

            ## Calculate ORBLv and add it to the output line
            relBLs = calc_ORF_relBLs(aliSeg, tree, fullBL = fullBL)
            orblv = relBLs[3]
            outputLine += '\t%s' % orblv

            ## Calculate ORBLq and add it to the output line, if requested
            if calcOrblQ :
                if biotypeWithFS == 'mixed' or orblv == 'NA' :
                    orblq = 'NA'
                else :
                    numCodonsWithStop = sum(
                        end - start + 1
                        for chrom, start, end in regionString_to_triples(intervalsStr)) // 3
                    orblq = orblqCalc(biotypeWithFS, orblv, numCodonsWithStop)
                outputLine += '\t%s' % orblq

            ## Add ORBLv components to the output line, if requested
            if outputComponents :
                atgRelBL, stopRelBL, frameRelBL = relBLs[:3]
                outputLine += '\t%s\t%s\t%s' % (atgRelBL, stopRelBL, frameRelBL)

            ## Print output line to stdout
            print(outputLine)
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

if __name__ == '__main__' :
    main()
