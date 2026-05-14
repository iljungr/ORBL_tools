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
import sys, os, subprocess, re
ThisDir = os.path.dirname(os.path.realpath(__file__))
sys.path.append(os.path.dirname(ThisDir)) # ORBL_tools's parent
from ORBL_tools.IntervalUtils import regionString_to_triples, get_intervals_length
from ORBL_tools.CmdLineUtils import check_arg, assert_num_args, get_associated_arg
from ORBL_tools.CalcORFrelBLs import calc_ORF_relBLs
from ORBL_tools.DownloadLocalAlignment import (
    download_local_alignment, download_alnset_tree)
from ORBL_tools.orblq import (
    ORBLqCalculator, BiotypesWithFS, get_untrans_ORBLv_dict_name)
from ORBL_tools.ShowAlignmentSets import (
    show_alignment_sets, get_alignment_set_assembly)

VersionStr = '1.0'

print_err = lambda *pArgs : print(*pArgs, file = sys.stderr)

UsageStr = (
    '(ALIGNMENT_SET (e.g., hg38_58) [INPUT_FILE [OUTPUT_FILE]] [--orblq] [--components]\n'
    '         [--CodAlignView] [--UCSCView] [--header] [--extraFields FIELD1,...]\n'
    '         [--writeBed BED_FILE [--ORFNameField FIELD_NAME]] |\n'
    '         -h,--help | -v,--version | --alignmentSets)'
)
CodAlignViewBaseURL = 'https://data.broadinstitute.org/compbio1/cav.php'
UCSCViewBaseURL = 'https://genome.ucsc.edu/cgi-bin/hgTracks'

def main() :
    """
    Parse and check arguments.
    If requested, print help, version, or list of Alignment Sets.
    Otherwise call calc_orbls to calculate ORBL and print requested information for
        input ORFs.
    """
    ## If requested, print user guide and exit.
    if (check_arg('-h', remove = True, silent = True) or
        check_arg('--help', remove = True, silent = True)) :
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

    ## If requested, print version and exit.
    if (check_arg('-v', remove = True, silent = True) or
        check_arg('--version', remove = True, silent = True)) :
        print_err('Version: %s' % VersionStr)
        return

    ## If requested, open browser with Alignment Set information and exit.
    if check_arg('--alignmentSets', remove = True, silent = True) :
        show_alignment_sets()
        return

    ## Process and check arguments
    calcOrblQ        = check_arg('--orblq',        remove = True, silent = True)
    outputComponents = check_arg('--components',   remove = True, silent = True)
    writeCAV         = check_arg('--CodAlignView', remove = True, silent = True)
    writeUCSC        = check_arg('--UCSCView',     remove = True, silent = True)
    writeHeader      = check_arg('--header',       remove = True, silent = True)
    extraFieldStr = get_associated_arg('--extraFields', None, remove = True, silent = True)
    extraFields = [] if extraFieldStr is None else extraFieldStr.split(',')
    bedFileName  = get_associated_arg('--writeBed', None, remove = True, silent = True)
    orfNameField = get_associated_arg('--ORFNameField', None, remove = True, silent = True)
    if orfNameField is not None :
        if bedFileName is None :
            print_err('--ORFNameField is only allowed with --writeBed.')
            raise SystemExit(1)
        if ((extraFields and orfNameField not in extraFields) or
            (not extraFields and re.match(r'^Extra[1-9][0-9]*$', orfNameField) is None)) :
                print_err('No field named %s.' % orfNameField)
                raise SystemExit(1)
    assert_num_args(1, UsageStr, maxNumAllowed = 3)
    if any(arg.startswith('-') for arg in sys.argv[1:]) : # Non-existent option
        assert_num_args(1, UsageStr, maxNumAllowed = 2, forceFailure = True)
    alnset = sys.argv[1]
    inFileName  = None if len(sys.argv) < 3 else sys.argv[2]
    outFileName = None if len(sys.argv) < 4 else sys.argv[3]

    ## Calculate ORBL for input ORFs.
    calc_orbls(alnset, inFileName, outFileName, calcOrblQ, outputComponents,
               writeCAV, writeUCSC, writeHeader, extraFields, bedFileName, orfNameField)

def calc_orbls(
        alnset, inFileName = None, outFileName = None,
        calcOrblQ = False, outputComponents = False,
        writeCAV = False, writeUCSC = False, writeHeader = False,
        extraFields = (), bedFileName = None, orfNameField = None) :
    """
    Calculate ORBLv and, if specified, ORBLq for each non-comment line in inFileName
        (or stdin if inFileName is None), writing resulta and other requestws information
        to outFileName (or stdout if that is None).
    """

    ## Get ORBLq calculator
    if calcOrblQ :
        untransORBLvFileName = get_untrans_ORBLv_dict_name(alnset)
        if not os.path.exists(untransORBLvFileName) :
            print_err('ERROR: ORBLq is not implemented for alignment set "%s".\n' % alnset
                      + 'Use --alignmentSets argument to see allowed alignment sets.')
            raise SystemExit(1)
        orblqCalc = ORBLqCalculator(alnset)
    else :
        orblqCalc = None # Suppress IDE warnings about possible uninitialized variables

    ## Download tree file; check for valid alnset
    try :
        tree = download_alnset_tree(alnset)
    except NotImplementedError :
        print_err('ERROR: Could not download tree file for "%s".\n' % alnset
                  + 'Use --alignmentSets argument to see allowed alignment sets.')
        raise SystemExit(1)
    tree.set_descendants()  # Improves speed of repeated calls to subtree_branch_length.
    fullBL = tree.subtree_branch_length()

    ## For interactive input, print instruction message
    interactive = inFileName is None and sys.stdin.isatty()
    if interactive :
        fieldsToEnter = ['Intervals', 'Strand'] + list(extraFields)
        exampleValues = (['chr10:1042727-1042762+chr10:1043301-1043321', '-'] +
                         ['value%d' % ii for ii in range(1, len(extraFields) + 1)])
        if calcOrblQ :
            fieldsToEnter.insert(2, 'Biotype+FrameShift')
            exampleValues.insert(2, 'intORF+1')
        print('Enter: %s (tab-separated)\n e.g., %s' %
              ('\t'.join(fieldsToEnter), '\t'.join(exampleValues)), file = sys.stderr)

    # These will be changed below if the corresponding file name is not None
    inFile  = sys.stdin
    outFile = sys.stdout
    bedFile = None

    try :
        ## Open output bed file if requested
        if inFileName is not None :
            inFile = open(inFileName, 'rt')
        if outFileName is not None :
            outFile = open(outFileName, 'wt')
        if bedFileName is not None :
            bedFile = open(bedFileName, 'wt')

        ## Write header if requested and if we already know extra field names
        if writeHeader and extraFields :
            print(_make_header(calcOrblQ, extraFields, outputComponents,
                               writeCAV, writeUCSC), file = outFile)
            outFile.flush()

        ## Iterate through input lines and print output lines
        extraFieldsUnset = not extraFields
        orfNameCounter = {}  # (chrom, strand, firstCoord) : countSoFar
        try :
            assembly = get_alignment_set_assembly(alnset) if writeUCSC else None
        except ValueError : # pragma: no cover
            # Should never happen; would have already failed to get tree.
            print_err('ERROR: Could not download information for "%s".\n' % alnset
                      + 'Use --alignmentSets argument to see allowed alignment sets.')
            raise SystemExit(1)
        while True : # Can't use "for line in ..." because it waits for EOF to start
            ## Process one input line
            inputLine = inFile.readline()
            if inputLine == '' : # End of file
                break
            inputLine = inputLine.strip('\n')
            if interactive and inputLine == '' : # Empty line terminates when interactive
                break

            ## If this is the first non-commment line and extraFields was not set,
            ## set extraFields and print header if requested.
            if extraFieldsUnset and not inputLine.startswith('#') and inputLine != '' :
                # Count fields, treating spaces among region, strand, and next field as
                #    field separators, but only tabs among later fields
                numFields = _fix_spaces(inputLine).count('\t') + 1
                numExtra = numFields - (3 if calcOrblQ else 2)
                extraFields = ['Extra%d' % nn for nn in range(1, numExtra + 1)]
                if orfNameField is not None and orfNameField not in extraFields :
                    # Only happens for --ORFNameField ExtraN where N > numExtra.
                    # Obscure case not worth extra call to download_alnset_tree in
                    # the unit test, so "pragma: no cover"
                    print_err('No field named %s.' % orfNameField) # pragma: no cover
                    raise SystemExit(1) # pragma: no cover
                extraFieldsUnset = False
                if writeHeader :
                    print(_make_header(calcOrblQ, extraFields, outputComponents,
                                       writeCAV, writeUCSC), file = outFile)
                    outFile.flush()

            ## Calculate and print corresponding output line
            outLine = _process_line(
                inputLine, interactive, bedFile, orfNameField, extraFields, calcOrblQ,
                orfNameCounter, alnset, tree, fullBL, orblqCalc, outputComponents,
                writeCAV, writeUCSC, assembly)
            print(outLine, file = outFile)
            outFile.flush()
    finally :
        # Check inFile, etc., not inFileName, etc., in case of exception when not all have
        # been opened.
        if inFile is not sys.stdin :
            inFile.close()
        if outFile is not sys.stdout :
            outFile.close()
        if bedFile is not None :
            bedFile.close()

"""
Potential future arguments: 
    [--nTermExt] [--secondORF] [--allowAltStart] [--allowAllNearCognate]
"""

def _fix_spaces(line) :
    """
    It's convenient for command line input to allow ' ' instead of tab between region,
        strand, and biotype, and to convert these to tabs for downstream processing,
        but don't convert spaces farther down in the input line to allow extra fields to
        contain spaces.
    """
    words = line.split('\t')
    words0 = words[0].split(None, 2)
    words[:1] = words0
    if len(words) > 1 : # Else bad input, but don't crash here
        words1 = words[1].split(None, 1)
        words[1:2] = words1
    return '\t'.join(words)

def _make_header(calcOrblQ, extraFields, outputComponents, writeCAV, writeUCSC) :
    """ Return a tab-separated header line for the output. """
    fieldNames = ['Intervals', 'Strand']
    if calcOrblQ :
        fieldNames.append('BiotypeWithFS')
    fieldNames.extend(extraFields)
    fieldNames.append('ORBLv')
    if calcOrblQ :
        fieldNames.append('ORBLq')
    if outputComponents :
        fieldNames.extend(['ATGRelBL', 'StopRelBL', 'FrameRelBL'])
    if writeCAV :
        fieldNames.append('CodAlignView')
    if writeUCSC :
        fieldNames.append('UCSCView')
    return '\t'.join(fieldNames)

def _process_line(line, interactive, bedFile, orfNameField, extraFields, calcOrblQ,
                  orfNameCounter, alnset, tree, fullBL, orblqCalc, outputComponents,
                  writeCAV, writeUCSC, assembly) :
    """Process input line and return the corresponding output line."""
    if line.startswith('#') :  # Comment line. Output as is.
        return line
    if line == '' :  # Blank line.
        assert not interactive  # Otherwise we should have already terminated
        return line  # Blank line. Output as is.
    outputLine = _fix_spaces(line)

    ## Parse and check region and strand
    words = outputLine.split('\t')
    numInputWords = len(words)
    intervalsStr = words[0]
    try :
        # intervalTriples = [(chrom, start, end), ...]
        intervalTriples = regionString_to_triples(intervalsStr)
    except ValueError :
        return outputLine + '\t' + 'InvalidRegion'
    chrom = intervalTriples[0][0]
    if any(c != chrom for c, start, end in intervalTriples) :
        return outputLine + '\t' + 'MoreThanOneChromosome'
    if not all(start <= end for c, start, end in intervalTriples) :
        return outputLine + '\t' + 'BackwardsInterval'
    if numInputWords < 2 :
        return outputLine + '\t' + 'MissingStrand'
    strand = words[1]
    if strand not in ['+', '-'] :
        return outputLine + '\t' + 'InvalidStrand'

    ## Print bed line if requested
    if bedFile is not None :
        intervals = [triple[1 :] for triple in intervalTriples]
        if orfNameField is not None :
            orfNameIndex = (extraFields.index(orfNameField) +
                            (3 if calcOrblQ else 2))
            name = outputLine.split('\t')[orfNameIndex]
        else :
            firstCoord = intervals[0][0]
            key = (chrom, strand, firstCoord)
            index = orfNameCounter.setdefault(key, 0) + 1
            name = '%s_%s_%d_%d' % (chrom, strand, firstCoord, index)
            orfNameCounter[key] += 1
        print(intervals_to_bed_line(chrom, intervals, strand, name),
              file = bedFile)

    # Check correct number of extra fields
    numExtraThisLine = numInputWords - (3 if calcOrblQ else 2)
    if numExtraThisLine > len(extraFields) :
        return outputLine + '\t' + 'TooManyExtraFields'
    if 0 <= numExtraThisLine < len(extraFields) :
        # Note: the only way we can get here with numExtraThisLine < 0, is if
        #    calcOrblQ and missing biotype. We'd rather report that than
        #    wrong number of extra fields because it's a more serious error;
        #    we catch and report it later.
        return outputLine + '\t' + 'TooFewExtraFields'

    ## Download alignment
    try :
        aliSeg = download_local_alignment(intervalsStr, strand, alnset)
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

    ## Add CodAlignView if requested
    if writeCAV :
        # Add 15 nt prologue and epilogue. Set one line per species (&os).
        CAVurl = CodAlignViewBaseURL + '?a=%s&i=%s&s=%s&p=15&e=15&os' % (
            alnset, intervalsStr, strand)
        outputLine += '\t' + CAVurl

    ## Add UCSCView if requested
    if writeUCSC :
        fromPos = intervalTriples[0][1]
        toPos = intervalTriples[-1][2]
        reverseStr = '&complement_%s=%d&hgt.revCmplDisp_%s=%d' % (
            assembly, strand == '-', assembly, strand == '-')
        UCSCurl = UCSCViewBaseURL + '?db=%s&position=%s:%d-%d%s' % (
            assembly, chrom, fromPos, toPos, reverseStr)
        outputLine += '\t' + UCSCurl

    return outputLine

def intervals_to_bed_line(chrom, intervals, strand, name = '') :
    """ Return a .bed formatted line with specified data. """
    chromStart = thickStart = intervals[0][0] - 1 # Bed counts from 0 instead of 1.
    chromEnd   = thickEnd   = intervals[-1][1]    # But chromEnd is position _after_ end.
    return '\t'.join(map(str, [
        chrom, chromStart, chromEnd, name, 0, strand,
        thickStart, thickEnd, 0, len(intervals),
        ','.join(map(str, [end - start + 1 for start, end in intervals])),
        ','.join(map(str, [start - chromStart - 1 for start, end in intervals])),
        ]))

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
