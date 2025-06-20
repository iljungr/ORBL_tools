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
import sys, os, itertools
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

UsageStr = ('(ALIGNMENT_SET [--orblq] [--components] [FILE_NAME] '
                '| -h,--help | -v,--version)')
AlnsetsURL = 'https://data.broadinstitute.org/compbio1/cav.php?Alnsets'

def main() :
    ## Parse arguments; print user guide if requested
    if (check_arg('-h', remove = True, silent = True) or
        check_arg('--help', remove = True, silent = True)) :
        assert_num_args(0, UsageStr, exact = True)
        print_err(UserGuide)
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

UserGuide = """
Overview:

ORBL is a tool for measuring evolutionary conservation and constraint on the "ORFness" 
of an open reading frame (ORF) in a clade using multi-species whole genome alignments.
It consists of two scores:

- ORBLv measures conservation of the ORF by calculating the relative branch length of 
  the phylogenetic tree of species in an alignment that have an intact orthologous ORF, 
  i.e., in which there is an aligned start codon, stop codon, and open reading frame. 
  It is a number between 0 and 1, with larger numbers indicating more conservation. 

- ORBLq measures evolutionary constraint on the ORFness of a non-canonical ORF (ncORF)
  by calculating the quantile of its ORBLv score among the ORBLv scores of untranslated 
  ORFs of the same biotype and similar length. It too is a number between 0 and 1, with 
  larger numbers indicating more constraint.

Usage Summary: 

ORBL_tools/orbl.py %s

Mandatory arguments:
    ALIGNMENT_SET Specifies the multispecies alignment as defined by CodAlignView.

Optional arguments:
    -h, --help     Print this message and exit.
    -v, --version  Print the orbl.py version number and exit.
    --orblq        Report the ORBLq constraint score as well as the orblv conservation 
                   score. This requires specifying the biotype-with-frameshift as the 
                   third field on each input line. Currently, --orblq is only implemented 
                   for certain alignment sets:
    --components   Report the relative branch lengths of the species satisfying each
                   or the three conditions used to consider the ORF to be conserved when
                   calculating the ORBLv score, namely having an aligned start codon, 
                   an aligned stop codon, and having an open reading frame (multiple of
                   3 nucleotides and no in-frame stop codons).
    FILE_NAME      Read input lines from FILE_NAME instead of standard input.


Details:

orbl.py takes input from the standard input or a specified file. Input consists of one 
or more lines, each representing an ORF in the reference species of a multispecies whole 
genome alignment, specified by the ALIGNMENT_SET mandory argument. Alignment sets are
defined by CodAlignView here: %s. 

Each line contains two or more tab-separates fields. The first field is one or nore
chromosomal intervals specifying the coordinates of an open reading frame in the 
reference species of the alignment. The second field is the strand, either + or -. 
The format of the first field consists of one or more intervals separated by plus signs: 
   chrom:start1-end1+chrom:start2-end2+... 
satisfying the requirement that:
   start1 <= end1 < start2 <= end2... 
(even if the region is on the minus strand). 
All segments must be on the same chromosome.

Example input line:
chr10:1042727-1042762+chr10:1043301-1043321 +

Coordinates follow the convention of GFF/GTF files and the UCSC genome browser, namely, 
intervals go from the first base of the region to the last base, inclusive, and 
coordinates are 1-based, i.e., the 1st base of the chromosome is position 1. Coordinates
may include commas, which are ignored.

If the region specified does not begin with ATG or end in a stop codon in the reference 
species then both ORBLv and ORBLq will be "NA".

If --orblq is specified, each input line must inlude a third field, containing the biotype
and frameshift relative to the main ORF for biotypes uoORF, intORF, and doORF. Frameshift 
is +1 or +2 depending on whether a ribosome reading the main frame would need to skip 
1 or 2 nucleotides to get in the frame of the ncORF. Valid values of
biotype-with-frameshift are:
   uORF, uoORF+1, uoORF+2, intORF+1, intORF+2, doORF+1, doORF+2, dORF, and lncRNA-ORF.
Biotype "mixed" is also allowed, but ORBLq is not implemented for this biotype and 
results in a value of "NA".

Currently, --orblq is only implemented for the following alignment sets:
    %s

Input lines may contain additional tab-separated fields, which are passed through to 
the output and otherwise ignored. 

Results are written to the standard output. There is one output line for each input line. 
Each output line contains the input line followed by the orblv score, the orblq score 
if the --orblq option is specified, and three additional fields if --components is 
specified, namely the relative branch lengths of species having an aligned start codon, 
an aligned stop codon, and an intact open reading frame.

Credits:

Questions should be directed to Irwin Jungreis at iljungr@csail.mit.edu.

Citing ORBL:

If you use ORBL, please cite, "High-quality peptide evidence for annotating 
non-canonical open reading frames as human proteins" by Deutsch et al 
(manuscript submitted). 

More information about ORBL can be found in that paper. The specific requirements 
for determining an ORF's biotype are reported in the supplementary materials.
""" % (UsageStr, AlnsetsURL, ', '.join(get_supported_alnsets()))

"""
Potential future arguments: 
    [--nTermExt] [--secondORF] [--allowAltStart] [--allowAllNearCognate]
"""

if __name__ == '__main__' :
    main()
