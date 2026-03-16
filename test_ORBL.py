#!/usr/bin/env python
# Copyright 2024 Irwin Jungreis
# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at http://www.apache.org/licenses/LICENSE-2.0
# Unless required by applicable law or agreed to in writing, software
# distributed under the License is distributed on an "AS IS" BASIS,
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
# See the License for the specific language governing permissions and
# limitations under the License.
"""
Test orbl.py and orblq.py.
"""
from __future__ import division, print_function
from __future__ import absolute_import
import sys, os, unittest, tempfile, subprocess
UnitTestDir = os.path.dirname(os.path.realpath(__file__))
sys.path.append(os.path.dirname(UnitTestDir))# ORBL_tools's parent
from ORBL_tools.orblq import ORBLqCalculator
import ORBL_tools.orbl
from ORBL_tools.orbl import main as orblMain
from ORBL_tools.orbl import UsageStr as orblUsageStr
from ORBL_tools.CalcORFrelBLs import _test_calc_seq_conservations
from ORBL_tools.ContextMgrsForTesting import (
    StreamCatcher, SysArgSaver, StdinSaver, InteractiveStringIO)
import ORBL_tools.ShowAlignmentSets as ShowAlignmentSets
if sys.version_info[0] < 3 :
    from itertools import izip_longest as zip_longest
else :
    from itertools import zip_longest as zip_longest

PrimateAlnset   = 'hg38_120mammals_primate'
PlacentalAlnset = 'hg38_120mammals_placental'

class TestORBLq(unittest.TestCase) :
    def test_ORBLqCalculator(self) :
        at = self.assertTrue
        ar = self.assertRaises
        aaeq = lambda x, y : at(abs(x - y) < 1e-11)  # "assertAlmostEqual"

        # Test for primate
        primateCalc = ORBLqCalculator(PrimateAlnset)
        primORBLv = 0.0139718396632
        primORBLq = 0.254531126872
        aaeq(primateCalc('lncRNA-ORF', primORBLv, 92), primORBLq)

        # Test for placental. Test matchedORBLvs.
        placentalCalc = ORBLqCalculator(PlacentalAlnset)
        placORBLv = 0.0642060273392
        placORBLq = 0.657398212512
        matchedORBLvs = []
        aaeq(placentalCalc('uoORF+2', placORBLv, 100, matchedORBLvs), placORBLq)
        at(len(matchedORBLvs) >= 1000)
        at(all(type(ov) == float for ov in matchedORBLvs))
        at(all(0 <= ov <= 1 for ov in matchedORBLvs))
        aaeq(1 - placORBLq,
             (sum(v >= placORBLv - 1e-11 for v in matchedORBLvs) + 1) / len(matchedORBLvs))

        # Invalid biotype
        ar(AssertionError, primateCalc, 'uoORF+3', primORBLv, 92)

        # ORBLq not implemented for mixed biotype
        ar(AssertionError, primateCalc, 'mixed',   primORBLv, 92)

class TestORBL(unittest.TestCase) :
    def test_orbl(self) :
        at = self.assertTrue
        ae = self.assertEqual
        ar = self.assertRaises
        ain = self.assertIn
        aaeq = lambda x, y : at(abs(x - y) < 1e-11)  # "assertAlmostEqual"
        currentVersion = 'Beta'

        # Get names for temporary files
        tempHandle, tempInName  = tempfile.mkstemp(suffix = '.in')
        os.close(tempHandle)
        tempHandle, tempOutName = tempfile.mkstemp(suffix = '.tsv')
        os.close(tempHandle)
        tempHandle, tempBedName = tempfile.mkstemp(suffix = '.bed')
        os.close(tempHandle)

        def compare_output_lines(l1, l2) :
            """
            Compare two lines allowing tiny differences in floats. This is needed because
                string conversion of floats is different in Python 2 and Python 3, so
                output of orbl will differ.
            """
            words1 = l1.split('\t')
            words2 = l2.split('\t')
            ae(len(words1), len(words2), str((words1, words2)))
            for word1, word2 in zip(words1, words2) :
                if word1 != word2 :
                    try : # Allow it if both are floats and they are almost equal
                        aaeq(float(word1), float(word2))
                    except ValueError :
                        self.fail('%s != %s' % (word1, word2))

        def compare_outputs(t1, t2) :
            ae(t1.count('\n'), t2.count('\n'))
            for l1, l2 in zip(t1.split('\n'), t2.split('\n')) :
                compare_output_lines(l1, l2)

        with StreamCatcher('both') as streamCatcher :
            with SysArgSaver() :
                ## Test help
                for helpArg in ['-h', '--help'] :
                    sys.argv = ['orbl.py', helpArg]
                    orblMain()
                    ae(streamCatcher.buffer('err'), '') # Nothing written to stderr
                    outStr = streamCatcher.buffer('out')
                    # Check approximate length, leaving room for future changes to README.md
                    at(13000 < len(outStr) < 26000, 'outStr length = %d' % len(outStr))
                    ain('constraint', outStr) # Likely regardless of changes to README.md
                    streamCatcher.clear('both')

                ## Test version
                for versionArg in ['-v', '--version'] :
                    sys.argv = ['orbl.py', versionArg]
                    orblMain()
                    ae(streamCatcher.buffer('err'), 'Version: %s\n' % currentVersion)
                    ae(streamCatcher.buffer('out'), '')
                    streamCatcher.clear('both')

                ## Test --alignmentSets
                saveShowAlnsetsAsText = ShowAlignmentSets.SHOW_ALNSETS_AS_TEXT
                ShowAlignmentSets.SHOW_ALNSETS_AS_TEXT = True
                sys.argv = ['orbl.py', '--alignmentSets']
                orblMain()
                ae(streamCatcher.buffer('out'), '')
                outputLines = streamCatcher.buffer('err').strip('\n').split('\n')
                firstLineStart, htmlFileName = outputLines[0][:-1].rsplit(None, 1)
                at(os.path.exists(htmlFileName))
                ae(firstLineStart, 'Alignment Set list was written to')
                ae(outputLines[1],
                   'However, could not open web browser. '
                   'Printing a text version instead.')
                ae(outputLines[2].split(), # Split to allow changes in number of spaces
                   ['Alignment', 'Set', 'ORBLq', 'Reference',
                    'Number', 'of', 'Source', 'Subset'])
                ae(outputLines[3].split(),
                   ['Assembly', 'Species', 'Alignment', 'Description'])
                ae(outputLines[-4], '')
                ae(outputLines[-3], 'Tree files in Newick and pdf format can be found in:')
                ae(outputLines[-2],
                   '    https://data.broadinstitute.org/compbio1/CodAlignViewFiles'
                   '/TreeNHs/{Alignment Set}.nh')
                ae(outputLines[-1],
                   '    https://data.broadinstitute.org/compbio1/CodAlignViewFiles'
                   '/TreePdfs/{Alignment Set}.pdf')
                tableLines = outputLines[4 : -4]
                streamCatcher.clear('both')
                with open(htmlFileName, 'rt') as htmlFile :
                    alnsetsHtmlEd = ShowAlignmentSets.AlnsetsHtmlEditor(htmlFile.read())
                    ae(alnsetsHtmlEd.get_column_names(),
                       ['Alignment Set', 'ORBLq', 'Reference<br>Assembly',
                        'Number of<br>Species', 'Source<br>Alignment',
                        'Subset<br>Description', 'Tree NH', 'Tree PDF'])
                    rowNames = alnsetsHtmlEd.get_row_names()
                    # Check that text table and html table have same Alignment Set names
                    ae([line.split()[0] for line in tableLines], rowNames)
                    # Check that text table and html table have the same ORBLq status
                    for rowName, line in zip(rowNames, tableLines) :
                        ae(alnsetsHtmlEd.get_values_from_row_name(rowName)[1] == 'yes',
                           line.split()[1] == 'yes')
                    # Hidden alnset with ORBLq defined
                    ain('hg38_447mammals_chimp', alnsetsHtmlEd.get_row_names())
                    ae(alnsetsHtmlEd.get_values_from_row_name('hg38_447mammals_chimp')[1],
                       'yes')
                    # Hidden alnset with ORBLq undefined
                    self.assertNotIn('hg38_7', alnsetsHtmlEd.get_row_names())
                    # Visible alnset with ORBLq defined
                    ae(alnsetsHtmlEd.get_values_from_row_name('hg38')[1], 'yes')
                    # Visible alnset with ORBLq undefined
                    ae(alnsetsHtmlEd.get_values_from_row_name('hg18')[1], '')
                    # Non-existent alnset
                    ar(ValueError, alnsetsHtmlEd.get_values_from_row_name, 'hg17')
                ShowAlignmentSets.SHOW_ALNSETS_AS_TEXT = saveShowAlnsetsAsText

                ## Test non-existent option
                sys.argv = ['orbl.py', '--q'] # No such option
                try :
                    orblMain()
                except SystemExit as ex :
                    ae(ex.code, 1)
                    ae(streamCatcher.buffer('err'), 'Usage: ' + orblUsageStr + '\n')
                    ae(streamCatcher.buffer('out'), '')
                else :
                    self.fail('orbl.py --q did not raise')
                streamCatcher.clear('both')

                ## Test non-existent alignment set
                fakeAlignmentSet = 'hg31'
                sys.argv = ['orbl.py', fakeAlignmentSet]
                try :
                    orblMain()
                except SystemExit as ex :
                    ae(ex.code, 1)
                    at(streamCatcher.buffer('err').startswith(
                        'ERROR: Could not download tree file for "%s".\n' %
                        fakeAlignmentSet), streamCatcher.buffer('err'))
                    ae(streamCatcher.buffer('out'), '')
                else :
                    self.fail('orbl.py %s did not raise' % fakeAlignmentSet)
                streamCatcher.clear('both')

                ## Test simple example with no options
                # Also tests placental, ' ' twixt region & strand, read from stdin, commas
                regionStr = 'chr1:3891002-3891057+chr1:3,892,883-3,892,889'
                sys.argv = ['orbl.py', PlacentalAlnset]
                inputLine = regionStr + ' -'
                with StdinSaver(inputLine) :
                    orblMain()
                ae(streamCatcher.buffer('err'), '') # Nothing written to stderr
                compare_output_lines(streamCatcher.buffer('out'),
                                     inputLine.replace(' ', '\t') + '\t0.648645306038\n')
                streamCatcher.clear('both')

                ## Test orblq example. Also tests tabs between region/strand/biotype;
                # extra field; header without extraFields; writeBed without ORFNameField
                inputLine = regionStr + '\t-\tintORF+2\tc1riboseqorf13'
                sys.argv = ['orbl.py', PlacentalAlnset, '--orblq', '--header',
                            '--writeBed', tempBedName]
                with StdinSaver(inputLine) :
                    orblMain() # Warning: removes the --orblq from sys.argv
                ae(streamCatcher.buffer('err'), '') # Nothing written to stderr
                compare_outputs(
                    streamCatcher.buffer('out'),
                    'Intervals\tStrand\tBiotypeWithFS\tExtra1\tORBLv\tORBLq\n' +
                    inputLine + '\t0.648645306038\t0.835130970724\n')
                with open(tempBedName, 'rt') as bedFile :
                    # Tests default ORF name chr1_-_3891002_1. Also 2-exon bed file.
                    compare_outputs(bedFile.read(),
                                    'chr1\t3891001\t3892889\tchr1_-_3891002_1\t0\t-\t'
                                    '3891001\t3892889\t0\t2\t56,7\t0,1881\n')
                streamCatcher.clear('both')

                ## Test real alnset for which orblq is not implemented (hg18)
                noOrblqInputLine = 'chr1:64831034-64831054\t+\tuORF\n'
                sys.argv = ['orbl.py', 'hg18']
                with StdinSaver(noOrblqInputLine) :
                    orblMain() # hg18 works OK without --orblq
                ae(streamCatcher.buffer('err'), '') # Nothing written to stderr
                compare_output_lines(
                    streamCatcher.buffer('out'),
                    'chr1:64831034-64831054\t+\tuORF\t0.00284527142986\n')
                streamCatcher.clear('both')
                sys.argv = ['orbl.py', 'hg18', '--orblq']
                with StdinSaver(noOrblqInputLine) :
                    try :
                        orblMain() # hg18 fails with --orblq
                    except SystemExit as ex :
                        ae(ex.code, 1)
                    else :
                        self.fail('orbl.py hg18 --orblq did not raise')
                ae(streamCatcher.buffer('err'),
                    'ERROR: ORBLq is not implemented for alignment set "hg18".\n'
                    'Use --alignmentSets argument to see allowed alignment sets.\n')
                ae(streamCatcher.buffer('out'), '')
                streamCatcher.clear('both')

                ## Test --components, reading from file, primate alignment, and many cases
                ##     of success and failure for individual input lines.
                # Note: orbl.py infers the number of extra fields (Extra1, Extra2, ...)
                # from the first non-comment input line.  Our test input intentionally
                # begins with a valid line that *does* include an extra field, but many
                # subsequent lines do not.  Split into two runs so the first line does
                # not force extra-field inference for the remaining cases.
                for InputLines, ExpectedOutput in zip(InputLinesList, ExpectedOutputsList) :
                    with open(tempInName, 'wt') as orblInFile :
                        orblInFile.write(InputLines)
                    sys.argv = ['orbl.py',  '--orblq', PrimateAlnset, tempInName,
                                '--components']
                    orblMain()
                    ae(streamCatcher.buffer('err'), '')  # Nothing written to stderr
                    outputLines = streamCatcher.buffer('out').split('\n')
                    for outputLine, expectedLine in zip_longest(
                            outputLines, ExpectedOutput.split('\n')) :
                        compare_output_lines(outputLine, expectedLine)
                    streamCatcher.clear('both')

                ## Test --CodAlignView, --UCSCView, --header, --extraFields, --writeBed,
                ##    --ORFNameField, and output to file.
                with open(tempInName, 'wt') as orblInFile :
                    orblInFile.write(InputForHeaderEtAl)
                sys.argv = [
                    'orbl.py', 'hg38_58', tempInName, tempOutName,
                    '--orblq', '--components', '--CodAlignView', '--UCSCView',
                    '--header', '--extraFields', 'Name,NumCodons',
                    '--writeBed', tempBedName, '--ORFNameField', 'Name']
                orblMain()
                ae(streamCatcher.buffer('err'), '')  # Nothing written to stderr
                ae(streamCatcher.buffer('out'), '')  # Nothing written to stdout
                with open(tempOutName, 'rt') as outFile :
                    compare_outputs(outFile.read(), ExpectedOutputForHeaderEtAl)
                with open(tempBedName, 'rt') as bedFile :
                    ae(bedFile.read(), ExpectedBedForHeaderEtAl)
                streamCatcher.clear('both')

                ## Test some errors
                for args, errStr in [
                    (['--ORFNameField'], 'Argument --ORFNameField requires a value.\n'),
                    ([ '--ORFNameField', 'xxx'],
                     '--ORFNameField is only allowed with --writeBed.\n'),
                    (['--writeBed', tempBedName, '--ORFNameField', 'xxx'],
                     'No field named xxx.\n'),
                ] :
                    sys.argv = ['orbl.py'] + args
                    try :
                        orblMain()
                    except SystemExit as ex :
                        ae(ex.code, 1)
                        ae(streamCatcher.buffer('err'), errStr)
                        ae(streamCatcher.buffer('out'), '')
                    else :
                        self.fail('orbl.py %s did not raise' % ''.join(args))
                    streamCatcher.clear('both')

                ## Test interactive mode
                inputLines = '# First line\nchr1:1\n# comment\n\n# skip this'
                sys.argv = ['orbl.py', PrimateAlnset]
                with StdinSaver(InteractiveStringIO(inputLines)) :
                    orblMain()
                ae(streamCatcher.buffer('err'),
                   'Enter: Intervals\tStrand (tab-separated)\n'
                   ' e.g., chr10:1042727-1042762+chr10:1043301-1043321\t-\n')
                # For interactive intput it should stop at the blank line
                ae(streamCatcher.buffer('out'),
                   '# First line\nchr1:1\tInvalidRegion\n# comment\n')
                streamCatcher.clear('both')
                sys.argv = ['orbl.py', '--orblq', PrimateAlnset]
                with StdinSaver(InteractiveStringIO('# Only line\n')) :
                    orblMain()
                ae(streamCatcher.buffer('err'),
                   'Enter: Intervals\tStrand\tBiotype+FrameShift (tab-separated)\n'
                   ' e.g., chr10:1042727-1042762+chr10:1043301-1043321\t-\tintORF+1\n')
                ae(streamCatcher.buffer('out'), '# Only line\n')
                streamCatcher.clear('both')

                ## Test calling as a command rather than calling main
                commandName = ORBL_tools.orbl.__file__.replace('.pyc', '.py')
                orblProc = subprocess.Popen([commandName, '-v'], stdout = subprocess.PIPE,
                                            stderr = subprocess.PIPE)
                out, error = orblProc.communicate()
                ae(error.decode('utf-8'), 'Version: %s\n' % currentVersion)
                ae(out.decode('utf-8'), '')

        for fileName in [tempInName, tempOutName, tempBedName] :
            if os.path.exists(fileName) :
                os.remove(fileName)


InputLinesList = ['''# Valid two-interval region; spaces between region/strand/biotype
chr1:3891002-3891057+chr1:3892883-3892889 - intORF+2	c1riboseqorf13
chr1:3891002-3891057+chr1:3892883-3892889 - intORF+2
chr1:3891002-3891057+chr1:3892883-3892889 - intORF+2	c1riboseqorf13	21
''',
'''# Invalid region
chr1:3891002-3891055+chr1:3892883	-
# start > end
chr1:3891051-3890992	+
# More than one chromosome
chr1:3891002-3891057+chr2:3892883-3892889	-
# Missing strand
chr1:3890992-3891051
# Invalid strand
chr1:3890992-3891051	.
# No alignment (no such chromosome)
chr25:3890992-3891051	+
# Blank line (don't stop for non-interactive input, just pass through to output)

# Missing biotype
chr1:3890992-3891051	+
# Invalid biotype (space/tab between region/strand/biotype)
chr1:3890992-3891051 +	uORF+1
# Mixed biotype 
chr1:3890992-3891051 +	mixed
# orblv is NA (tab/space between region/strand/biotype)
chr1:3890995-3891051 +	uORF
''']

ExpectedOutputsList = ['''# Valid two-interval region; spaces between region/strand/biotype
chr1:3891002-3891057+chr1:3892883-3892889	-	intORF+2	c1riboseqorf13	\
0.974690114626	0.902157164869	0.974690114626	1.0	0.974690114626
chr1:3891002-3891057+chr1:3892883-3892889	-	intORF+2	TooFewExtraFields
chr1:3891002-3891057+chr1:3892883-3892889	-	intORF+2	c1riboseqorf13	21	\
TooManyExtraFields
''',
'''# Invalid region
chr1:3891002-3891055+chr1:3892883	-	InvalidRegion
# start > end
chr1:3891051-3890992	+	BackwardsInterval
# More than one chromosome
chr1:3891002-3891057+chr2:3892883-3892889	-	MoreThanOneChromosome
# Missing strand
chr1:3890992-3891051	MissingStrand
# Invalid strand
chr1:3890992-3891051	.	InvalidStrand
# No alignment (no such chromosome)
chr25:3890992-3891051	+	UnableToDownloadAlignment
# Blank line (don't stop for non-interactive input, just pass through to output)

# Missing biotype
chr1:3890992-3891051	+	0.432983873668	MissingBiotypeWithFS	0.698355707319	\
0.856181535801	0.432983873668
# Invalid biotype (space/tab between region/strand/biotype)
chr1:3890992-3891051	+	uORF+1	0.432983873668	InvalidBiotypeWithFS	\
0.698355707319	0.856181535801	0.432983873668
# Mixed biotype 
chr1:3890992-3891051	+	mixed	0.432983873668	NA	0.698355707319	0.856181535801	\
0.432983873668
# orblv is NA (tab/space between region/strand/biotype)
chr1:3890995-3891051	+	uORF	NA	NA	NA	0.856181535801	0.432983873668
''']

InputForHeaderEtAl = 'chr1:16395157-16395234	+	dORF	c1riboseqorf29	26'
ExpectedOutputForHeaderEtAl = (
    'Intervals	Strand	BiotypeWithFS	Name	NumCodons	ORBLv	ORBLq	ATGRelBL	'
    'StopRelBL	FrameRelBL	CodAlignView	UCSCView\n'
    'chr1:16395157-16395234	+	dORF	c1riboseqorf29	26	'
    '0.819463199838	0.995604904002	0.871454537388	0.98745916947	0.95279827777	'
    'https://data.broadinstitute.org/compbio1/cav.php?'
    'a=hg38_58&i=chr1:16395157-16395234&s=+&p=15&e=15&os	'
    'https://genome.ucsc.edu/cgi-bin/hgTracks?'
    'db=hg38&position=chr1:16395157-16395234&complement_hg38=0&hgt.revCmplDisp_hg38=0\n')
ExpectedBedForHeaderEtAl = (
    'chr1	16395156	16395234	c1riboseqorf29	0	+	16395156	16395234	'
    '0	1	78	0\n')

class TestCalcSeqConservations(unittest.TestCase) :
    """
    Test the checks for various kinds of conservation used in ORBLv.
    """
    def test_calc_seq_conservations(self) :
        self.assertIs(_test_calc_seq_conservations(), None)

if __name__ == '__main__' :
    unittest.main() # Allows running unit test with "python test_TESTNAME.py"