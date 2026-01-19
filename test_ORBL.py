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

        def compare_output_lines(l1, l2) :
            words1 = l1.split('\t')
            words2 = l2.split('\t')
            ae(len(words1), len(words2), str((words1, words2)))
            for word1, word2 in zip(words1, words2) :
                if word1 != word2 :
                    try : # Allow it if both are floats and they are almost equal
                        aaeq(float(word1), float(word2))
                    except ValueError :
                        self.fail('%s != %s' % (word1, word2))

        with StreamCatcher('both') as streamCatcher :
            with SysArgSaver() :
                ## Test help
                for helpArg in ['-h', '--help'] :
                    sys.argv = ['orbl.py', helpArg]
                    orblMain()
                    ae(streamCatcher.buffer('err'), '') # Nothing written to stderr
                    outStr = streamCatcher.buffer('out')
                    # Check approximate length, leaving room for future changes to README.md
                    at(4000 < len(outStr) < 8000)
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
                savePrintHtmlFileName = ShowAlignmentSets.PRINT_HTML_FILE_NAME
                ShowAlignmentSets.PRINT_HTML_FILE_NAME = True
                sys.argv = ['orbl.py', '--alignmentSets']
                orblMain()
                ae(streamCatcher.buffer('out'), '')
                htmlFileName = streamCatcher.buffer('err')[:-1] # Remove \n
                at(os.path.exists(htmlFileName))
                streamCatcher.clear('both')
                with open(htmlFileName, 'rt') as htmlFile :
                    alnsetsHtmlEd = ShowAlignmentSets.AlnsetsHtmlEditor(htmlFile.read())
                    ae(alnsetsHtmlEd.get_column_names(),
                       ['Alignment Set', 'ORBLq', 'Reference<br>Assembly',
                        'Number of<br>Species', 'Source<br>Alignment',
                        'Subset<br>Description', 'Tree NH', 'Tree PDF'])
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
                ShowAlignmentSets.PRINT_HTML_FILE_NAME = savePrintHtmlFileName

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

                ## Test orblq example
                # Also tests tabs between region/strand/biotype; extra field.
                inputLine = regionStr + '\t-\tintORF+2\tc1riboseqorf13'
                sys.argv = ['orbl.py', PlacentalAlnset, '--orblq']
                with StdinSaver(inputLine) :
                    orblMain() # Warning: removes the --orblq from sys.argv
                ae(streamCatcher.buffer('err'), '') # Nothing written to stderr
                compare_output_lines(streamCatcher.buffer('out'),
                                     inputLine + '\t0.648645306038\t0.835130970724\n')
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
                tempHandle, tempName = tempfile.mkstemp(prefix = 'temp', suffix = '.in')
                os.close(tempHandle)
                with open(tempName, 'wt') as orblInFile :
                    orblInFile.write(InputLines)
                sys.argv = ['orbl.py',  '--orblq', PrimateAlnset, tempName,
                            '--components']
                orblMain()
                ae(streamCatcher.buffer('err'), '')  # Nothing written to stderr
                outputLines = streamCatcher.buffer('out').split('\n')
                for outputLine, expectedLine in zip_longest(outputLines,
                                                            ExpectedOutput.split('\n')) :
                    compare_output_lines(outputLine, expectedLine)
                os.remove(tempName)
                streamCatcher.clear('both')

                ## Test interactive mode
                inputLines = '# First line\nchr1:1\n# comment\n\n# skip this'
                sys.argv = ['orbl.py', PrimateAlnset]
                with StdinSaver(InteractiveStringIO(inputLines)) :
                    orblMain()
                ae(streamCatcher.buffer('err'),
                   'Enter: intervals strand (tab-separated) '
                   '(e.g., chr10:1042727-1042762+chr10:1043301-1043321   -)\n')
                # For interactive intput it should stop at the blank line
                ae(streamCatcher.buffer('out'),
                   '# First line\nchr1:1	InvalidRegion\n# comment\n')
                streamCatcher.clear('both')
                sys.argv = ['orbl.py', '--orblq', PrimateAlnset]
                with StdinSaver(InteractiveStringIO('# Only line\n')) :
                    orblMain()
                ae(streamCatcher.buffer('err'),
                   'Enter: intervals strand biotype+frame_shift (tab-separated) '
                   '(e.g., chr10:1042727-1042762+chr10:1043301-1043321   -   intORF+1)\n')
                ae(streamCatcher.buffer('out'), '# Only line\n')
                streamCatcher.clear('both')

                ## Test calling as a command rather than calling main
                commandName = ORBL_tools.orbl.__file__.replace('.pyc', '.py')
                orblProc = subprocess.Popen([commandName, '-v'], stdout = subprocess.PIPE,
                                            stderr = subprocess.PIPE)
                out, error = orblProc.communicate()
                ae(error.decode('utf-8'), 'Version: %s\n' % currentVersion)
                ae(out.decode('utf-8'), '')

InputLines = ('''# Valid two-interval region; spaces between region/strand/biotype
chr1:3891002-3891057+chr1:3892883-3892889 - intORF+2	c1riboseqorf13
# Invalid region
chr1:3891002-3891055+chr1:3892883	-
# start > end
chr1:3891051-3890992	+
# Missing strand
chr1:3890992-3891051
# Invalid strand
chr1:3890992-3891051	.
# No alignment (no such chromosome)
chr25:3890992-3891051	+
# Blank line (don't stop for non-intractive input)

# Missing biotype
chr1:3890992-3891051	+
# Invalid biotype (space/tab between region/strand/biotype)
chr1:3890992-3891051 +	uORF+1
# Mixed biotype 
chr1:3890992-3891051 +	mixed
# orblv is NA (tab/space between region/strand/biotype)
chr1:3890995-3891051 +	uORF
''')

ExpectedOutput = ('''# Valid two-interval region; spaces between region/strand/biotype
chr1:3891002-3891057+chr1:3892883-3892889	-	intORF+2	c1riboseqorf13	\
0.974690114626	0.902157164869	0.974690114626	1.0	0.974690114626
# Invalid region
chr1:3891002-3891055+chr1:3892883	-	InvalidRegion
# start > end
chr1:3891051-3890992	+	BackwardsInterval
# Missing strand
chr1:3890992-3891051	MissingStrand
# Invalid strand
chr1:3890992-3891051	.	InvalidStrand
# No alignment (no such chromosome)
chr25:3890992-3891051	+	UnableToDownloadAlignment
# Blank line (don't stop for non-intractive input)
	InvalidRegion
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
''')


class TestCalcSeqConservations(unittest.TestCase) :
    """
    Test the checks for various kinds of conservation used in ORBLv.
    """
    def test_calc_seq_conservations(self) :
        self.assertIs(_test_calc_seq_conservations(), None)

if __name__ == '__main__' :
    unittest.main() # Allows running unit test with "python test_TESTNAME.py"