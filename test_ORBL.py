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
import sys, os, unittest, tempfile
UnitTestDir = os.path.dirname(os.path.realpath(__file__))
sys.path.append(os.path.dirname(UnitTestDir))# ORBL_tools's parent

from ORBL_tools.orblq import ORBLqCalculator
from ORBL_tools.orbl import main as orblMain
from ORBL_tools.CalcORFrelBLs import _test_calc_seq_conservations
from ORBL_tools.ContextMgrsForTesting import StreamCatcher, SysArgSaver

class TestORBLq(unittest.TestCase) :
    def test_ORBLqCalculator(self) :
        at = self.assertTrue
        ar = self.assertRaises
        aaeq = lambda x, y : at(abs(x - y) < 1e-11)  # "assertAlmostEqual"

        primateCalc = ORBLqCalculator('hg38_120mammals_primate')
        primORBLv = 0.0139718396632
        primORBLq = 0.254531126872
        aaeq(primateCalc('lncRNA-ORF', primORBLv, 92), primORBLq)

        placentalCalc = ORBLqCalculator('hg38_120mammals_placental')
        placORBLv = 0.0642060273392
        placORBLq = 0.657398212512
        matchedORBLvs = []
        aaeq(placentalCalc('uoORF+2', placORBLv, 100, matchedORBLvs), placORBLq)
        at(len(matchedORBLvs) >= 1000)
        at(all(type(ov) == float for ov in matchedORBLvs))
        at(all(0 <= ov <= 1 for ov in matchedORBLvs))

        ar(AssertionError, primateCalc, 'uoORF+3', primORBLv, 92)
        ar(AssertionError, primateCalc, 'mixed',   primORBLv, 92)

class TestORBL(unittest.TestCase) :
    def test_orbl(self) :
        at = self.assertTrue
        ae = self.assertEqual
        aaeq = lambda x, y : at(abs(x - y) < 1e-11)  # "assertAlmostEqual"

        tempHandle, tempName = tempfile.mkstemp(prefix = 'temp', suffix = '.in', dir = None)
        os.close(tempHandle)

        inLine = 'chr1:3891002-3891057+chr1:3,892,883-3,892,889 -\tintORF+2\tc1riboseqorf13'
        outLine1 = inLine.replace(' ', '\t') + '\t0.974690114626'
        outLine2 = outLine1 + '\t0.902157164869'
        outLine3 = outLine2 + '\t0.974690114626\t1.0\t0.974690114626'

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

        with open(tempName, 'wt') as orblInFile :
            print(inLine, file = orblInFile)
        with StreamCatcher('both') as streamCatcher :
            alnset = 'hg38_120mammals_primate'
            with SysArgSaver() :
                sys.argv = ['orbl.py', alnset, tempName]
                orblMain()
                ae(streamCatcher.buffer('err'), '') # Nothing written to stderr
                compare_output_lines(streamCatcher.buffer('out').rstrip('\n'), outLine1)
                streamCatcher.clear('both')

                sys.argv = ['orbl.py', alnset, '--orblq', tempName]
                orblMain() # Warning: removes the --orblq from sys.argv
                ae(streamCatcher.buffer('err'), '') # Nothing written to stderr
                compare_output_lines(streamCatcher.buffer('out').rstrip('\n'), outLine2)
                streamCatcher.clear('both')

                sys.argv = ['orbl.py', alnset, '--orblq', tempName, '--components']
                orblMain()
                ae(streamCatcher.buffer('err'), '') # Nothing written to stderr
                compare_output_lines(streamCatcher.buffer('out').rstrip('\n'), outLine3)
                streamCatcher.clear('both')

        os.remove(tempName)

class TestCalcSeqConservations(unittest.TestCase) :
    """
    Test the checks for various kinds of conservation used in ORBLv.
    """
    def test_calc_seq_conservations(self) :
        self.assertIs(_test_calc_seq_conservations(), None)

if __name__ == '__main__' :
    unittest.main() # Allows running unit test with "python test_TESTNAME.py"