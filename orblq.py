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
Tool for computing ORBLq, a measure of evolutionary constraint on the "ORFness" of an ORF.
"""
from __future__ import division, print_function
import os, pickle

BiotypesWithFS = ['uORF', 'uoORF+1', 'uoORF+2', 'intORF+1', 'intORF+2', 
                  'doORF+1', 'doORF+2', 'dORF', 'lncRNA-ORF']

UntransORBLvDictsDir = os.path.join(os.path.dirname(os.path.realpath(__file__)),
                                    'UntransORBLvDicts')

def get_untrans_ORBLv_dict_name(alnset) :
    return os.path.join(UntransORBLvDictsDir, 'UntransORBLvDict.%s.cp' % alnset)

def get_supported_alnsets() :
    untransORBLvDictFileNames = os.listdir(UntransORBLvDictsDir)
    supportedAlnsets = [fname[len('UntransORBLvDict.') : -3]
                        for fname in untransORBLvDictFileNames
                        if fname.startswith('UntransORBLvDict.') and
                        fname.endswith('.cp')]
    return supportedAlnsets

class ORBLqCalculator(object) :
    """
    Class for computing ORBLq, a measure of evolutionary constraint on the "ORFness"
        of an ORF. It is the quantile of the ORBLv (ORF Relative Branch Length
        conservation score) among 1000 or more ORFs of the same biotype (and frame shift
        relative to the main frame for overlapping biotypes) and similar length that are
        not believed to be translated.
    Requires previous creation of a dictionary mapping each biotypeWithFS to an
        OrderedDict mapping ORFlength to a list of ORBLv values of untranslated ORFs of
        that biotypeWithFS and length. The OrderedDict is sorted by increasing ORFlength.
    ORBLv and ORBLq are defined relative to a particular "Alignment Set" as defined by
        CodAlignView, https://data.broadinstitute.org/compbio1/cav.php. The predefined
        dictionary mapping biotypeWithFS to the OrderedDict of ORBLv values by length
        must be stored in a pickle file as specified by get_untrans_ORBLv_dict_name.
    """
    def __init__(self, alnset) :
        with open(get_untrans_ORBLv_dict_name(alnset), 'rb') as infile :
            self.untransORBLvDict = pickle.load(infile)
            # biotypeWithFS : OrderedDict((orfLen, [orblv, orblv, ...]), ...)
    def __call__(self, biotypeWithFS, ORBLv, numCodonsWithStop, matchedORBLvs = None) :
        """
        Return ORBLq of an ORF with specified biotypeWithFS, ORBLv, and numCodonsWithStop.
        If matchedORBLvs is not None, it should be an empty list, which will be used to
            return the ORBLv's of the matched transcripts.
        """
        assert biotypeWithFS in BiotypesWithFS
        if matchedORBLvs is not None :
            assert len(matchedORBLvs) == 0

        untransORBLvListsByLen = self.untransORBLvDict[biotypeWithFS]
            
        # Sort ORBLv lists by absolute value of difference from numCodonsWithStop.
        sortedItems = sorted(untransORBLvListsByLen.items(), # item is orfLen, ORBLvList
                             key = lambda item : (abs(item[0] - numCodonsWithStop),
                                                      item[0] - numCodonsWithStop))
        
        # Make a list of lists of ORBLv values, where each of the sublists consists of
        #   ORBLv values from an orf length with the same absolute difference from
        #   numCodonsWithStop.
        sortedORBLvlists = []
        ind = 0
        while ind < len(sortedItems) - 1 :
            orfLen     = sortedItems[ind    ][0]
            nextOrfLen = sortedItems[ind + 1][0]
            if abs(orfLen - numCodonsWithStop) == abs(nextOrfLen - numCodonsWithStop) :
                sortedORBLvlists.append(sortedItems[ind][1] + sortedItems[ind + 1][1])
                ind += 2
            else :
                sortedORBLvlists.append(sortedItems[ind][1])
                ind += 1
        if ind == len(sortedItems) - 1 :
            sortedORBLvlists.append(sortedItems[ind][1])

        # Compare ORBLv to at least 1000 ORBLv's from untranslated ORFs of length as close
        # as possible to numCodonsWithStop. If taking _any_ ORFs whose length is a given
        # absolute distance from numCodonsWithStop, take _all_ such ORFs.
        minNumToCheck = 1000
        numChecked = numORBLvNotLess = 0
        for sameDistORBLvs in sortedORBLvlists :
            if numChecked >= minNumToCheck :
                break
            # Subtract 1e-11 to avoid instability due to floating point errors in ORBLv
            # Note: smallest difference between distinct values in untransORBLv's for
            #      any biotype for hg38_120mammals_placental is about 1.11e-9, so 1e-11
            #      should be small enough to avoid affecting any real differences.
            numORBLvNotLess += sum(orblv >= ORBLv - 1e-11 for orblv in sameDistORBLvs)
            numChecked += len(sameDistORBLvs)
            if matchedORBLvs is not None :
                matchedORBLvs.extend(sameDistORBLvs)
        assert numChecked >= minNumToCheck, numChecked
        # Add a pseudocount of 1 to avoid ORBLq == 1
        pseudocount = 1
        if numORBLvNotLess < numChecked : # Check this to avoid ORBLq < 0
            numORBLvNotLess += pseudocount
        return 1 - numORBLvNotLess / numChecked
