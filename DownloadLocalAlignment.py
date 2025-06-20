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
Download a local alignment or species tree from CodAlignView.
"""
from __future__ import division, print_function

import sys, os, re
from .NewickTreeUtil import parse_nh_str
if sys.version_info[0] < 3 : # pragma: no cover
    from urllib import urlretrieve
else : # pragma: no cover
    from urllib.request import urlretrieve

CodAlignViewBaseURL = 'https://data.broadinstitute.org/compbio1/cav.php'
TreeBaseURL         = 'https://data.broadinstitute.org/compbio1/CodAlignViewFiles/TreeNHs/'

def download_local_alignment(intervalsStr, strand, alnset, removeRefGaps = False,
                             maxCodons = None, debug = False) :
    """
    Download a local alignment from CodAlignView and return an aliSeg (a list of pairs
        (species, bases).
    If removeRefGaps is true, hide inserts and jumps so reference sequence has no gaps.
    If maxCodons is not None, allow regions up to maxCodons codons, otherwise region will
        be limited to CodAlignView's default (currently 1000 codons).
    """
    url = '%s?i=%s&s=%s&a=%s&fo=' % (CodAlignViewBaseURL, intervalsStr, strand, alnset)
    if removeRefGaps :
        url += '&h=on&u=on'
    if maxCodons is not None :
        url += '&m=%d' % maxCodons
    tempFileName = urlretrieve(url)[0]

    try :
        aliSeg = read_fasta(tempFileName)
    except NotImplementedError :
        if debug :
            print('Error in %s' % tempFileName, file = sys.stderr)

        # Probably CodAlignView detected a problem. See if we can extract it.
        fileContents = open(tempFileName).read()
        # Probably html content. Remove all html tags.
        fileContents = re.sub(r'<[^>]*>', '', fileContents)
        # Remove "CodAlignView Error"
        fileContents = fileContents.replace('CodAlignView Error\n', '')

        raise NotImplementedError(fileContents)
    os.remove(tempFileName)
    return aliSeg

def download_alnset_tree(alnset) :
    """
    Download and parse tree file representing the phylogenetic tree of species in alnset.
    Return a NewickTreeUtil.NewickNode containing the root of the tree.
    """
    url = TreeBaseURL + alnset + '.nh'
    tempFileName = urlretrieve(url)[0]

    with open(tempFileName, 'rt') as tempFile :
        treeString = tempFile.read().rstrip()
        if 'not found on this server' in treeString :
            raise NotImplementedError('No tree file for %s.' % alnset)

    tree = parse_nh_str(treeString)
    os.remove(tempFileName)
    return tree

def read_fasta(faFileName) :
    """Read a fasta file returning list of pairs: (sequence name, sequence)"""
    pairs = []
    with open(faFileName) as faFile :
        seqName = None
        seq = ''
        for line in faFile :
            line = line.rstrip()
            if len(line) == 0 :
                continue
            if line[0] == '>' :
                if seqName is not None :
                    pairs.append((seqName, seq))
                seqName = line[1:]
                seq = ''
            elif seqName is None :
                raise NotImplementedError('Not a valid fasta file: %s' % faFileName)
            else :
                seq += line
        if seqName is not None :
            pairs.append((seqName, seq))
    return pairs
