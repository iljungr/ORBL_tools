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
if sys.version_info[0] < 3 : # pragma: no cover
    from urllib2 import urlopen, HTTPError, URLError
else : # pragma: no cover
    from urllib.request import urlopen
    from urllib.error import HTTPError, URLError
from .NewickTreeUtil import parse_nh_str

CodAlignViewBaseURL = 'https://data.broadinstitute.org/compbio1/cav.php'
TreeBaseURL         = 'https://data.broadinstitute.org/compbio1/CodAlignViewFiles/TreeNHs/'

def get_from_url(url) :
    """Visit the URL and return the resulting string, ensuring socket is closed."""
    response = urlopen(url)
    try:
        result = response.read().decode('utf-8') # Read returns bytes. Convert to str.
    finally:
        response.close()  # ensure socket is closed to avoid ResourceWarning
    return result

def download_local_alignment(intervalsStr, strand, alnset, removeRefGaps = False,
                             maxCodons = None) :
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
    try :
        fastaStr = get_from_url(url)
    except (HTTPError, URLError) :
        raise NotImplementedError('Could not download alignment for %s %s.' %
                                  (intervalsStr, strand))
    try :
        aliSeg = parse_fasta_str(fastaStr)
    except NotImplementedError :
        # Probably CodAlignView detected a problem. See if we can extract it.
        # Probably html content. Remove all html tags.
        fileContents = re.sub(r'<[^>]*>', '', fastaStr)
        # Remove "CodAlignView Error"
        fileContents = fileContents.replace('CodAlignView Error\n', '')
        raise NotImplementedError(fileContents)
    return aliSeg

def download_alnset_tree(alnset) :
    """
    Download and parse tree file representing the phylogenetic tree of species in alnset.
    Return a NewickTreeUtil.NewickNode containing the root of the tree.
    """
    url = TreeBaseURL + alnset + '.nh'
    try :
        treeString = get_from_url(url).rstrip()
    except (HTTPError, URLError) :
        raise NotImplementedError('Could not download tree file for %s.' % alnset)
    tree = parse_nh_str(treeString)
    return tree

def parse_fasta_str(fastaStr) :
    """
    Read string containing the contents of a file in fasta format, and return a list of
        pairs: (sequence name, sequence)
    """
    pairs = []
    seqName = None
    seq = ''
    for line in fastaStr.split('\n') :
        line = line.rstrip()
        if len(line) == 0 :
            continue
        if line[0] == '>' :
            if seqName is not None :
                pairs.append((seqName, seq))
            seqName = line[1:]
            seq = ''
        elif seqName is None :
            raise NotImplementedError('Not valid fasta format.')
        else :
            seq += line
    if seqName is not None :
        pairs.append((seqName, seq))
    return pairs
