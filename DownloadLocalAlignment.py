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
ThisDir = os.path.dirname(os.path.realpath(__file__)) # Does abspath. Never includes ~.
sys.path.append(os.path.dirname(ThisDir)) # ORBL_tools's parent
from ORBL_tools.NewickTreeUtil import parse_nh_str
from ORBL_tools.CmdLineUtils import assert_num_args

MaxMaxCodons = 40000

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

def download_local_alignment(intervalsStr, strand, alnset, returnFasta = False,
                             removeRefGaps = False) :
    """
    Download a local alignment from CodAlignView and return either a string containing the
        alignment in fasta format (if returnFasta) or a list of pairs (species, bases).
    If removeRefGaps is true, hide inserts and jumps so reference sequence has no gaps.
    """
    url = '%s?i=%s&s=%s&a=%s&m=%d&fo=' % (
        CodAlignViewBaseURL, intervalsStr, strand, alnset, MaxMaxCodons)
    if removeRefGaps :
        url += '&h=on&u=on'
    try :
        fastaStr = get_from_url(url)
    except (HTTPError, URLError) : # pragma: no cover
        raise NotImplementedError('Could not download alignment for %s %s.' %
                                  (intervalsStr, strand))
    # Convert the fasta string to pairs even if returnFasta in order to check validity
    try :
        pairs = parse_fasta_str(fastaStr)
    except NotImplementedError :
        # Probably CodAlignView detected a problem. See if we can extract it.
        # Probably html content. Remove all html tags.
        fileContents = re.sub(r'<[^>]*>', '', fastaStr)
        # Remove "CodAlignView Error"
        fileContents = fileContents.replace('CodAlignView Error\n', '')
        raise NotImplementedError(fileContents)
    return fastaStr if returnFasta else pairs

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

def main() :
    """
    Given an input file containing a list of regions, use CodAlignView to extract the
        local alignment for each region with respect to some Alignment Set (as defined
        in CodAlignView), and save the results as fasta files.
    Input file has one line for each region. That line must contain three tab separated
        fields containing the intervals string, strand, and a unique region name.
        Additional tab separated fields are allowed and will be ignored.
    The format of the intervals string is one or more intervals separated by plus signs:
        chrom:start1-end1+chrom:start2-end2+... satisfying the requirement that:
        start1 <= end1, start2 <= end2, ...  (even if the region is on the minus strand).
        All segments must be on the same chromosome.
    The strand must be + or -.
    Example input line:  chr10:1042727-1042762+chr10:1043301-1043321    +   Name1
    Blank lines and lines that start with a # character will be ignored.
    The output file will be REGION_NAME.fa in the specified output directory.
    """
    assert_num_args(3, 'ALIGNMENT_SET (e.g., hg38_58) IN_FILE OUT_DIR', exact = True)
    alnset, inFileName, outDirName = sys.argv[1:]

    ## Download tree file as a check for valid alnset
    try :
        download_alnset_tree(alnset)
    except NotImplementedError :
        print('ERROR: Invalid Alignment Set "%s".\n' % alnset
              + 'Choose alignment set from list at: %s' %
              'https://data.broadinstitute.org/compbio1/cav.php?Alnsets',
              file = sys.stderr)
        raise SystemExit(1)

    inFileName = os.path.abspath(os.path.expanduser(inFileName))
    outDirName = os.path.abspath(os.path.expanduser(outDirName))
    if not os.path.isdir(outDirName) :
        os.mkdir(outDirName) # Won't handle case when need to add more than one level...
    with open(inFileName, 'rt') as inFile :
        for line in inFile :
            if line.startswith('#') : # Skip comment lines
                continue
            line = line.strip()
            if line == '' :
                continue
            words = line.split('\t')
            if len(words) < 3 :
                print('Invalid input line: %s' % line, file = sys.stderr)
                raise SystemExit(1)
            intervalsStr, strand, regionName = words[:3]
            try :
                fastaStr = download_local_alignment(intervalsStr, strand, alnset,
                                                    returnFasta = True)
            except NotImplementedError :
                print('UnableToDownloadAlignment for line: %s' % line, file = sys.stderr)
                raise SystemExit(1)
            with open(os.path.join(outDirName, regionName + '.fa'), 'wt') as outFile :
                outFile.write(fastaStr)

if __name__ == '__main__' : # pragma: no cover
    main()