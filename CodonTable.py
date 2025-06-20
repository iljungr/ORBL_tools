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

STOP = '*'

def translate(baseString, allowInvalid = False, sec = False, mitoGroup = None) :
    """
    Translate multiple-of-3-length string of bases. Stops translate to STOP (*).
    Ignore white space in the input string.
    Bases may be any of A, C, G, T, U, a, c, g, t, u.
    If allowInvalid then also allow other characters in bases and
      translate any codon including another character to '#'.
    If sec, translate TGA as U rather than *.
    Allowed mitoGroups: 'vertebrate', 'yeast', 'protozoan', 'invertebrate', 'echinoderm',
        or None. If not None, translate according to appropriate mitochondrial code.
    """
    baseString = ''.join(baseString.upper().replace('U', 'T').split())
    aas = ''
    mitoTable = get_mito_table(mitoGroup)
    for i in range(0, len(baseString) - len(baseString) % 3, 3) :
        codon = baseString[i : i + 3]
        if sec and codon == 'TGA' :
            aa = 'U'
        elif codon in mitoTable :
            aa = mitoTable[codon]
        else :
            aa = _CodonTable.get(codon, '#')
        aas += aa
        assert allowInvalid or aa != '#', codon
    return aas

def transltOld(baseString) :
    return translate(baseString).replace(STOP, 'Stop')
    
def is_stop(baseString, sec = False, mitoGroup = None):
    return (len(baseString) == 3 and
            translate(baseString, allowInvalid = True,
                      sec = sec, mitoGroup = mitoGroup) == '*')

def is_ORF(bases, sec = False, mitoGroup = None) :
    """
    Return True if the bases form an ORF: multiple of 3, starts with ATG, ends with stop,
        no other in-frame stop codons, no unknown bases including '.' except ignore -'s.
    """
    if '.' in bases :
        return False # No way to tell if in-frame
    bases = bases.replace('-', '')
    if len(bases) % 3 :
        return False # Stop is not in frame
    aas = translate(bases, allowInvalid = True, sec = sec, mitoGroup = mitoGroup)
    if '#' in aas :
        return False # Some bases are unknown; could be stop
    if aas[0] != 'M' :
        return False # Start is not start codon
    if aas[-1] != '*' :
        return False # End is not stop codon
    if '*' in aas[:-1] :
        return False # Early stop codon
    return True

def get_AAs() :
    ## Return the set of amino acid 1-letter codes
    return set([val for val in _CodonTable.values() if val != STOP])
    
def get_aa_info(nameOrAbbreviation = None) :
    """
    Input: amino acid 1 letter code, 3 letter code, name, codon, any subset of name, or
        None.
    Tries to find a (case-insensitive) exact match to amino acid 1 letter code, 3 letter
        code, name, or codon. This match will be unique if it exists, because no amino
        acid 3 letter codes are also codons. If there is no exact match, then looks for
        subset of amino acid name.
    Output: [Amino Acid name, 3 letter code, 1 letter code, [all codons]]
             or, if input is None, list for all (non-degenerate) amino acids.
    Example: get_aa_info('k') -> ['Lysine', 'Lys', 'K', ['AAA', 'AAG']]
    """
    abbreviations = AminoAcidAbbreviations + DegenerateAminoAcidAbbreviations
    if nameOrAbbreviation is None :
        return [get_aa_info(a) for a in sorted(get_AAs())]
    if nameOrAbbreviation.upper().replace('U', 'T') in _CodonTable :
        searchString = translate(nameOrAbbreviation)
    else :
        searchString = nameOrAbbreviation[0].upper() + \
                   nameOrAbbreviation[1 :].lower()
    for triple in abbreviations :
        if searchString in triple :
            return list(triple) + [sorted([codon 
                                           for codon in _CodonTable
                                           if _CodonTable[codon] == triple[2]])]
    for triple in abbreviations :
        if searchString.lower() in triple[0].lower() :
            return list(triple) + [sorted([codon 
                                           for codon in _CodonTable
                                           if _CodonTable[codon] == triple[2]])]
    return ''
aa_info = get_aa_info

# Standard Codon Table:
_CodonTable = {
    'TTT': 'F', 
    'TTC': 'F', 
    'TTA': 'L', 
    'TTG': 'L', 
    'CTT': 'L', 
    'CTC': 'L', 
    'CTA': 'L', 
    'CTG': 'L',
    'ATT': 'I', 
    'ATC': 'I', 
    'ATA': 'I', 
    'ATG': 'M', 
    'GTT': 'V', 
    'GTC': 'V', 
    'GTA': 'V', 
    'GTG': 'V', 
    'TCT': 'S', 
    'TCC': 'S', 
    'TCA': 'S', 
    'TCG': 'S', 
    'CCT': 'P', 
    'CCC': 'P', 
    'CCA': 'P', 
    'CCG': 'P', 
    'ACT': 'T', 
    'ACC': 'T', 
    'ACA': 'T', 
    'ACG': 'T', 
    'GCT': 'A', 
    'GCC': 'A', 
    'GCA': 'A', 
    'GCG': 'A', 
    'TAT': 'Y', 
    'TAC': 'Y', 
    'TAA': STOP, # Ochre
    'TAG': STOP, # Amber
    'CAT': 'H', 
    'CAC': 'H', 
    'CAA': 'Q', 
    'CAG': 'Q', 
    'AAT': 'N', 
    'AAC': 'N', 
    'AAA': 'K', 
    'AAG': 'K', 
    'GAT': 'D', 
    'GAC': 'D', 
    'GAA': 'E', 
    'GAG': 'E', 
    'TGT': 'C', 
    'TGC': 'C', 
    'TGA': STOP, # Opal (will translate to 'U' if sec)
    'TGG': 'W', 
    'CGT': 'R', 
    'CGC': 'R', 
    'CGA': 'R', 
    'CGG': 'R', 
    'AGT': 'S', 
    'AGC': 'S', 
    'AGA': 'R', 
    'AGG': 'R', 
    'GGT': 'G', 
    'GGC': 'G', 
    'GGA': 'G', 
    'GGG': 'G' 
    }

# Mitochondrial codon tables from https://www.ncbi.nlm.nih.gov/Taxonomy/Utils/wprintgc.cgi
def get_mito_table(group) :
    MitoTables = {
        'vertebrate' : _CodonTableVertebrateMito,
        'yeast' : _CodonTableYeastMito,
        'protozoan' : _CodonTableProtozoanMito,
        'invertebrate' : _CodonTableInvertebrateMito,
        'echinoderm' : _CodonTableEchinodermMito,
        None : {}
        }
    assert group in MitoTables, 'Name of mito table must be in %s' % list(MitoTables)
    return MitoTables[group]

_CodonTableVertebrateMito = {
    'AGA' : STOP,
    'AGG' : STOP,
       # Actually, AGA and AGG don't stop; instead they cause a -1 frame shift and since
       # previous base is T in the genes that have them (1 each in human), that results
       # in a stop.
    'ATA' : 'M',
    'TGA' : 'W',
    }
_CodonTableYeastMito = {
    'ATA' : 'M',
    'CTA' : 'T',
    'CTC' : 'T',
    'CTG' : 'T',
    'CTT' : 'T',
    'TGA' : 'W',
    }
_CodonTableProtozoanMito = { # Mold, Protozoan, and Coelenterate
    'TGA' : 'W',
    }
_CodonTableInvertebrateMito = { # Includes Drosophila, Apis, and Caenorhabditis, at least
    'AGA' : 'S',
    'AGG' : 'S', # In some arthropods it is 'K'
    'ATA' : 'M',
    'TGA' : 'W',
    }
_CodonTableEchinodermMito = { # Echinoderm and Flatworm
    'AAA' : 'N',
    'AGA' : 'S',
    'AGG' : 'S',
    'TGA' : 'W',
    }
# Notes:
# - There are other codon tables, including one for Candida albicans nuclear genes.
# - Bacteria, archaea, and plant plastids/chloroplasts use the standard genetic code,
# except they sometimes have more non-cognate start codons (which I guess is due to
# incomplete annotation of start codons in the standard code)
# https://www.ncbi.nlm.nih.gov/Taxonomy/Utils/wprintgc.cgi#SG11

AminoAcidAbbreviations = [
    ('Alanine', 'Ala', 'A'),
    ('Arginine', 'Arg', 'R'),
    ('Asparagine', 'Asn', 'N'),
    ('Aspartic acid', 'Asp', 'D'),
    ('Cysteine', 'Cys', 'C'),
    ('Glutamic acid', 'Glu', 'E'),
    ('Glutamine', 'Gln', 'Q'),
    ('Glycine', 'Gly', 'G'),
    ('Histidine', 'His', 'H'),
    ('Isoleucine', 'Ile', 'I'),
    ('Leucine', 'Leu', 'L'),
    ('Lysine', 'Lys', 'K'),
    ('Methionine', 'Met', 'M'),
    ('Phenylalanine', 'Phe', 'F'),
    ('Proline', 'Pro', 'P'),
    ('Pyrrolysine', 'Pyl', 'O'),
    ('Selenocysteine', 'Sec', 'U'),
    ('Serine', 'Ser', 'S'),
    ('Threonine', 'Thr', 'T'),
    ('Tryptophan', 'Trp', 'W'),
    ('Tyrosine', 'Tyr', 'Y'),
    ('Valine', 'Val', 'V'),
    ]
DegenerateAminoAcidAbbreviations = [ # from https://www.genome.jp/kegg/catalog/codes1.html
    ('Asn or Asp', 'Asx', 'B'),
    ('Gln or Glu', 'Glx', 'Z'),
    ('Unknown', 'Unk', 'X'),
]


    
