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
from __future__ import division, print_function
import sys, os
from .CodonTable import translate

def calc_ORF_relBLs(aliSeg, tree, secondORF = False, nTermExt = False, fullBL = None,
                    mitoGroup = None, returnSpecLists = False) :
    """
    Determine the relative branch length of the species present in the local alignment
        of an ORF (or secondORF, or N-terminal extension) that satisfy various conditions.
    Return a tuple whose members depend on secondORF and nTermExt as specified below.
    aliSeg: An alignment segment, which is a list of two element lists representing a
        multispecies genomic alignment. Each two element list consists of the species name
        and a string of bases (the same number for each species). In addition to
        the four DNA nucleotides (in upper case), the base string can include '-'
        representing a gap relative to the other species, '.' meaning that the there is no
        alignment there (i.e., the homologous sequence is unknown or non-existent), or '|'
        meaning there is a "jump", i.e., the coordinates on each side are not consecutive.
        The first species is the reference species. The region representing the ORF should
        include the start and stop codons (or both stop codons for second ORF, or both
        start codons for nTermExt).
    tree: a NewickTreeUtil.NewickNode representing the root of a phylogenetic tree that
        includes the species in the alignment. Returned values will be relative to the
        full branch length of the tree, so typically it should be the full tree of the
        whole genome alignment, not just the species present in aliSeg.
        (Call tree.set_descendants() once before passing in tree multiple times to improve
        speed of repeated calculations of subtree branch length.)
    secondORF: Region is a second ORF including stop codons at both ends.
    nTermExt: Region is an N-terminal extension including both first and second start.
        Not allowed if secondORF.
    fullBL: Branch length of the full tree. If None, will be determined from tree, but may
        be passed in to improve performance of repeated calls.
    mitoGroup: Codon table group for calls to translate if mitochondrial else None. We are
        assuming the codon table is the same for all rows in the alignment, which it won't
        be if a mitochondrial gene migrates to the nucleus.
    In default case (i.e., not secondORF or nTermExt), return a tuple with the following:
        ATGRelBL   # RelBL of ATG in 1st codon; NA if ref species is not ATG
        StopRelBL  # RelBL of stop in last codon; NA if ref species is not stop;
        FrameRelBL # RelBL of length being a multiple of 3 and no premature in-frame stops
        ORFRelBL   # RelBL of species satisfying all three (NA if ref species doesn't)
        NearCognateRelBL # RelBL of ATG or near cognate in 1st codon; NA if ref spec isn't
        NearCognateORFRelBL # RelBL of species satisfying all 3; NA if ref species doesn't
        Note: by "ATG" we mean any codon that translates to M, taking mitoGroup into
              account, so for vertebrate mitochondria it includes ATA.
    If secondORF, return a tuple with the following:
        Stop1RelBL # RelBL of stop in 1st codon of region; NA if ref species is not stop
        SameStop1RelBL # RelBL of same 1st stop as reference; NA if ref not stop
        Stop2RelBL # RelBL of stop in codon at end; NA if ref species is not stop
        FrameRelBL # RelBL of species with length multiple of 3 and no in-frame stops
        ORFRelBL # RelBL of species with both stops and frame (NA if ref species doesn't)
    If nTermExt, return a tuple with the following:
        ATGRelBL         # RelBL of ATG in 1st codon; NA if ref species is not ATG
        SameStartRelBL   # RelBL of same (1st) start codon as ref species
        NearCognateRelBL # RelBL of ATG or near cognate in 1st codon; NA if ref spec isn't
        FrameRelBL      # RelBL of species with length multiple of 3 and no in-frame stops
        ATGExtRelBL         # RelBL of FrameRelBL and ATGRelBL, or NA
        SameStartExtRelBL   # RelBL of FrameRelBL and SameStartRelBL, or NA
        NearCognateExtRelBL # RelBL of FrameRelBL and NearCognateRelBL, or NA
    returnSpecLists: if true, append to the tuple comma-separated lists of the conserved
        species for each of the relative branch lengths being returned (for debugging).
    Notes:
    - If a species occurs more than once in the alignment (one-to-many alignment), each
        row of the alignment is processed separately for each RelBL calculation and
        a species is considered to satisfy the condition if any of the rows
        for that species satisfies the condition. That means it is possible that, for
        example, one row could have an aligned ATG, a different row an aligned stop,
        and a third row a conserved frame, so the species would be counted in
        ATGRelBL, StopRelBL, and FrameRelBL, but not ORFRelBL because no individual row
        satisfied all three conditions.
    - See comment in _calc_seq_conservations header about where frame calculation starts
        and ends, particularly in the case of first and last codons with more than three
        nucleotides in the aligned region; also for handling of ATG and stop for
        non-standard genetic codes.
    """
    if fullBL is None :
        fullBL = tree.subtree_branch_length()

    conservedSpecies = {consType : [] for consType in ConsTypes}

    # Calculate the 0-based column index of the last base of the first codon and the
    # first base of the last codon
    refBases = aliSeg[0][1]
    firstCodonLastInd = _index_of_kth_non_gap(refBases, 3, fromEnd = False)
    lastCodonFirstInd = _index_of_kth_non_gap(refBases, 3, fromEnd = True)

    refSpec = aliSeg[0][0]
    refStartCodon = refBases[: firstCodonLastInd + 1].replace('-', '')
    assert len(refStartCodon) == 3
    for aliRow in aliSeg :
        bases = aliRow[1]
        seqCons = _calc_seq_conservations(bases, firstCodonLastInd, lastCodonFirstInd,
                                          refStartCodon, secondORF, nTermExt, mitoGroup)
        for consType in ConsTypes :
            if consType in seqCons :
                conservedSpecies[consType].append(aliRow[0])

    # Note that subtree_branch_length converts the species list to a set, so if a species
    # is included more than once (in a 1-to-many alignment), it will be counted only once.
    def _calc_rel_bl(specList) :
        return ('NA' if refSpec not in specList else
                tree.subtree_branch_length(specList) / fullBL)
    relBLs = {consType : _calc_rel_bl(conservedSpecies[consType])
              for consType in ConsTypes}

    returnTypes = get_return_types(secondORF, nTermExt)
    result = tuple(relBLs[consType] for consType in returnTypes)
    if returnSpecLists :
        result += tuple(','.join(conservedSpecies[consType]) for consType in returnTypes)

    return result

def get_return_types(secondORF, nTermExt) :
    """
    Return the list of kinds of relative branch length calc_ORF_relBLs will return.
    """
    if secondORF :
        # Note that we are not returning SameStop1SecondORF, but it can be calculated
        #   as SameStop1 and SecondORF.
        return ['Stop1', 'SameStop1', 'StopAtEnd', 'Stop1Frame', 'SecondORF']
    elif nTermExt :
        # Note that there is loss of information here, because ATGFrame, NearCognateFrame,
        # and SameStartFrame could be different, but we are only returning
        # NearCognateFrame. The case that they differ should be sufficiently rare that
        # it is not worth changing this, for now.
        return ['ATG', 'SameStart', 'NearCognate', 'NearCognateFrame',
                'ATGExt', 'SameStartExt', 'NearCognateExt']
    else :
        # Same note as above regarding loss of information.
        return ['ATG', 'StopAtEnd', 'NearCognateFrame', 'ORF',
                'NearCognate', 'NearCognateORF']

def _index_of_kth_non_gap(bases, k, fromEnd = False) :
    """Calculate the 0-based index of the kth non-gap base from the beginning or end."""
    assert k > 0
    numNonGap = 0
    numBases = len(bases)
    columnInd = 0 # To suppress dumb warning
    for columnInd in (range(numBases - 1, -1, -1) if fromEnd else range(numBases)) :
        if bases[columnInd] in '.-|' :
            continue
        numNonGap += 1
        if numNonGap >= k :
            break
    assert numNonGap == k, 'Fewer than %s non-gap bases' % k
    return columnInd

ConsTypes = [ # Types of conservation related to ORF conservation
    # Each makes sense only for certain values of secondORF and nTermExt
    'ATG',         # Is there an ATG in sequence aligned to ref start codon
    'NearCognate', # Is there a near cognate in bases aligned to ref start codon
    'SameStart',   # Is there same codon aligned to ref start, needn't be near cognate
    'Stop1',       # Is there a stop in bases aligned to ref 1st codon
    'SameStop1',   # Is there same stop in bases aligned to ref 1st codon; must be a stop
    'StopAtEnd',   # Is there a stop in bases aligned to ref final codon

    # The following report if base string from position in the region aligned to the 1st
    #     codon to position in region aligned to last codon is multiple of 3, without
    #     .'s, |'s, or in-frame stop codons. Base string starts at specified codon if
    #     present, otherwise start of first codon. Base string ends at stop codon if
    #     present and not nTermExt, otherwise end of last codon.
    'ATGFrame',         # From ATG
    'NearCognateFrame', # From near cognate codon
    'SameStartFrame',   # From codon matching reference codon
    'Stop1Frame',       # From stop codon
    # Note: there's no need for SameStop1Frame because we require region before stop1 to
    #     be a multiple of 3, so all stops are in the same frame

    'ORF',             # ATGFrame         and ATG         and StopAtEnd
    'NearCognateORF',  # NearCognateFrame and NearCognate and StopAtEnd
    'SameStartORF',    # SameStartFrame   and SameStart   and StopAtEnd

    'SecondORF',          # Stop1     and Stop1Frame and StopAtEnd
    'SameStop1SecondORF', # SameStop1 and Stop1Frame and StopAtEnd

    # The following report if valid N-terminal extension; they do NOT require ATG in last
    #    codon.
    'ATGExt',          # ATGFrame         and ATG
    'SameStartExt',    # SameStartFrame   and SameStart
    'NearCognateExt',  # NearCognateFrame and NearCognate

# See comment in _calc_seq_conservations header about where frame calculation starts and
#     ends, and meaning of ATG and stop for non-standard codon tables.
]

OrderConsTypes = lambda result : sorted(result, key = lambda s : ConsTypes.index(s))

ConsTypesForFirstORF  = ['ATG', 'NearCognate', 'SameStart',
                         'StopAtEnd',
                         'ATGFrame', 'NearCognateFrame', 'SameStartFrame',
                         'ORF', 'NearCognateORF', 'SameStartORF']
ConsTypesForSecondORF = ['Stop1', 'SameStop1',
                         'StopAtEnd',
                         'Stop1Frame',
                         'SecondORF', 'SameStop1SecondORF']
ConsTypesForNTermExt  = ['ATG', 'NearCognate', 'SameStart',
                         'ATGFrame', 'NearCognateFrame', 'SameStartFrame',
                         'ATGExt', 'SameStartExt', 'NearCognateExt']
assert all(consType in ConsTypes for consType in (ConsTypesForFirstORF +
                                                  ConsTypesForSecondORF +
                                                  ConsTypesForNTermExt))

def is_near_cognate(codon) :
    return (len(codon) == 3 and
            sum(b1 != b2 for b1, b2 in zip(codon, 'ATG')) <= 1 and
            codon[1] not in '.|' and codon[2] not in '.|')
    # NTG is always near cognate, so it's OK if codon[0] is '.' or '|'
    # In principle, the same is true for ATN, but it is easier to check for frame
    #     conservation if we don't allow AT. and AT|
    # Note that is_near_cognate includes all "ATG" codons, i.e., ones that translate to
    #     M, even with alternate codon tables, because (at least among the currently
    #     handled ones), the only other codon that translates to M is ATA, which is
    #     near cognate.

def _calc_seq_conservations(bases, firstCodonLastInd, lastCodonFirstInd,
                            refStartCodon = 'ATG', secondORF = False, nTermExt = False,
                            mitoGroup = None) :
    """
    Determine the conservation of bases with regard to each of the conservation types
        in ConsTypes.
    Input:
        bases: string of bases (or ., -, or |), aligned to reference ORF (including start
            and stop codon, if present), or to second ORF (including stops at both ends),
            or to N-terminal extension (including both ATGs).
        numBasesInFirstCodon: the number of bases aligned to the reference first codon
        numBasesInLastCodon:  the number of bases aligned to the reference final codon
        refStartCodon: first three bases of reference sequence (no gaps)
        secondORF: ORF conservation requires conservation of Stop1 rather than ATG
        nTermExt: consider frame conservation up to the end of bases, rather than to
            a stop codon aligned to the last codon
        mitoGroup:  as defined in translate(). "ATG" and "stop" mean whatever translates
            to M and * with this mitoGroup.
    Return a set with the ConsTypes that are conserved. Only return the ones that make
        sense based on whether secondORF or nTermExt are True.

    If more than three nucleotides are aligned to an ATG, Stop, or NearCognate codon
        in the reference species, consider it to have an aligned ATG etc if the
        aligned string includes an ATG etc. Consider the frame to be conserved if the
        length from this ATG (or any of them if more than one) to this stop (or any if
        more than one) is a multiple of three (and no intervening .'s, |'s, or
        in-frame stops). If there is no aligned ATG (or stop), consider frame from the
        first base aligned to the reference ATG (or to the last base aligned to the
        reference stop codon). (This means that frame conservation might differ for ATG,
        near cognate, and same start.) For stop1, only consider stops in the region
        aligned to the reference stop1 for which the earlier portion is a multiple of 3
        with no .'s, |'s, or in-frame stops. In the case of nTermExt, since there is no
        requirement for the final codon to be ATG, for frame conservation go all the way
        to the end of the final codon.
    Note that by ATG we mean any codon that translates to M, which for vertebrate
        mitochondria includes ATA, and stop means any codon that translates to *.
    """
    assert not (secondORF and nTermExt)
    firstORF = not secondORF and not nTermExt

    def local_translate(baseStr) :
        return translate(baseStr, allowInvalid = True, mitoGroup = mitoGroup)

    def conserved_frame(bs) :
        """ Is bs a multiple of 3 and doesn't contain ., |, or in-frame stops? """
        assert '-' not in bs
        return (len(bs) % 3 == 0 and '.' not in bs and '|' not in bs and
                '*' not in local_translate(bs))

    # Remove dashes, and adjust firstCodonLastInd and lastCodonFirstInd
    firstCodonLastInd = sum(b != '-' for b in bases[: firstCodonLastInd + 1]) - 1
    lastCodonFirstInd = sum(b != '-' for b in bases[: lastCodonFirstInd])
    bases = bases.replace('-', '')
    refStartCodon = refStartCodon.replace('-', '')

    startBases = bases[: firstCodonLastInd + 1]
    endBases   = bases[lastCodonFirstInd :]

    atgOffets = [ii for ii in range(len(startBases) - 2)
                 if 'M' == local_translate(startBases[ii : ii + 3])]
    hasATG = len(atgOffets) > 0

    sameStartOffsets = [ii for ii in range(len(startBases) - 2)
                        if startBases[ii : ii + 3] == refStartCodon]
    hasSameStart = len(sameStartOffsets) > 0

    nearCognateOffsets = [ii for ii in range(len(startBases) - 2)
                          if is_near_cognate(startBases[ii : ii + 3])]
    hasNearCognate = len(nearCognateOffsets) > 0

    stop1offsets = [ii for ii in range(0, len(startBases) - 2, 3)
                    if '*' == local_translate(startBases[ii : ii + 3]) and
                    conserved_frame(startBases[: ii])]
                    # No need to check for .|* in [ii+3:]; will be in check for frame
    hasStop1 = len(stop1offsets) > 0

    sameStop1offsets = [ii for ii in stop1offsets if ii in sameStartOffsets]
    hasSameStop1 = len(sameStop1offsets) > 0

    stopAtEndOffsets = [ii for ii in range(len(endBases))
                        if '*' == local_translate(endBases[ii : ii + 3])]
    hasStopAtEnd = len(stopAtEndOffsets) > 0

    """
    Frame conservation depends on where the ORF/2nd-ORF/extension starts and ends, and
        if the region aligned to the first codon has more than 3 bases that can have 
        different answers depending on whether we are looking for an ATG, same start, or
        any near cognate-initiated ORF. We therefore need several consFrames. Then
        the various ORF conservations will have to check the appropriate consFrame.
    """
    firstPosesAfterOrfs = ([len(bases)] if nTermExt or not hasStopAtEnd else
                            # All positions just before final stop codon
                            [lastCodonFirstInd + ii for ii in stopAtEndOffsets])
    def any_conserved_frame(startOffsets) :
        frameStarts = [offset + 3 for offset in startOffsets] if startOffsets else [0]
        return any(conserved_frame(bases[start : end])
                   for start in frameStarts
                   for end in firstPosesAfterOrfs)
    consFrameATG         = any_conserved_frame(atgOffets)
    consFrameNearCognate = any_conserved_frame(nearCognateOffsets)
    consFrameSameStart   = any_conserved_frame(sameStartOffsets)
    consFrameStop1       = any_conserved_frame(stop1offsets)

    result = set()
    if firstORF or nTermExt :
        if hasATG :               result.add('ATG')
        if hasNearCognate :       result.add('NearCognate')
        if hasSameStart :         result.add('SameStart')
        if consFrameATG :         result.add('ATGFrame')
        if consFrameNearCognate : result.add('NearCognateFrame')
        if consFrameSameStart :   result.add('SameStartFrame')
    if secondORF :
        if hasStop1           : result.add('Stop1')
        if hasSameStop1       : result.add('SameStop1')
        if consFrameStop1     : result.add('Stop1Frame')
    if firstORF or secondORF :
        if hasStopAtEnd :   result.add('StopAtEnd')

    if hasStopAtEnd :
        if firstORF :
            if hasATG and consFrameATG :                 result.add('ORF')
            if hasNearCognate and consFrameNearCognate : result.add('NearCognateORF')
            if hasSameStart and consFrameSameStart :     result.add('SameStartORF')
        if secondORF :
            if hasStop1 and consFrameStop1 :             result.add('SecondORF')
            if hasSameStop1 and consFrameStop1 :         result.add('SameStop1SecondORF')
    if nTermExt :
        if hasATG and consFrameATG :
            result.add('ATGExt')
        if hasNearCognate and consFrameNearCognate :
            result.add('NearCognateExt')
        if hasSameStart and consFrameSameStart :
            result.add('SameStartExt')

    assert all(ct in ConsTypes for ct in result), result
    assert all(ct in (ConsTypesForFirstORF if firstORF else
                      ConsTypesForSecondORF if secondORF else
                      ConsTypesForNTermExt) for ct in result), result
    return result

def _test_calc_seq_conservations() :
    for testArgs in [
        # firstORF
        ('ATG', 'TAG', 'CCC', 'ATG', False, False, None, # Good ORF
         ConsTypesForFirstORF),
        ('ATT', 'TAG', 'CCC', 'ATG', False, False, None, # Near cognate start
         ['NearCognate', 'StopAtEnd', 'ATGFrame', 'NearCognateFrame', 'SameStartFrame',
          'NearCognateORF']),
        ('ATT', 'TAG', 'CCC', 'ATT', False, False, None, # Same near cognate start
         ['NearCognate', 'SameStart', 'StopAtEnd', 'ATGFrame',
          'NearCognateFrame', 'SameStartFrame', 'NearCognateORF', 'SameStartORF']),
        ('.TG', 'TAG', 'CCC', 'ATT', False, False, None, # dot at start of near cognate
         ['NearCognate', 'StopAtEnd', 'NearCognateFrame', 'NearCognateORF']),
        ('|TG', 'TAG', 'CCC', 'ATT', False, False, None, # jump at start of near cognate
         ['NearCognate', 'StopAtEnd','NearCognateFrame', 'NearCognateORF']),
        ('ATG', 'TAC', 'CCC', 'ATG', False, False, None, # No stop
         ['ATG', 'NearCognate', 'SameStart',
          'ATGFrame', 'NearCognateFrame', 'SameStartFrame']),
        ('ATG', 'TAG', 'CCCA', 'ATG', False, False, None, # Not multiple of 3
         ['ATG', 'NearCognate', 'SameStart', 'StopAtEnd']),
        ('ATG', 'TAG', 'CCCAAA', 'ATG', False, False, None, # Good ORF
         ConsTypesForFirstORF),
        ('ATG', 'TAG', 'CCC.AA', 'ATG', False, False, None, # Dot
         ['ATG', 'NearCognate', 'SameStart', 'StopAtEnd']),
        ('ATG', 'TAG', 'CCC|AA', 'ATG', False, False, None, # Jump
         ['ATG', 'NearCognate', 'SameStart', 'StopAtEnd']),
        ('ATG', 'TAG', 'CCCTGA', 'ATG', False, False, None, # Interior in-frame stop
         ['ATG', 'NearCognate', 'SameStart', 'StopAtEnd']),
        ('ATG', 'TAG', 'CCTGAC', 'ATG', False, False, None, # Interior out-of-frame stop
         ConsTypesForFirstORF),
        ('ATA', 'TAG', 'CCCAAA', 'ATG', False, False, 'yeast', # Yeast mito start
         ['ATG', 'NearCognate', 'StopAtEnd',
          'ATGFrame', 'NearCognateFrame', 'SameStartFrame',
          'ORF', 'NearCognateORF']),
        ('ATA', 'TGA', 'CCCAAA', 'ATG', False, False, 'yeast', # Yeast mito (not stop
         ['ATG', 'NearCognate', 'ATGFrame', 'NearCognateFrame', 'SameStartFrame']),
        ('ATG', 'ATAG', 'CCC', 'ATG', False, False, None, # extra at start of last codon
         ['ATG', 'NearCognate', 'SameStart', 'StopAtEnd']),
        ('ATG', 'TAGA', 'CCC', 'ATG', False, False, None,  # extra at end of last codon
         ConsTypesForFirstORF),
        ('AATG', 'TAG', 'CCC', 'ATG', False, False, None,  # Extra at start of 1st codon
         ConsTypesForFirstORF),
        ('ATGA', 'TAG', 'CCC', 'ATG', False, False, None,  # Extra at end of 1st codon
         ['ATG', 'NearCognate', 'SameStart', 'StopAtEnd']),
        ('AATT', 'TAG', 'CCC', 'ATG', False, False, None,  # Near cognate, extra at start
         ['NearCognate', 'StopAtEnd', 'NearCognateFrame', 'NearCognateORF']),
        ('ATGCATT', 'TAG', 'CCC', 'ATT', False, False, None,  # out-of-frame ATG with
         ['ATG', 'NearCognate', 'SameStart', 'StopAtEnd',     # in-frame same near cognate
          'NearCognateFrame', 'SameStartFrame', 'NearCognateORF', 'SameStartORF']),
        ('ATGCATT', 'TAG', 'CCCGG', 'ATT', False, False, None, # in-frame ATG with
         ['ATG', 'NearCognate', 'SameStart', 'StopAtEnd',      # out-of-frame near cognate
          'ATGFrame', 'NearCognateFrame', 'ORF', 'NearCognateORF']),
        ('ATG', 'TAGCTGA', 'CCC', 'ATG', False, False, None,  # Two stops either frame
         ConsTypesForFirstORF),
        ('ATG', 'TAGCTGA', 'CCCG', 'ATG', False, False, None,  # Two stops not 3rd frame
         ['ATG', 'NearCognate', 'SameStart', 'StopAtEnd']),
        ('ATG', 'TAGCTGA', 'CCCGG', 'ATG', False, False, None,  # Two stops either frame
         ConsTypesForFirstORF),
        ('A-TG-', 'TAG', 'CCC', 'AT--G', False, False, None, # Dashes in 1st codon
         ConsTypesForFirstORF),
        ('ATG', 'TAG', 'C-CC', 'ATG', False, False, None,    # Dashes in mid-bases
         ConsTypesForFirstORF),
        ('ATG', '-T-AG', 'CCC', 'ATG', False, False, None,   # Dashes in last codon
         ConsTypesForFirstORF),

        # secondORF
        ('TGA', 'TAG', 'CCCAAA', 'TGA', True, False, None, # Good 2nd ORF with same stop1
         ConsTypesForSecondORF),
        ('TGA', 'TAG', 'CCCAAA', 'TAA', True, False, None, # Good 2nd ORF with diff stop1
         ['Stop1', 'StopAtEnd', 'Stop1Frame', 'SecondORF']),
        ('ATG', 'TAG', 'CCCAAA', 'ATG', True, False, None, # Same non-stop in 1st codon
         ['StopAtEnd', 'Stop1Frame']),
        ('ATGC', 'TAG', 'CCCAAA', 'TGA', True, False, None, # Frame for no stop1 is start
         ['StopAtEnd']),
        ('TGA', 'TAC', 'CCCAAA', 'TGA', True, False, None,  # No 2nd stop
         ['Stop1', 'SameStop1', 'Stop1Frame']),
        ('TGA', 'ATAG', 'CCCAAA', 'TGA', True, False, None,  # Out of frame 2nd stop
         ['Stop1', 'SameStop1', 'StopAtEnd']),

        # nTermExt
        ('ATG', 'ATG', 'CCC', 'ATG', False, True, None,  # Good ATG-initiated nTermExt
         ConsTypesForNTermExt),
        ('ATG', 'TAC', 'CCC', 'ATG', False, True, None,  # Good ATG-initiated nTermExt
         ConsTypesForNTermExt),                          # Final codon needn't be ATG
        ('ATT', 'TAC', 'CCC', 'ATT', False, True, None, # Good same near cognate nTermExt
         ['NearCognate', 'SameStart', 'ATGFrame', 'NearCognateFrame', 'SameStartFrame',
          'SameStartExt', 'NearCognateExt']),
        ('CTG', 'TAC', 'CCC', 'ATT', False, True, None, # Good diff near cognate nTermExt
         ['NearCognate', 'ATGFrame', 'NearCognateFrame', 'SameStartFrame',
          'NearCognateExt']),
        ('ATT', 'TAG', 'CCC', 'ATT', False, True, None, # Stop in final codon breaks frame
         ['NearCognate', 'SameStart']),
        ('ATT', 'GAGA', 'CCC', 'ATT', False, True, None,  # 4 base final codon breaks frame
         ['NearCognate', 'SameStart']),
        ('ATT', 'GA.', 'CCC', 'ATT', False, True, None,  # . in final codon breaks frame
         ['NearCognate', 'SameStart']),
        ('ATT', 'GA|', 'CCC', 'ATT', False, True, None,  # | in final codon breaks frame
         ['NearCognate', 'SameStart']),
        ('ATC', 'GATGCC', 'CCC', 'ATT', False, True, None,  # ATG out-of-frame in final
         ['NearCognate', 'ATGFrame', 'NearCognateFrame',    # codon does not break frame
          'SameStartFrame', 'NearCognateExt']),
    ] :
        def test(firstCodon, lastCodon, midBases, refStartCodon, secondORF, nTermExt,
                 mitoGroup, expectedResult) :
            bases = firstCodon + midBases + lastCodon
            firstCodonLastInd = len(firstCodon) - 1
            lastCodonFirstInd = len(firstCodon) + len(midBases)
            args = (bases, firstCodonLastInd, lastCodonFirstInd, refStartCodon,
                    secondORF, nTermExt, mitoGroup)
            assert OrderConsTypes(_calc_seq_conservations(*args)) == expectedResult, (
                firstCodon, lastCodon, midBases, refStartCodon, secondORF, nTermExt,
                mitoGroup, expectedResult, OrderConsTypes(_calc_seq_conservations(*args)))
        test(*testArgs)

    # Make sure it disallows secondORF and nTermExt both True
    try :
        _calc_seq_conservations('ATGCCCAAATAG', 2, 9, 'ATG', True, True, None)
    except AssertionError :
        pass
    else :
        raise NotImplementedError('Should not have allowed both secondORF and nTermExt.')
