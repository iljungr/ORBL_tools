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
Utilties for working with intervals on DNA strands.
"""
from __future__ import division, print_function

def antistrand(strand) :
    return '-' if strand == '+' else '+' if strand == '-' else strand

def regionString_to_triples(regionStr) :
    """
    Parse a string of form chrom:START-END[+chrom:START-END]*,
        where the numbers can include commas.
    Return [(chrom, START, END)...]
    """
    triples = []
    for intervalStr in regionStr.split('+') :
        try :
            chrom, intervalPart = intervalStr.split(':')
            intervalPart = intervalPart.replace(',', '')
            uminus = intervalPart.startswith('-')
            if uminus :
                intervalPart = intervalPart[1 :]
            start, end = map(int, intervalPart.split('-', 1))  # Keeps unary minus on end
            if uminus :
                start = -start
            triples.append((chrom, start, end))
        except ValueError :
            raise ValueError('Invalid interval string: %s' % intervalStr)
    return triples

def get_intervals_length(intervals, includesChrom = True) :
    """
    Return number of bases in a list of intervals: [(CHROM, START, END)...]
        (or [(START, END)...] if not includesChrom)
    """
    startIndex = 1 if includesChrom else 0
    return sum(interval[startIndex + 1] - interval[startIndex] + 1
               for interval in intervals)

def intervals_prefix(intervals, strand, numBases, includesChrom = True) :
    """
    Return subset of intervals numBases long starting at beginning relative to strand.
    Intervals is a list: [(CHROM, START, END)]
        (or [(START, END)...] if not includesChrom)
    """
    if numBases < 0 or numBases > get_intervals_length(intervals, includesChrom) :
        raise ValueError('intervals_prefix: invalid numBases (%s) for length %s' %
                         (numBases, get_intervals_length(intervals, includesChrom)))
    result = []
    reverse = strand == '-'
    startInd = 1 if includesChrom else 0
    for interval in intervals[: :[1, -1][reverse]] :
        intervalLen = interval[startInd + 1] - interval[startInd] + 1
        if numBases >= intervalLen :
            result.insert([len(result), 0][reverse], interval)
            numBases -= intervalLen
            continue
        if numBases > 0 :
            if reverse :
                result.insert(0, interval[:startInd] +
                              (interval[startInd + 1] - numBases + 1,
                               interval[startInd + 1]))
            else :
                result.insert(len(result), interval[:startInd] +
                              (interval[startInd],
                               interval[startInd] + numBases - 1))
        return result
    return result

def intervals_suffix(intervals, strand, numBases, includesChrom = True) :
    """
    Return subset of intervals numBases long ending at end relative to strand.
    Intervals is a list: [(CHROM, START, END)]
        (or [(START, END)...] if not includesChrom)
    """
    return intervals_prefix(intervals, ['-', '+'][strand == '-'], numBases, includesChrom)

def remove_final_codon(intervals, strand) :
    """ Return intervals with final codon removed """
    intervals = list(intervals) # To avoid changing original.
    for ii in range(3) :
        if len(intervals) == 0 :
            break
        if strand == '-' :
            intervals[0] = [intervals[0][0] + 1, intervals[0][1]] # Don't change original
            if intervals[0][0] > intervals[0][1] :
                del intervals[0]
        else :
            intervals[-1] = [intervals[-1][0], intervals[-1][1] - 1] # Don't change orig
            if intervals[-1][1] < intervals[-1][0] :
                del intervals[-1]
    return intervals

def anti_interval(interval, strandOfInterval) :
    """
    Return interval on other strand that shares 3rd codon positions, and other strand.
    Note that it always shifts 2, which makes sense if this interval has length a
        multiple of three, or is part of a set of intervals that have overall length
        a multiple of 3, but does not make sense if this is an isolated interval that
        is not a multiple of 3.
    Explanation: if interval is coding, then the anti-interval often gets a "ghost"
        evolutionary signature of coding. Comparison might distinguish which is real.
    We could instead return an interval that shares the first two codon positions, i.e.,
        shift back 1 instead of forward 2. At first glance, that would make sense,
        because the coding evolutionary signature comes from suppression of substitutions
        in the first two positions rather than excess substitutions in the third position.
        However, when carefully considering the purpose of looking at the anti-interval,
        it seems like the right shift is slightly better. The goal is to distinguish
        the real signal from the ghost by comparing evolutionary signatures of the two. We
        want an antisense interval that will get a bad score if it is a ghost and a good
        one if it is real. In most cases it's a wash, but for the common case that the
        original interval ends just before a stop, there's a slight advantage to the
        right shift. That's because in that case the anti-interval will be enriched
        for stops if it is a ghost (because a TAA or TAG codon on the original strand
        becomes TAx on that anti-strand), whereas these will have been suppressed if the
        anti-interval is coding.
    """
    antiOffset = 2 if strandOfInterval == '+' else -2
    return (type(interval)([interval[0] + antiOffset, interval[1] + antiOffset]),
            antistrand(strandOfInterval))

def anti_intervals(intervals, strandOfIntervals) :
    """
    Return intervals on other strand that share 3rd codon positions, and other strand.
    Only works if total intervals length is a multiple of 3.
    See explanation in anti_interval.
    Results at splice boundaries are silly.
    """
    assert sum(e - s + 1 for s, e in intervals) % 3 == 0, intervals
    return (type(intervals)([anti_interval(interval, strandOfIntervals)[0]
                             for interval in intervals]),
            antistrand(strandOfIntervals))
