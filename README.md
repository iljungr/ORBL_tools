# ORBL_tools

## Overview

ORBL measures cross-species evolutionary conservation and constraint on the "ORFness" of
an open reading frame (ORF), without regard to conservation of the encoded amino acid
sequence. It is intended to distinguish ORFs, such as regulatory uORFs, whose translation
is functional, but which do not necessarily encode a functional peptide. It uses
multi-species whole genome alignments to obtain the local alignment of the ORF in a
particular clade, and then computes two scores:

- ORBLv measures conservation of ORFness by calculating the relative branch length of the
  phylogenetic tree of species in the alignment that have an intact orthologous ORF,
  i.e., in which there is an aligned ATG start codon, stop codon, and open reading frame.
  It is a number between 0 and 1, with larger numbers indicating more conservation.

- ORBLq measures evolutionary constraint on ORFness by calculating the quantile of its
  ORBLv score among the ORBLv scores of untranslated ORFs of the same biotype and similar
  length. Comparison to these matched ORFs corrects for conservation due to chance or to
  constraint on an overlapping CDS. ORBLq is also a number between 0 and 1, with larger
  numbers indicating more constraint. The number 1 - ORBLq can be thought of as a
  p-value, since it approximates the probability that a similar ORF would get the same or
  higher ORBLv score under the null hypothesis that its ORFness were not constrained.

## Installation

Install with the following line on a unix shell. It should only take a few seconds.
```
git clone https://github.com/iljungr/ORBL_tools.git
```
Test by running the following, which should take less than a minute and not report errors.
```
cd ORBL_tools
python -m unittest test_ORBL
```
ORBL has been tested in python 2.7.16 and python 3.12.4 but should work using 
other versions of python 2.7 and 3.

Running orbl.py requires an internet connection for downloading alignments.

## Usage Summary
```
python orbl.py (ALIGNMENT_SET [--orblq] [--components] [FILE] |
                -h,--help | -v,--version | --alignmentSets)

Mandatory arguments:
   ALIGNMENT_SET Specifies the multi-species whole-genome alignment. Use --alignmentSets 
                 for the list of available alignment sets. 
                 Example: hg38_120mammals_primate.    

Optional arguments:
   -h, --help      Print this message and exit.
   -v, --version   Print the orbl.py version number and exit.
   --alignmentSets Open a browser window showing information about each allowed alignment
                   set including whether ORBLq is implemented.
   --orblq         Report the ORBLq constraint score as well as the orblv conservation 
                   score. This requires specifying the biotype-with-frameshift as the 
                   third field on each input line. Currently, --orblq is only 
                   implemented for certain alignment sets.
   --components    Report the relative branch lengths of the species satisfying each
                   or the three conditions used to consider the ORF to be conserved when
                   calculating the ORBLv score, namely having an aligned start codon, 
                   an aligned stop codon, and having an open reading frame (multiple of
                   3 nucleotides and no in-frame stop codons).
   FILE            Read input lines from FILE instead of standard input.
```
## Details

orbl.py takes input from the standard input or a specified file. Input consists of one 
or more lines, each representing an ORF in the reference species of a multi-species whole 
genome alignment, specified by the ALIGNMENT_SET mandory argument. Alignment sets include
all those defined by CodAlignView [here](https://data.broadinstitute.org/compbio1/cav.php?Alnsets), 
as well as some subclades of those. Use the --alignmentSets argument to see the complete 
list. Input is terminated by the end of the input file, or, for interactive input, by an 
empty line. Input lines beginning with a '#' character are treated as comments; they are 
passed through to the output and otherwise ignored.

Each line contains two or more tab-separated fields. The first field is one or more
chromosomal intervals specifying the coordinates of an open reading frame in the 
reference species of the alignment. The second field is the strand, either + or -. 
The format of the first field consists of one or more intervals separated by plus signs: 
```
chrom:start1-end1+chrom:start2-end2+... 
```
satisfying the requirement that:
```
start1 <= end1, start2 <= end2, ... 
```
(even if the region is on the minus strand). 

All segments must be on the same chromosome.

Example input line:
```
chr10:1042727-1042762+chr10:1043301-1043321 +
```
Coordinates follow the convention of GFF/GTF files and the UCSC genome browser, namely, 
intervals go from the first base of the region to the last base, inclusive, and 
coordinates are 1-based, i.e., the 1st base of the chromosome is position 1. Coordinates
may include commas, which are ignored. 

If the region specified is not an ORF in the reference species (i.e., first codon is
not ATG, last codon is not a stop codon, length is not a multiple of 3 nucleotides, 
or contains a premature in-frame stop codon) then both ORBLv and ORBLq will be "NA",
as well as the corresponding components if the --components option is specified.

If --orblq is specified, each input line must include a third field, containing the
biotype and frameshift relative to the main ORF for biotypes uoORF, intORF, and doORF. 
Frameshift is +1 or +2 depending on whether a ribosome reading the main frame would 
need to skip 1 or 2 nucleotides to get in the frame of the ORF. Valid values of
biotype-with-frameshift are:
```
uORF, uoORF+1, uoORF+2, intORF+1, intORF+2, doORF+1, doORF+2, dORF, and lncRNA-ORF.
```
Biotype "mixed" is also allowed, but ORBLq is not implemented for this biotype and 
results in a value of "NA". If the frameshift for a uoORF, intORF, or doORF is 
unknown then we suggest computing ORBLq for both +1 and +2 (using two input lines), 
with further investigation of the frameshift only needed if one but not the other 
ORBLq value exceeds some threshold of interest.

Currently, --orblq is only implemented for the human hg38/GRCh38 assembly, and for
only some of its alignment sets. Use the --alignmentSets option to see which ones.

Input lines may contain additional tab-separated fields, which are passed through to 
the output and otherwise ignored. 

Results are written to the standard output. There is one output line for each input line. 
Each output line contains the input line followed by the orblv score, the ORBLq score 
if the --orblq option is specified, and three additional fields if --components is 
specified, namely the relative branch lengths of species having an aligned start codon, 
an aligned stop codon, and an intact open reading frame.

## Usage Notes

This section discusses considerations when choosing which tool to use and what to expect.

ORBLv measures ORF conservation whereas ORBLq measures evolutionary constraint on
ORFness. ORF conservation can result from constraint on ORFness, indicating that
translation is functional, but it can also result from other constraints on the start,
stop, or reading frame unrelated to translation, or simply be due to chance. ORBLq is
intended to more directly test if the conservation is due to constraint on ORFness, so
ORBLq is preferrable to ORBLv when the goal is to determine if translation of an ORF is
functional. However, ORBLq is currently implemented only for some Alignment Sets. Some
things to consider when using ORBLv for alignments for which ORBLq is not implemented are
that a high ORBLv score is more likely to be due to chance if the ORF is short or if the
phylogenetic branch length of the alignment is small, and is more likely to be due to
"free" conservation from an overlapping CDS for intORFs, uoORFs, and doORFs.

ORBL uses a strict definition of ORF conservation, requiring aligned start and stop
codons and full length open reading frame, whereas a functional ORF might still be
functional if the start and stop codons move a little. This strict definition gives ORBLq
high specificity but might lower the sensitivity. The main consideration that can lower
specificity is that an ORF can be conserved due to constraint on features other than
ORFness. Although ORBLq considers and adjusts for overlap with a CDS, there are many
functional elements it does not adjust for, such as enhancers. To minimize ORBLq false
positives, we recommend using CodAlignView to investigate the alignment of any candidate
functional ORF; if conservation extends to the regions flanking the ORF then it is
probably due to something other than constraint on ORFness.

Two things to consider when choosing which clade to use are that ORBLq scores with
respect to smaller clades might be able to capture lineage-specific constraint, but these
clades do not provide the high statistical power needed for ORBLq to have high
sensitivity. For example, Hominoidea (apes) have so little phylogenetic branch length
that it is easy for untranslated ORFs to be almost perfectly conserved due to chance, and
some will have equal or higher ORBLv scores than even well conserved functional ORFs,
lowering the ORBLq scores of the latter.

There are many tools available for calculating conservation or constraint using a
multi-species alignment, such as PhyloP [1], PhyloCSF [2], dN/dS (for example, PAML
implementation [3]), and BLS [4], but they have different purposes, strengths, and
weaknesses. ORBLq is specifically measuring constraint on ORFness, indicating functional
translation, which can be due to production of a functional protein or to a regulatory
effect of translation itself, as is common among uORFs. In contrast, PhyloCSF and dN/dS
are testing for signatures of protein-coding regions and would not be expected to detect
ORFs whose translation has a regulatory effect but does not produce a functional protein.
PhyloP measures nucleotide-level constraint, which can be due to many functions other
than translation. BLS is conceptually similar to ORBLv but uses a much looser definition
of ORF conservation; consequently it would be expected to have higher sensitivity but
lower specificity than ORBLv. Also BLS measures conservation and has no equivalent of
ORBLq to measure constraint.

## Credits

Questions should be directed to [Irwin Jungreis](mailto:iljungr@csail.mit.edu)

## Citing ORBL

If you use ORBL, please cite, "High-quality peptide evidence for 
annotating non-canonical open reading frames as human proteins" by Deutsch et al 
(manuscript submitted). 

More information about ORBL can be found in that paper. The specific requirements 
for determining an ORF's biotype are reported in the supplementary materials.

## References

[1] **Pollard KS et al.** Detection of nonneutral substitution rates on mammalian
    phylogenies. *Genome Research*. 2010.  
    https://doi.org/10.1101/gr.097857.109  

[2] **Lin MF et al.** PhyloCSF: a comparative genomics method to distinguish protein
    coding and non-coding regions. *Bioinformatics*. 2011.  
    https://doi.org/10.1093/bioinformatics/btr209  

[3] **Yang Z.** PAML 4: Phylogenetic analysis by maximum likelihood. *Molecular
    Biology and Evolution*. 2007.  
    https://doi.org/10.1093/molbev/msm088  

[4] **Chang et al.** Evolutionary remodeling of non-canonical ORF translation
    in mammals. *eLife (reviewed preprint)*. 
    https://elifesciences.org/reviewed-preprints/109128  