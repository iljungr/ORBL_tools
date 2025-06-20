# ORBL_tools

## Overview

ORBL is a tool for measuring evolutionary conservation and constraint on the "ORFness"
of an open reading frame (ORF) in a clade using multi-species whole genome alignments.
It consists of two scores:
- ORBLv measures conservation of the ORF by calculating the relative branch length of 
  the phylogenetic tree of species in an alignment that have an intact orthologous ORF, 
  i.e., in which there is an aligned start codon, stop codon, and open reading frame. 
  It is a number between 0 and 1, with larger numbers indicating more conservation.


- ORBLq measures evolutionary constraint on the ORFness of a non-canonical ORF (ncORF)
  by calculating the quantile of its ORBLv score among the ORBLv scores of untranslated 
  ORFs of the same biotype and similar length. It too is a number between 0 and 1, with 
  larger numbers indicating more constraint.

## Installation

```
git clone https://github.com/iljungr/ORBL_tools.git
```

Test by running:
```
cd ORBL_tools
python -m unittest test_ORBL
```

ORBL has been tested in python 2.7.16 and python 3.12.4 but should work using other versions of
python 2.7 and 3.

Requires an internet connection for downloading alignments.

## Usage Summary

```
ORBL_tools/orbl.py (ALIGNMENT_SET [--orblq] [--components] [FILE_NAME] | -h,--help | -v,--version)

Mandatory arguments:
    ALIGNMENT_SET Specifies the multispecies alignment as defined by CodAlignView.

Optional arguments:
    -h, --help     Print this message and exit.
    -v, --version  Print the orbl.py version number and exit.
    --orblq        Report the ORBLq constraint score as well as the orblv conservation 
                   score. This requires specifying the biotype-with-frameshift as the 
                   third field on each input line. Currently, --orblq is only implemented 
                   for certain alignment sets:
    --components   Report the relative branch lengths of the species satisfying each
                   or the three conditions used to consider the ORF to be conserved when
                   calculating the ORBLv score, namely having an aligned start codon, 
                   an aligned stop codon, and having an open reading frame (multiple of
                   3 nucleotides and no in-frame stop codons).
    FILE_NAME      Read input lines from FILE_NAME instead of standard input.
```

## Details

orbl.py takes input from the standard input or a specified file. Input consists of one 
or more lines, each representing an ORF in the reference species of a multispecies whole 
genome alignment, specified by the ALIGNMENT_SET mandory argument. Alignment sets are
defined by CodAlignView [here](https://data.broadinstitute.org/compbio1/cav.php?Alnsets). 

Each line contains two or more tab-separates fields. The first field is one or nore
chromosomal intervals specifying the coordinates of an open reading frame in the 
reference species of the alignment. The second field is the strand, either + or -. 
The format of the first field consists of one or more intervals separated by plus signs: 
```
chrom:start1-end1+chrom:start2-end2+... 
```
satisfying the requirement that:
```
start1 <= end1 < start2 <= end2... 
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

If the region specified does not begin with ATG or end in a stop codon in the reference 
species then both ORBLv and ORBLq will be "NA".

If --orblq is specified, each input line must inlude a third field, containing the biotype 
and frameshift relative to the main ORF for biotypes uoORF, intORF, and doORF. Frameshift
is +1 or +2 depending on whether a ribosome reading the main frame would need to skip 
1 or 2 nucleotides to get in the frame of the ncORF. Valid values of
biotype-with-frameshift are:
```
   uORF, uoORF+1, uoORF+2, intORF+1, intORF+2, doORF+1, doORF+2, dORF, and lncRNA-ORF.
```
Biotype "mixed" is also allowed, but ORBLq is not implemented for this biotype and 
results in a value of "NA".

Currently, --orblq is only implemented for the following alignment sets:
```
     hg38_120mammals_placental, hg38_120mammals_primate
```
Input lines may contain additional tab-separated fields, which are passed through to 
the output and otherwise ignored. 

Results are written to the standard output. There is one output line for each input line. 
Each output line contains the input line followed by the orblv score, the orblq score 
if the --orblq option is specified, and three additional fields if --components is 
specified, namely the relative branch lengths of species having an aligned start codon, 
an aligned stop codon, and an intact open reading frame.

## Credits

Questions should be directed to [Irwin Jungreis](mailto:iljungr@csail.mit.edu)

## Citing ORBL

If you use ORBL, please cite, "High-quality peptide evidence for 
annotating non-canonical open reading frames as human proteins" by Deutsch et al 
(manuscript submitted). 

More information about ORBL can be found in that paper. The specific requirements 
for determining an ORF's biotype are reported in the supplementary materials.
