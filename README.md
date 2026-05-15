# ORBL_tools

## Overview

ORF Relative Branch Length (ORBL) measures cross-species evolutionary conservation and
constraint on the "ORFness" of an open reading frame (ORF), without regard to
conservation of the encoded amino acid sequence. It is intended to detect ORFs encoding
poorly conserved peptides, as well as ORFs whose translation is functional but which do
not necessarily encode a functional peptide, such as regulatory uORFs. This is
particularly relevant for non-canonical ORFs (ncORFs). ORBL uses multispecies
whole-genome alignments to obtain the local alignment of the ORF in a particular clade,
and then computes two scores:

- ORBLv measures conservation of ORFness by calculating the relative branch length of the
  phylogenetic tree of species in the alignment that have an intact orthologous ORF,
  i.e., in which there is an aligned ATG start codon, stop codon, and open reading frame.
  It is a number between 0 and 1, with larger numbers indicating more conservation.


- ORBLq measures evolutionary constraint on ORFness by calculating the quantile of its
  ORBLv score among the ORBLv scores of untranslated ORFs of the same biotype and similar
  length. Comparison to these matched ORFs corrects for conservation due to chance or to
  constraint on an overlapping coding sequence (CDS). ORBLq is also a number between 0
  and 1, with larger numbers indicating more constraint. The number 1 - ORBLq can be
  interpreted as a p-value, since it approximates the probability that a similar ORF
  would get the same or higher ORBLv score under the null hypothesis that its ORFness
  were not constrained.

By detecting evolutionary constraint on translation itself, ORBL expands the 
scope of comparative genomics to detect functional ORFs that would be missed
by conventional protein conservation analyses.

## Installation

Install with the following line in a unix shell. It should only take a few seconds.
```
git clone https://github.com/iljungr/ORBL_tools.git
```
Test by running the following, which should take less than a minute and not report errors.
```
cd ORBL_tools
python -m unittest test_ORBL
```
ORBL has been tested in Python 2.7.16 and Python 3.12.4 but should work using 
other versions of Python 2.7 and 3.

Running orbl.py requires an internet connection for downloading alignments.


## Quick Example

Compute ORBLv for a single ORF:
```
echo "chr1:16395157-16395234    +" | python orbl.py hg38_120mammals_placental
```
Example output:
```
chr1:16395157-16395234    +    0.856779640335
```

## Usage Summary
```
python orbl.py (ALIGNMENT_SET [INPUT_FILE [OUTPUT_FILE]] [--orblq] [--components]
         [--CodAlignView] [--UCSCView] [--header] [--extraFields FIELD1,...]
         [--writeBed BED_FILE [--ORFNameField FIELD_NAME]] |
         -h,--help | -v,--version | --alignmentSets)

Mandatory arguments:
   ALIGNMENT_SET Specifies the multispecies whole-genome alignment. Use --alignmentSets 
                 for the list of available Alignment Sets. 
                 Examples: hg38_120mammals_placental, hg38_243primates.    

Optional arguments:
   INPUT_FILE      Read input lines from INPUT_FILE instead of standard input.
   OUTPUT_FILE     Write results to OUTPUT_FILE instead of standard output.
   --orblq         Report the ORBLq constraint score as well as the ORBLv conservation 
                   score. This requires specifying the biotype-with-frameshift as the 
                   third field on each input line. Currently, --orblq is only 
                   implemented for certain Alignment Sets.
   --components    Report the relative branch lengths of the species satisfying each
                   of the three conditions used to consider the ORF to be conserved when
                   calculating the ORBLv score, namely having an aligned start codon, 
                   an aligned stop codon, and an aligned open reading frame (multiple of
                   3 nucleotides and no in-frame stop codons).
   --CodAlignView  Add to each line a URL to show color-coded alignment in CodAlignView.
   --UCSCView      Add to each line a URL to show the interval containing the ORF in the
                   UCSC Genome Browser.
   --header        Preprend a header line describing the output columns.
   --extraFields   Comma-separated names of additional input fields that appear after
                   the required fields on each input line.
   --writeBed      Write a file representing all of the input ORFs in BED12 format.
   --ORFNameField  Use the specified input field to name each ORF in the BED file.
   -h, --help      Print this message and exit.
   -v, --version   Print the orbl.py version number and exit.
   --alignmentSets Open a browser window showing information about each allowed alignment
                   set, including whether ORBLq is implemented, and exit.
```
## Details

orbl.py takes input from the standard input or a specified file. Input consists of one 
or more lines, each representing an ORF in the reference species of a multispecies whole 
genome alignment, specified by the ALIGNMENT_SET mandatory argument. Alignment Sets include
all those defined by CodAlignView [here](https://data.broadinstitute.org/compbio1/cav.php?Alnsets), 
as well as some subclades of those. Use the --alignmentSets argument to see the complete 
list. Input is terminated by the end of the input file or, for interactive input, by an 
empty line. Input lines beginning with a '#' character are treated as comments; they are 
passed through to the output and otherwise ignored. For non-interactive input, blank 
lines are also passed through to the output and otherwise ignored.

Each line contains two or more tab-separated fields. The first field contains the 
coordinates of one or more chromosomal intervals comprising an open reading frame in the 
reference species of the alignment. The second field is the strand, either + or -. 
The format of the first field is one or more intervals separated by plus signs:
```
chrom:start1-end1+chrom:start2-end2+... 
```
Intervals must satisfy start <= end even if the region is on the minus strand. All
segments must be on the same chromosome. Coordinates follow the convention of GFF/GTF
files and the UCSC genome browser, namely, intervals go from the first base of the region
to the last base, inclusive, and coordinates are 1-based, i.e., the 1st base of the
chromosome is position 1. Coordinates may include commas, which are ignored.

See Examples section below for examples of input lines.

If the region specified is not an ORF in the reference species (i.e., first codon is
not ATG, last codon is not a stop codon, length is not a multiple of three nucleotides, 
or contains a premature in-frame stop codon) then both ORBLv and ORBLq will be "NA",
as well as the corresponding components if the --components option is specified.

If --orblq is specified, each input line must include a third field containing the ORF
biotype, and each output line will report the ORBLq score after the ORBLv score. Biotypes
for non-canonical ORFs are explained
[here](https://github.com/iljungr/ORBL_tools/blob/main/docs/ncORF_Biotypes.pdf). For
biotypes uoORF, intORF, and doORF, which indicate overlap with a known CDS, the
frameshift relative to the CDS must be appended. Frameshift is +1 or +2 depending on
whether a ribosome reading the main frame would need to skip 1 or 2 nucleotides to get in
the frame of the ORF. Valid values of biotype-with-frameshift are:
```
uORF, uoORF+1, uoORF+2, intORF+1, intORF+2, doORF+1, doORF+2, dORF, and lncRNA-ORF.
```
Biotype "mixed" is also allowed, but ORBLq is not implemented for this biotype and 
results in a value of "NA". If the frameshift for a uoORF, intORF, or doORF is 
unknown, we suggest computing ORBLq for both +1 and +2 (using two input lines), 
with further investigation of the frameshift only needed if one but not the other 
ORBLq value exceeds some threshold of interest.

Currently, ORBLq is only implemented for the human hg38/GRCh38 assembly, and for
only some of its Alignment Sets. Use the --alignmentSets option to see which ones.

Input lines may contain additional tab-separated fields, which are passed through to 
the output and otherwise ignored. All input lines must contain the same number of
additional fields.

Results are written to OUTPUT_FILE if specified, otherwise to standard output. There is
one output line for each input line. Each output line contains the input line followed by
the ORBLv score; the ORBLq score if the --orblq option is specified; three additional
fields if --components is specified, namely the relative branch lengths of species having
an aligned start codon, an aligned stop codon, and an intact open reading frame; a URL to
show the color-coded alignment of the ORF with 15 nt of context on each end if 
--CodAlignView is specified; and a URL to open the UCSC Genome Browser to the interval
containing the ORF if --UCSCView is specified.

If --header is specified, orbl.py prepends a header line describing the output columns,
making the output into a tab-delimited file that can be read into a spreadsheet program.
If each input line contains fields after the required ones, the --extraFields argument
followed by a comma-separated list of names can be used to name these fields; otherwise
default names will be assigned.

If --writeBed is specified, a file will be written containing the coordinates of all the
ORFs in BED12 format, which can be used to view the ORFs in a genome browser. If
--ORFNameField is supplied, it specifies a field in each input line that will be used
as the name of each ORF in the BED12 file; otherwise a default name will be assigned.

## Typical Workflow

1. Prepare a file containing the genomic coordinates of candidate ORFs.
2. Choose an Alignment Set for the desired reference species and clade. To see
   available Alignment Sets use --alignmentSets.
3. If ORBLq is implemented for specified Alignment Set (currently only for human 
   hg38/GRCh38), include ORF biotypes in input file. Include frameshift (+1 or +2)
   relative to main ORF for overlapping biotypes (uoORF, intORF, and doORF).
4. Run orbl.py. Include option --orblq if it is implemented. Include option
   --CodAlignView to get links for viewing color-coded alignment. 
5. Investigate high-scoring candidates using CodAlignView.

## Examples

Run these in the ExampleFiles directory. Expected output files are in 
ExampleFiles/ExpectedOutput.

### Calculate placental mammal ORBLv for two ORFs

Input file: SimpleExample.in
```
chr1:16395157-16395234	+
chr10:1042727-1042762+chr10:1043301-1043321	-
```

Calculate ORBLv conservation score for the placental mammal clade using the 116
placental mammal subset of the 120 mammals alignment:
```
../orbl.py hg38_120mammals_placental SimpleExample.in SimpleExample.116placentals.out
```

Expected output:
```
chr1:16395157-16395234	+	0.856779640335
chr10:1042727-1042762+chr10:1043301-1043321	-	0.8129781736
```

### Calculate primate ORBLq and add CodAlignView links

Input file, which includes biotype-with-frameshift needed for ORBLq plus two extra
fields and comments: ORBLqExample.in. 

```
## Note: two extra fields represent ORF name and number of codons
# A one-exon dORF:
chr1:16395157-16395234	+	dORF	c1riboseqorf29	26
# A two-exon intORF shifted +1 nt from main frame:
chr10:1042727-1042762+chr10:1043301-1043321	-	intORF+1	c10norep2	19
```
Calculate both ORBLv conservation score and ORBLq constraint score for primate clade
using the 243 primate alignment, and add CodAlignView link:
```
../orbl.py hg38_243primates --orblq --CodAlignView ORBLqExample.in ORBLqExample.243primates.out
```
Each output line has input line plus ORBLv score, ORBLq score, and
URL that can be pasted in browser address bar to view color-coded alignment
in CodAlignView.

Expected output:
```
## Note: two extra fields represent ORF name and number of codons
# A one-exon dORF:
chr1:16395157-16395234	+	dORF	c1riboseqorf29	26	0.994444706512	0.996297153437	https://data.broadinstitute.org/compbio1/cav.php?a=hg38_243primates&i=chr1:16395157-16395234&s=+&p=15&e=15&os
# A two-exon intORF shifted +1 nt from main frame:
chr10:1042727-1042762+chr10:1043301-1043321	-	intORF+1	c10norep2	19	0.833713772092	0.56202698256	https://data.broadinstitute.org/compbio1/cav.php?a=hg38_243primates&i=chr10:1042727-1042762+chr10:1043301-1043321&s=-&p=15&e=15&os
```
### Example with 58 mammals and most optional arguments

Input file: ORBLqExample.in (same as above)

Calculate scores relative to 58-placental-mammal subset of 100-vertebrates alignment.
Report ORBLv, ORBLq, and relative branch length of ATG, stop, and frame; add links to
show ORF in CodAlignView and to open UCSC Genome Browser to interval containing ORF; add
header line using user-specified names for extra fields; write a file in BED12 format
using specified field as ORF name.
```
../orbl.py hg38_58 ORBLqExample.in ORBLqExample.58mammals.WithHeader.tsv \
    --orblq --components --CodAlignView --UCSCView  \
    --header --extraFields Name,NumCodons \
    --writeBed ORBLqExample.bed --ORFNameField Name
```
See ORBLqExample.58mammals.WithHeader.tsv and ORBLqExample.bed in ExpectedOutput
directory for expected output.


## Usage Notes

This section discusses considerations when choosing which tool to use and what to expect.

ORBLv measures ORF conservation whereas ORBLq measures evolutionary constraint on
ORFness. ORF conservation can result from constraint on ORFness, indicating that
translation is functional, but it can also result from other constraints on the start,
stop, or reading frame unrelated to translation, or simply be due to chance. ORBLq is
intended to more directly test if the conservation is due to constraint on ORFness, so
ORBLq is preferable to ORBLv when the goal is to determine if translation of an ORF is
functional. However, ORBLq is currently implemented only for some Alignment Sets. Some
things to consider when using ORBLv for alignments for which ORBLq is not implemented are
that a high ORBLv score is more likely to arise by chance if the ORF is short or if the
phylogenetic branch length of the alignment is small, and may also reflect "free"
conservation from an overlapping CDS for intORFs, uoORFs, and doORFs.

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
multispecies alignment, such as PhyloP [1], PhyloCSF [2], dN/dS (for example, PAML
implementation [3]), and BLS [4], but they have different purposes, strengths, and
weaknesses. ORBLq is specifically measuring constraint on ORFness, indicating functional
translation, which can be due to production of a functional protein or to a regulatory
effect of translation itself, as is common among uORFs. In contrast, PhyloCSF and dN/dS
are testing for signatures of protein-coding regions and would not be expected to detect
ORFs whose translation has a regulatory effect but does not produce a functional protein.
PhyloP measures nucleotide-level constraint, which can be due to many functions other
than translation. BLS is conceptually similar to ORBLv, but uses a much looser definition
of ORF conservation; consequently it would be expected to have higher sensitivity but
lower specificity than ORBLv. Also BLS measures conservation and has no equivalent of
ORBLq to measure constraint.

## Credits

Questions should be directed to [Irwin Jungreis](mailto:iljungr@csail.mit.edu) and
[Manolis Kellis](mailto:manoli@mit.edu).

## Citing ORBL

If you use ORBL, please cite:
- **Deutsch et al.** Expanding the human proteome with microproteins and peptideins. 
*Nature* 2026 https://doi.org/10.1038/s41586-026-10459-x.
- **Jungreis, I. and Kellis, M.** ORBL_tools: tools for measuring evolutionary conservation 
  and constraint of 'ORFness' of an open reading frame. *Zenodo*. (2026).
  [![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.18749291.svg)](https://doi.org/10.5281/zenodo.18749291)
  (This DOI refers to all versions. For reproducibility, please cite the
specific version you used (e.g., v1.0), available by running "orbl.py --version".)

More information about ORBL can be found in Deutsch et al. The specific requirements
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