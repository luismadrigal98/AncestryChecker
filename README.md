# Ancestry Checker

This program will check the ancestry of an F2 individual and its correspondence with the expectations under a NAM design. If the genome of the founder line per family are known, the expectations are that the F2 is a mosaic of the founder genomes. The program will also make quality checking over the vcf, preserving only for the ancestry analysis SNPs and calls where we have data completeness for all founders and at least one actual sample. We defined the filtering in that regard because the samples could have been sequenced at low coverage. In that case, it is likely to have missing data in at least one of those samples per site. [A flag allow to consider no missing data at all in case of high coverage sequencing for the samples]

# Input

Input includes:

1) The name of the founders as they appear in the vcf file.
2) The name of the F2 samples as they appear in the vcf file.
3) Relationship map (txt), based on the crosses performed.

## Example of relationship map

```txt

17  F2  155c    664c
123 F2  444c    664c
657 F2  767c    664c

...

```

AS can be seen, the format is a plain txt file, with '\t' as separator, and the first column states the ID of the sample, the second one, the filial position in the crosses, the third and fourth columns states the founders of the family.
