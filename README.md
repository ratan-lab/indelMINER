indelMINER
==========

Identify indels from BAM files.

## REQUIREMENTS
indelMINER should work on any standard 64 bit Linux environment with 

- GCC
- zlib (http://zlib.net)

## ACKNOWLEDGEMENTS
indelMINER uses Paul Hsieh's Superfasthash (http://www.azillionmonkeys.com/qed/hash.html) hash function. The hash function is encapsulated in the source file superfasthash.c within this distribution.
indelMINER also includes the 0.1.19 version of SAMtools by Heng Li as part of
the distribution.

I would like to thank the following users for bug reports and other
contributions towards the on-going development of indelMINER.
- Trent Walradt, Yale School of Medicine
- Pieter Busschaert, University of Leuven
- Joao Rietra, Recife, Brazil for updates to the hash functions 

## SUMMARY
indelMINER refers to a set of algorithms to identify indels from whole genome 
resequencing datasets using paired-end reads. indelMINER uses a split-read 
approach to identify the precise breakpoints for indels of size less than a user
specified threshold, and supplements that with a paired-end approach to identify 
larger variants that are frequently missed with the split-read approach. 

## INSTALLATION
indelMINER uses the SAMtools (http://samtools.sourceforge.net) API to 
manipulate files in the BAM (Binary Alignment/Mapping) format. Current version
of indelMINER comes with the version 0.1.19 of SAMtools included as part of the
distribution.

Type:
```
make 
```
This should create an binary named indelminer in the "src" folder. Add the "src"
folder to your PATH or use the full path to the binary when running indelMINER

## TEST-DATASET
The test\_data directory contains the following files

- reference.fa : A reference sequence in FASTA format.
- alignments.bam : An alignment file in the BAM format.
- alignments.bai : Index for the provided alignment.bam.
- indelminer.config : An optional file that can be used to speed up a run.
- indelminer.expected.vcf : A list of expected variants in the BAM file in VCF format.

The alignment file in the BAM format for 100x100 bp paired end sequences with 
an mean insert length of 500 bps. The sequences were aligned to the reference 
sequence in reference.fa using BWA, then cleaned using GATK IndelRealigner. The 
indel alignments.bai was generated on a x86\_64 system, so please reindex if 
you are attempting to run on a different platform.

To run indelMINER on the provided BAM file, do the following (in the test\_data
folder)

```
./../src/indelminer reference.fa sample=alignments.bam > indelminer.flt.vcf
```

The expected output is provided in the file "indelminer.expected.vcf".  You can 
do a simple "diff" to compare the output in "indelminer.expected.vcf" to the 
output generated by your run in "indelminer.flt.vcf".

The first step in indelMINER involves calculation of the minimum and maximum
proper insert lengths for all the read-groups in the BAM, and the expected
average coverage on all the contigs/chromosomes in the reference sequence. The
user can speed up indelMINER by providing these information in a configuration
file to indelMINER. The configuration file should contain lines in the following
format:

```
IL read_group min_insert max_insert
```

to specify that a read group (read\_group) contains pairs with the minimum and
maximum proper insert lengths of min\_insert and max\_insert.

Similarly the chromosomal coverage can be specified using the following format:

```
RC chromosome_name average_coverage
```

There is a file called indelminer.config in your test\_data directory which
specifies these values for the BAM file included with the distribution. In order
to use it in a run of indelMINER, please do the following:

```
./../src/indelminer -i indelminer.config reference.fa sample=alignments.bam > indelminer.flt.vcf
```

## APPLICATION NOTES
### TUMOR/NORMAL PAIRS
In cases where you have a BAM file of alignments from a tumor sample, the first
step is to run indelMINER on the tumor sample.

```
indelminer reference.fa tumor=tumor.bam > tumor.vcf
```

This identifies all the indels in the tumor sample. You might want to change a
few of the arguments to indelMINER in this command based on yoru dataset.
Specially you might want to supply a configuration file, which might make
indelMINER run faster. You might want to change "-b" if you have smaller size
reads (The defaults work well for 100 bp reads).

One you have the indels in the tumor, run indelMINER to annotate those indels
based on their presence/absence in the normal.

```
indelminer -q 0 -a -e 1 reference.fa tumor.vcf normal=normal.bam > tumor.normal.vcf
```

This adds an annotation "normal" to each variant in tumor.vcf if it is also seen
in the normal. Using "-q 0 -a -e 1" ensures that any sign of its presence in the
tumor is considered.
