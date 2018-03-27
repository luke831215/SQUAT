# Read Labeling

<br>

## Summary

In this section, we go through each process of classifying sequencing reads into seven different categories. 

<br>

## BWA Mapping

SQUAT maps sequencing reads against the reference assemblies using two different algorithms: **BWA-MEM and BWA-backtrack**. The former algorithm performs local alignment with multiple alignments for different part of a query sequence, while the latter is more suitable for short reads and tries to map the whole sequence. For more explanation, see the [manual of BWA](http://bio-bwa.sourceforge.net/bwa.shtml).

<br>

## Reading Labeling Workflow

The basic workflow looks as following,

![read labeling workflow](imgs/label_workflow.png)

First, we label the sequences containing N to the N class. Next, based on the Flag and MAP fields in SAM file, reads indicated as unmapped and multi-mapped are assigned to F and M class respectively. The rest of reads are uniquely-mapped. They an be further divided into P(Perfect), S(substitution errors), C(contain clips), and O(other errors) class according to the CIGAR field of their primary alignment.

| Label | Description|
|:-:|:-:|
| P | Perfectly-matched reads |
| S | Reads with substitution errors |
| C | Reads that contain clips |
| O | Reads with other errors|
| M | Multi-mapped reads |
| F | Unmapped reads |
| N | Reads that contain N |

For more details of read labeling, see [SAM Spec](http://samtools.github.io/hts-specs/SAMv1.pdf).

<br>

## Label Distribution

### Table

For each strategy SQUAT adopts for alignment, it lists the label distribution with an icon to represent each read class in a table.

![Label Distribution Table](imgs/label_dis_table.png)

<br>

### Barchart

In addition, we plot a barchart to display the **percentage of high and poor quality reads** for each read class. Reads in P and M class are always considered high quality while unmapped (F) reads will always be of poor quality. For reads in S, C, O and N class, we compute certain ratio and specify a threshold to determine the sequencing quality of each read. For computation of the ratio, see [here](ratio.md)

![Label Distribution Barchart](imgs/label_dis_bar.png)
