# Output
---

Here we use the [example data](installation.md#example-data) in SQUAT to demonstrate the output structure in hierarchical fashion. Please execute the following commands first:

	cd SQUAT
	./squat.sh example/SEQ.fastq -o example -r example/ASSEMBLY.fasta

- Output directory: example
- Dataset: SEQ
- Assembly: ASSEMBLY

```bash
example/
│
├── SEQ.html //index page link to other report
│
├── SEQ.fastq //input sequencing reads
│
├── ASSEMBLY.fasta //input assembly
│
├── config //config file stating the threshold values
│
├── SEQ //directory containing reports and assessment result of the dataset SEQ
│   │
│   ├── SEQ.fastq //(sampled) input sequencing reads with modified id
│   │
│   ├── SEQ.ids //file to record the original and correspondent id   
│   │
│   ├── post-assembly_report.html //HTML report based on mapping reads to assemblies (read mapping)
│   │
│   ├── pre-assembly_report.htm //HTML report based on quality scores before genome assembly
│   │
│   ├── report.pdf //Post-assembly report in PDF version
│   │
│   ├── SEQ.log //log file of the whole SEQ quality assessment process
│   │
│   ├── bwa-mem //The mapping algorithms which performs local alignment
│   │   │
│   │   ├── align_info //A pickle-dumped file saving alignment fields of SAM file using Python
│   │   │
│   │   ├── ids
│   │   │   │
│   │   │   ├── SEQ_0_reads.cnt //Record the count of reads for each read label
│   │   │   │
│   │   │   ├── SEQ_0_reads.info //Record the label and repeat count of each read
│   │   │   │
│   │   │   ├── SEQ_1_mappable_unique_noerror.ids //Record read ids for type P reads
│   │   │   │
│   │   │   └── ...
│   │   │
│   │   ├── plot
│   │   │   │
│   │   │   ├── aln_score_C.png //Alignment score distribution graph of type C reads
│   │   │   │
│   │   │   ├── aln_score_P.png //Alignment score distribution graph of type P reads
│   │   │   │
│   │   │   ├── aln_score_S.png //Alignment score distribution graph of type S reads
│   │   │   │
│   │   │   ├── clip_ratio.png //Clip ratio distribution graph of type C reads
│   │   │   │
│   │   │   └── mismatch_ratio.png //Mismatch ratio distribution graph of type S reads
│   │   │
│   │   ├── index //Intermediate data and log file while indexing ASSEMBLY.fasta
│   │   │   │
│   │   │   └── ...
│   │   │
│   │   └── log // Log file of BWA-mem read mapping
│   │       │
│   │       └── ...
│   │
│   ├── bwa-backtrack //The mapping algorithms which adopts end-to-end strategy (similar structure as BWA-MEM)
│   │   ├── align_info 
│   │   ├── ids
│   │   │   └── ...
│   │   ├── plot
│   │   │   └── ...
│   │   ├── index
│   │   │   └── ...
│   │   └── log
│   │       └── ...
│   │
│   ├── images
│   │   │
│   │   ├── basic_stats.png //Basic statistics of sequencing reads (SEQ.fastq)
│   │   │
│   │   ├── eval_table.png //QUAST evaluation of the assembly (ASSEMBLY.fasta)
│   │   │
│   │   ├── label_dis_bar.png //Label distribution barchart
│   │   │
│   │   └── piechart.png //Label distribution piechart
│   │
│   ├── quast
│   │   │
│   │   ├── report.txt/report.pdf/report.html //Reports of genome assemblies evaluated by QUAST
│   │   │
│   │   ├── quast.log Log file of QUAST evaluation
│   │   │
│   │   └── ...
│   │
│   └── subset //The directory containing subsets of SEQ (exists if --subset specified)

```