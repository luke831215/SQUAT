# Usage

## For impatient users,
	./run.sh seq1 seq2 ... seqN
			 -o <output_dir>  
			 -r <ref_seq> 

Map sequences to assemblies

Map sequences to reference genomes

## Output
**[seq/] [seq.html]**

SQUAT will generate an HTML report named after the sequencing dataset and a directory containing all the analysis information.

## Command Options
SQUAT runs from a command line with the following options:  

**-o < path >**
> The path to output directory.

**-r < str >**
>Path to the assembly file as reference for alignment. SQUAT accepts FASTA format.

**-g < path >**
>Assembly file. The tool accepts assemblies or reference genomes in FASTA format.

**-t < int >** (or --thread < int >)  
>Number of threads to use. The default value is 1/3 of the number of CPUs of the current machine.

**-k** (or --keep)  
>Keep the sam file after each mapping experiment in `{output_dir}/{seq}/{aligner}`

**-s < str >**
>Return the subset of sequencing reads with labels specified in capitals. For ex., **-s PSC** means only selecting reads labeled with P, S, and C). The subset of sequencing reads will be stored in `{output_dir}/{seq}/subset`
    
**-g < str >**
> Path to the reference genome file for GAGE benchmark tool

**--gage**
>Activate gage mode for assembly evaluation, must specify reference genome (-g)