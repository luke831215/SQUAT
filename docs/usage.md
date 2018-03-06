# Usage

For impatient users,

	./run.sh -o <output_dir>  
			 -r <read_name> 
			 -i <fastq_path>  
			 -g <genome_path>

SQUAT runs from a command line with the following options:  

>**-o \<path>**
>>Path to the output directory.

>**-r \<str>**
>>Specify the name of the sequencing reads.

>**-i \<path>**
>>Sequencing dataset. The tool accepts sequencing reads in FASTQ format.

>**-g \<path>**
>>Reference genome file. The tool accepts assemblies and reference genomes in FASTA format.

>**-t \<int>** (or --thread <int>)  
>>Number of threads to use. The default value is 1/3 of the number of CPUs of the current machine.

>**-k** (or --keep)  
>>Keep the sam file after each mapping experiment in `{output_dir}/{aligner}/Pauto`

>**-s \<str>**
>>Return the subset of sequencing reads according to the labels (in capitals, e.g. PSCO). The subset in FASTA format will be stored in `{output_dir}/subset`



    
