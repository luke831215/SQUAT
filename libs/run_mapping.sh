#!/bin/bash

#Automatic exit from bash shell script on error
set -e

while getopts ":h" opt; do
  case $opt in
    h)
      echo -e "arg1: source directory"
      echo -e "arg2: output directory"
      echo -e "arg3: name of the data"
      echo -e "arg4: number of reads"
      echo -e "arg5: path of reference genome"
      echo -e "arg6: path of reads fastq file" 
      echo -e "arg7: number of threads" 
      echo -e "arg8(optional): bwaOnly or bowtie2Only\n" >&2
      ;;
    :)
      echo "Option -$OPTARG requires an argument." >&2
      exit 1
      ;;
  esac
done

if [ "$#" -lt "7" ]; then
	echo "Using:" $0 " [SRCDIR] [OUTDIR] [dataname] [READSIZE] [REFLOC] [ECVLOC] [MAXPROCESS] [bwaOnly|bowtie2Only]"
	exit 1
fi

function to_abs	{
	case $1 in
  		/*) absolute=$1;;
  		*) absolute=$PWD/$1;;
	esac

	echo $absolute
}

SRCDIR=$1
EXECDIR="$( cd "$(dirname "$0")" ; pwd)"
SCRIPTDIR=${EXECDIR}/map_proc

ORGDIR=${PWD}
OUTDIR=$2
DATA=$3
READSIZE=$4
REFLOC=$5
ECVLOC=$6
MAXPROCESS=$7

#mkdir -p ${OUTDIR}/imgs

#if [[ $8 != "bwaOnly" ]]; then

	#bowtie2 - local alignment
	#The log output of Bowtie is directed to stderr
	#echo 'bowtie2 - local alignment'
	#BOWTIEDIR=${OUTDIR}/bowtie2-local
	#BOWTIEDIR1=${BOWTIEDIR}
	#if [ -d ${BOWTIEDIR} ]; then
	#	rm -rf ${BOWTIEDIR}
	#fi
	#mkdir -p ${BOWTIEDIR}/index &> /dev/null
	#cd ${BOWTIEDIR};
	#ln -s ${ECVLOC} ${DATA}.fastq 2>&1 > /dev/null
	#ln -s ${REFLOC} scaffolds.fasta 2>&1 > /dev/null
	#${SRCDIR}/bowtie2/bowtie2-build scaffolds.fasta index/${DATA} > build.log 2>&1
	#${SRCDIR}/bowtie2/bowtie2 -p ${MAXPROCESS} --local -x index/${DATA} -U ${ECVLOC} -S ${DATA}_ecv_all.sam -k 5 --omit-sec-seq --reorder 2> ${DATA}_ecv_all.log
	#bash ${SCRIPTDIR}/bowtie2_benchmark_v2.sh ${DATA} scaffolds auto ${SRCDIR} I 2>&1 > /dev/null
	#rm ${DATA}.fastq scaffolds.fasta &> /dev/null
#
	##bowtie2 - end to end
	##The log output of Bowtie is directed to stderr
	#echo 'bowtie2 - end to end'
	#BOWTIEDIR=${OUTDIR}/bowtie2-endtoend
	#BOWTIEDIR2=${BOWTIEDIR}
	#if [ -d  ${BOWTIEDIR} ]; then
	#	rm -rf ${BOWTIEDIR}
	#fi
	#mkdir -p ${BOWTIEDIR}/index &> /dev/null
	#cd ${BOWTIEDIR};
	#ln -s ${ECVLOC} ${DATA}.fastq 2>&1 > /dev/null
	#ln -s ${REFLOC} scaffolds.fasta 2>&1 > /dev/null
	#${SRCDIR}/bowtie2/bowtie2-build scaffolds.fasta index/${DATA} > build.log 2>&1
	#${SRCDIR}/bowtie2/bowtie2 -p ${MAXPROCESS} -x index/${DATA} -U ${ECVLOC} -S ${DATA}_ecv_all.sam -k 5 --omit-sec-seq --reorder 2> ${DATA}_ecv_all.log
	#bash ${SCRIPTDIR}/bowtie2_benchmark_v2.sh ${DATA} scaffolds auto ${SRCDIR} I 2>&1 > /dev/null
	#rm ${DATA}.fastq scaffolds.fasta &> /dev/null
#fi

if [[ $8 != "bowtie2Only" ]]; then
	#bwa-mem
	echo 'BWA - MEM'
	BWADIR=${OUTDIR}/bwa-mem
	BWADIR1=${BWADIR}
	if [ -d  ${BWADIR} ]; then
		rm -rf ${BWADIR} &> /dev/null
	fi
	mkdir -p ${BWADIR}/index &> /dev/null
	cd ${BWADIR};
	ln -s ${ECVLOC} . 2>&1 > /dev/null
	ln -s ${REFLOC} scaffolds.fasta 2>&1 > /dev/null
	bash ${SCRIPTDIR}/bwa_mem_v1.sh ${DATA} scaffolds auto ${SRCDIR} P ${MAXPROCESS} > /dev/null 2>&1
	rm ${DATA}.fastq scaffolds.fasta &> /dev/null

	#bwa-end to end
	echo 'BWA - backtrack'
	BWADIR=${OUTDIR}/bwa-backtrack
	BWADIR2=${BWADIR}
	if [ -d  ${BWADIR} ]; then
		rm -rf ${BWADIR}
	fi
	mkdir -p ${BWADIR}/index &> /dev/null
	cd ${BWADIR};
	ln -s ${ECVLOC} . 2>&1 > /dev/null
	ln -s ${REFLOC} scaffolds.fasta 2>&1 > /dev/null
	bash ${SCRIPTDIR}/bwa_backtrack_v4.sh ${DATA} scaffolds auto ${SRCDIR} P ${MAXPROCESS} 2>&1 > /dev/null
	rm ${DATA}.fastq scaffolds.fasta &> /dev/null
fi