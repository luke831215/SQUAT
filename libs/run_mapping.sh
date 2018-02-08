#!/bin/bash

while getopts ":h" opt; do
  case $opt in
    h)
      echo -e "arg1: output directory"
      echo -e "arg2: name of the data"
      echo -e "arg3: number of reads"
      echo -e "arg4: path of reference genome"
      echo -e "arg5: path of reads fastq file" 
      echo -e "arg6: number of threads" 
      echo -e "arg7(optional): bwaOnly or bowtie2Only\n" >&2
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

mkdir -p ${OUTDIR}/imgs

if [[ $8 != "bwaOnly" ]]; then

	#bowtie2 - local alignment
	echo 'bowtie2 - local alignment'
	BOWTIEDIR=${OUTDIR}/bowtie2-local
	BOWTIEDIR1=${BOWTIEDIR}
	if [ -d  ${BOWTIEDIR} ]; then
		rm -rf ${BOWTIEDIR}
	fi
	mkdir -p ${BOWTIEDIR}
	cd ${BOWTIEDIR}; mkdir -p Pauto
	ln -s ${ECVLOC} .
	ln -s ${REFLOC} scaffolds.fasta
	${SRCDIR}/bowtie2/bowtie2-build scaffolds.fasta ${DATA} > build.log
	${SRCDIR}/bowtie2/bowtie2 -p ${MAXPROCESS} --local -x ${DATA} -U ${DATA}_ecv.fastq -S Pauto/${DATA}_ecv_all.sam -k 20 --omit-sec-seq --reorder > ${DATA}_ecv_all.log
	bash ${SCRIPTDIR}/bowtie2_benchmark_v2.sh ${DATA}_ecv scaffolds auto I

	#bowtie2 - end to end
	echo 'bowtie2 - end to end'
	BOWTIEDIR=${OUTDIR}/bowtie2-endtoend
	BOWTIEDIR2=${BOWTIEDIR}
	if [ -d  ${BOWTIEDIR} ]; then
		rm -rf ${BOWTIEDIR}
	fi
	mkdir -p ${BOWTIEDIR}
	cd ${BOWTIEDIR}; mkdir -p Pauto
	ln -s ${ECVLOC} .
	ln -s ${REFLOC} scaffolds.fasta
	${SRCDIR}/bowtie2/bowtie2-build scaffolds.fasta ${DATA} > build.log
	${SRCDIR}/bowtie2/bowtie2 -p ${MAXPROCESS} -x ${DATA} -U ${DATA}_ecv.fastq -S Pauto/${DATA}_ecv_all.sam -k 20 --omit-sec-seq --reorder > ${DATA}_ecv_all.log
	bash ${SCRIPTDIR}/bowtie2_benchmark_v2.sh ${DATA}_ecv scaffolds auto I

fi

if [[ $8 != "bowtie2Only" ]]; then
	#bwa-mem
	echo 'bwa - mem'
	BWADIR=${OUTDIR}/bwa-mem
	BWADIR1=${BWADIR}
	if [ -d  ${BWADIR} ]; then
		rm -rf ${BWADIR}
	fi
	mkdir -p ${BWADIR}
	cd ${BWADIR};
	ln -s ${ECVLOC} .
	ln -s ${REFLOC} scaffolds.fasta
	bash ${SCRIPTDIR}/bwa_mem_v1.sh ${DATA}_ecv scaffolds auto P
	
	#bwa-end to end
	echo 'bwa - end to end'
	BWADIR=${OUTDIR}/bwa-endtoend
	BWADIR2=${BWADIR}
	if [ -d  ${BWADIR} ]; then
		rm -rf ${BWADIR}
	fi
	mkdir -p ${BWADIR}
	cd ${BWADIR};
	ln -s ${ECVLOC} .
	ln -s ${REFLOC} scaffolds.fasta
	bash ${SCRIPTDIR}/bwa_endtoend_v2.sh ${DATA}_ecv scaffolds auto P

fi