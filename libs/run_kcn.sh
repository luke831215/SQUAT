#!/bin/bash
# kmer counting & indexing
#[out_fpath, dataset, fastq_path]: from raw reads to all-reads_joinP_rn_sample.csv

if [ "$#" -lt "2" ]; then
	echo "Using:" $0 "[output_dir] [dataname] [fastq_path]"
	exit
fi

EXECDIR="$( cd "$(dirname "$0")" ; pwd)"
OUTDIR=$1
DATA=$2
ECVLOC=$3
UTILDIR=${EXECDIR}

mkdir -p ${OUTDIR}
#cd ${UTILDIR}/kcnLandscape
${UTILDIR}/KMC/bin/kmc -k25 -r -cs4294967295 -ci1 $ECVLOC ${OUTDIR}/${DATA}-k25.res ${OUTDIR}

# generate histo.txt
${UTILDIR}/KMC/bin/kmc_tools transform ${OUTDIR}/${DATA}-k25.res histogram ${OUTDIR}/${DATA}-k25_histo.txt

# generate all-reads_class.txt
#mkdir -p ${OUTDIR}/${2}/classify

#time ./merge-ids3.pl ${OUTDIR}/${2} ${1}

# generate all-class2.csv
#time ./MergeLabels3 ${OUTDIR}/${2}/classify/all-reads_class.txt ${OUTDIR}/${2}/ids/${1}_ecv_0_repeats.stock ${OUTDIR}/${2}/classify/all-reads_class2.csv

# use histo.htm to select t1 t2 (manual)
${UTILDIR}/kcnLandscape/histoAnalyze ${OUTDIR}/${DATA}-k25_histo.txt ${OUTDIR}/${DATA}-k25_histo.htm

# generate all-reads_joinP_sample.csv
#time ./kcnLandscapeNewPtile ${OUTDIR}/${2}/${1}_ecv.fastq ${UTILDIR}/repeat-detect/${1}/${1}-k25.res ${OUTDIR}/${2}/classify/all-reads_class2.csv ${OUTDIR}/${2}/all-reads_joinP_sample.csv 0
