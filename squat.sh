#!/bin/bash
#set -e

usage()
{
	echo -e "SQUAT: Sequencing Quality Assessment Tool"
	echo -e "Usage: $0 seq1 seq2 ...  seqN [-o <output_dir>] [-r <ref_seq>]\n"
	echo "Optional args:"
	echo "-t	--thread	<int>	Number of threads to use" 
	echo "-k	--keep	Don't flush The sam file after alignment" 
	echo "-s 	--subset 	<str>	Return the subset of sequencing reads according to the labels (in capitals, e.g. PSCO)" 
	echo "-g   <str>   Path to the reference genome file for GAGE benchmark tool" 
    echo "--gage    Activate gage mode, must specify reference genome (-R)"
    echo "--noextract    extract the zip file containing the report information"
}

function to_abs	{
	case $1 in
  		/*) absolute=$1;;
  		*) absolute=$PWD/$1;;
	esac
	echo $absolute
}

if [[ $# -eq 0 ]]; then
    usage
    exit 1
fi

SEQ_LIST=()
NUM_SEQ=0
while true; do
    if [[ $1 != "-"* && $# -gt 0 ]];
    then
        SEQ_LIST+=($( to_abs $1 ))
        NUM_SEQ=$(($NUM_SEQ+1))
        shift #pass seq files
    else
        break
    fi
done


while [[ $# -gt 0 ]]
do
key="$1"

case $key in
    -h|--help)
	usage
	exit 0
	;;
    -k|--keep)
	KEEP_SAM=YES
    shift
    ;;
    -o)
    OUTDIR="$( to_abs $2 )"
    shift # past argument
    shift # past value
    ;;
    -r)
    REFLOC="$( to_abs $2 )"
    shift # past argument
    shift # past value
    ;;
    -t|--thread)
    MAXPROC="$2"
    shift # past argument
    shift # past value
    ;;
    -s|--subset)
    SUBSET="$2"
    shift # past argument
    shift # past value
    ;;
    -g)
    GAGELOC="$( to_abs $2 )"
    shift # past argument
    shift # past value
    ;;
    --gage)
    GAGE=YES
    shift # past argument
    ;;
    --extract)
    EXTRACT=NO
    shift # past argument
    ;;
    *)    # unknown option
	echo "Unknown option: "$1 >&2
	exit 1
    ;;
esac
done

if [[ -z "$OUTDIR" || -z "$REFLOC" ]]; then
    usage
	exit 1
fi

if [[ -z "$MAXPROC" ]]; then
	MAXPROC=$(($(grep -c ^processor /proc/cpuinfo)/3))
fi

if [[ -z "$KEEP_SAM" ]]; then
	KEEP_SAM=NO
fi

if [[ -z "$EXTRACT" ]]; then
    EXTRACT=YES
fi


function do_squat {
    ECVLOC="$( to_abs $1 )"
    #echo $ECVLOC
    xbase=${ECVLOC##*/}
    DATA=${xbase%.*}
    SEQDIR=${OUTDIR}/${DATA}
    #ECVLOC=$(dirname "$0")/${DATA}_ecv.fastq
    
    echo "Calculate number of reads"
    READSIZE=$(($(wc -l $ECVLOC | cut -d ' ' -f 1) /4))

    EXECDIR="$( cd "$(dirname "$0")" ; pwd)"
    #shift $((OPTIND-1))
    #echo "$@"

    #delete if output dir already exists
    if [ -d ${SEQDIR} ]; then
        rm -rf ${SEQDIR} &> /dev/null
        mkdir ${SEQDIR} &> /dev/null
    fi

    #map reads to genome using alignment tools
    echo "map reads to genome using alignment tools"
    bash ${EXECDIR}/libs/run_mapping.sh $EXECDIR $SEQDIR $DATA $READSIZE $REFLOC $ECVLOC $MAXPROC 

    #quast evaluation
    echo "Evaluate genome assemblies"
    if [[ -z "$GAGELOC" && -z "$GAGE" ]]; then
        python ${EXECDIR}/quast/quast.py ${REFLOC} -o ${SEQDIR}/quast --min-contig 200 -t ${MAXPROC} 2>&1 > /dev/null
    else
        python ${EXECDIR}/quast/quast.py ${REFLOC} -o ${SEQDIR}/quast --min-contig 200 -t ${MAXPROC} -R ${GAGELOC} --gage 2>&1 > /dev/null 
    fi

    #analysis modules
    echo "Generate reports"
    mkdir -p ${SEQDIR}/subset &> /dev/null
    mkdir -p ${SEQDIR}/images &> /dev/null
    python ${EXECDIR}/analysis.py ${OUTDIR} ${ECVLOC} ${DATA} ${READSIZE} ${SUBSET}

    #pre-Q report
    #${EXECDIR}/libs/peQdist ${ECVLOC} ${SEQDIR}/${Data}_peQ --interleave

    #flush sam files
    if [ "$KEEP_SAM" == "NO" ]; then
        for tool in bowtie2-local bowtie2-endtoend bwa-mem bwa-endtoend; do
            rm ${SEQDIR}/${tool}/${DATA}_ecv_all.sam
        done
    fi

    #zip the file
    echo "Compress files"
    zip -r ${SEQDIR} ${SEQDIR}/* 2>&1 > /dev/null
    if [ "$EXTRACT" == "NO" ]; then
        rm -r ${SEQDIR} &> /dev/null
    fi
}

for ((i=0;i<$NUM_SEQ;i++)); do
    do_squat ${SEQ_LIST[$i]}
done