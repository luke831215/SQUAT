#!/bin/bash

#Automatic exit from bash shell script on error
set -e

usage()
{
    echo -e "SQUAT: a Sequencing Quality Assessment Tool"
    echo -e "Usage: $0 seq1 seq2 ...  seqN [-o <output_dir>] [-r <ref_asm>]\n"
    echo "Optional args:"
    echo "--sample-size <int>    the read size for random sampling, default 1M"
    echo "--seed   <int>   Specify the seed for random sampling, default 0" 
}

function change_id {
    cat $2 | \
    awk -v OUT=$1 \
        -v MATCH=$(dirname "$1")/$3.ids \
    'BEGIN{
        id = 0
    }
    {
        if(NR%4==1){
            print "@"id > OUT
            print id"\t@"$0 > MATCH
            id+=1
        }
        else{
            print $0 > OUT
        }
    }
    END{
        print id
    }
    '
}


function to_abs {
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

#default parameter value
NUM_SAMPLE=1000000
SEED=0

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
    -o)
    OUTDIR="$( to_abs $2 )"
    shift # past argument
    shift # past value
    ;;
    --sample-size)
    NUM_SAMPLE=$2
    shift # past argument
    shift # past argument
    ;;
    --seed)
    SEED=$2
    shift # past argument
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

EXECDIR="$( cd "$(dirname "$0")" ; pwd)"

#write config file
mkdir -p ${OUTDIR} &> /dev/null

function do_squat {
    ORGECV="$( to_abs $1 )"
    #echo $ECVLOC
    xbase=${ORGECV##*/}
    DATA=${xbase%.*}
    SEQDIR=${OUTDIR}/${DATA}
    ECVLOC=${SEQDIR}/${DATA}.fastq
    #ECVLOC=$(dirname "$0")/${DATA}_ecv.fastq

    #delete if output dir already exists
    if [ -d ${SEQDIR} ]; then
        rm -r ${SEQDIR} &> /dev/null
    fi
    mkdir ${SEQDIR}
    touch ${SEQDIR}/${DATA}.log
    echo "Start examining ${DATA}" | tee -a ${SEQDIR}/${DATA}.log

    #echo "Calculate number of reads"  | tee -a ${SEQDIR}/${DATA}.log
    if [ "$FULLSET" == "YES" ]; then
        READSIZE=$( change_id ${ECVLOC} ${ORGECV} ${DATA} )
        NUM_SAMPLE=${READSIZE}
        echo "No. of reads: ${READSIZE}" | tee -a ${SEQDIR}/${DATA}.log
    else
        READSIZE=$(($(wc -l $ORGECV | cut -d ' ' -f 1)/4))
        if [ "$NUM_SAMPLE" -gt "$READSIZE" ]; then
            NUM_SAMPLE=$( change_id ${ECVLOC} ${ORGECV} ${DATA} )
            echo "No. of reads: ${NUM_SAMPLE}" | tee -a ${SEQDIR}/${DATA}.log
        else
            echo "sampling ${NUM_SAMPLE} out of ${READSIZE} records" | tee -a  ${SEQDIR}/${DATA}.log
            python3 ${EXECDIR}/library/rand_sample.py ${ECVLOC} ${ORGECV} ${READSIZE} ${NUM_SAMPLE} ${SEED}
        fi
    fi

}

for ((i=0;i<$NUM_SEQ;i++)); do
    do_squat ${SEQ_LIST[$i]}
done
