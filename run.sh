#!/bin/bash

usage()
{
	echo -e "SQUAT: Sequencing Quality Assessment Tool"
	echo -e "Usage: $0 [-o <output_dir>] [-r <read_name>] [-l <read_length>] [-i fastq_path] [-g <genome_path>]\n"
	echo "Optional args:"
	echo "-t	--thread	<int>	Number of threads to use" 
}

if [[ $# -eq 0 || $(echo $1 | cut -c1) != "-" ]];then
    usage
    exit 1
fi

while getopts ":o:t:h" opt; do
  case $opt in
    h)
	  usage
	  exit 0
      ;;
    o)
	  OUTDIR=$OPTARG
	  ;;
	r)
	  DATA=$OPTARG
	  ;;
	l)
	  READSIZE=$OPTARG
	  ;;
	g)
	  REFLOC=$OPTARG
	  ;;
	i)
	  ECV_LOC=$OPTARG
	  ;;
	t)
	  MAXPROCESS=$OPTARG
	  ;;
	:)
	  echo "Missing option argument for -$OPTARG" >&2; exit 1
      ;;
    ?)
      echo "Unknown option: -"${OPTARG} >&2; exit 1
      ;;
    *) 
	  echo "Unimplemented option: -$OPTARG" >&2; exit 1
	  ;;
  esac
done

if [[ -z "$OUTDIR" || -z "$DATA" || -z "$READSIZE" || -z "$REFLOC" || -z "$ECV_LOC" ]]; then
	usage
	exit 1
fi

EXECDIR="$( cd "$(dirname "$0")" ; pwd)"
#shift $((OPTIND-1))
#echo "$@"


bash ${CRTDIR}/squat_libs/run_mapping.sh $OUTDIR $DATA $READSIZE $REFLOC $ECV_LOC
#bash ${CRTDIR}/squat_libs/run_kcn.sh