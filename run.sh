#!/bin/bash
set -e

usage()
{
	echo -e "SQUAT: Sequencing Quality Assessment Tool"
	echo -e "Usage: $0 [-o <output_dir>] [-r <read_name>] [-i fastq_path] [-g <genome_path>]\n"
	echo "Optional args:"
	echo "-t	--thread	<int>	Number of thread to use" 
}

function to_abs	{
	case $1 in
  		/*) absolute=$1;;
  		*) absolute=$PWD/$1;;
	esac

	echo $absolute
}

if [[ $# -eq 0 || $(echo $1 | cut -c1) != "-" ]];then
    usage
    exit 1
fi

while getopts ":o:r:g:i:t:h" opt; do
  case $opt in
    h)
	  usage
	  exit 0
      ;;
    o)
	  OUTDIR=$( to_abs ${OPTARG} )
	  ;;
	r)
	  DATA=$OPTARG
	  ;;
	g)
	  REFLOC=$( to_abs ${OPTARG} )
	  ;;
	i)
	  ECVLOC=$( to_abs ${OPTARG} )	
	  ;;
	t)
	  MAXPROC=$OPTARG
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

if [[ -z "$OUTDIR" || -z "$DATA" || -z "$REFLOC" || -z "$ECVLOC" ]]; then
	usage
	exit 1
fi

if [[ -z "$MAXPROC" ]]; then
	MAXPROC=$(($(grep -c ^processor /proc/cpuinfo)/3))
fi

READSIZE=$(($(wc -l $ECVLOC | cut -d ' ' -f 1) /4))
EXECDIR="$( cd "$(dirname "$0")" ; pwd)"
#shift $((OPTIND-1))
#echo "$@"

#delete if output dir already exists
if [ -d  ${OUTDIR} ]; then
	rm -rf ${OUTDIR}
fi

#map reads to genome using alignment tools
echo "map reads to genome using alignment tools"
bash ${EXECDIR}/libs/run_mapping.sh $EXECDIR $OUTDIR $DATA $READSIZE $REFLOC $ECVLOC $MAXPROC

#generate kmer info
echo "generate kmer stats"
bash ${EXECDIR}/libs/run_kcn.sh $OUTDIR/kcn_histo $DATA $ECVLOC

#concatenate tables and imgs to report.pdf
echo "Generate reports"
cp -r ${EXECDIR}/template/link ${OUTDIR}/
mkdir -p ${OUTDIR}/label_dis/imgs
python -i ${EXECDIR}/gen_report.py ${OUTDIR} ${ECVLOC} ${DATA} ${READSIZE}