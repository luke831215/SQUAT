#!/bin/bash
set -e

usage()
{
	echo -e "SQUAT: Sequencing Quality Assessment Tool"
	echo -e "Usage: $0 [-o <output_dir>] [-r <read_name>] [-i fastq_path] [-g <genome_path>]\n"
	echo "Optional args:"
	echo "-t	--thread	<int>	Number of thread to use" 
	echo "-k	--keep	Don't flush the sam file after alignment" 
	echo "-s 	<str>	return the subset of sequencing reads according to the labels (in capitals, e.g. PSCO)" 

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

while getopts ":o:r:g:i:t:s:hk" opt; do
  case $opt in
    h)
	  usage
	  exit 0
      ;;
    k)
	  KEEP_SAM=true
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
	s)
	  SUBSET=$OPTARG
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

if [[ -z "$KEEP_SAM" ]]; then
	KEEP_SAM=false
fi

echo "Calculate number of reads"
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

#analysis modules
echo "Generate reports"
#mkdir -p ${OUTDIR}/label_dis
#mkdir -p ${OUTDIR}/subset
python ${EXECDIR}/analysis.py ${OUTDIR} ${ECVLOC} ${DATA} ${READSIZE} ${SUBSET}

#flush sam files
if [ "$KEEP_SAM" = false ]; then
	for tool in bowtie2-local bowtie2-endtoend bwa-mem bwa-endtoend; do
		rm ${OUTDIR}/${tool}/Pauto/${DATA}_ecv_all.sam
	done
fi