#1/bin/bash

if [ "$#" -lt "5" ]; then
	echo "Using:" $0 "[DataName] [RefGenome] [PARASET] [SRCDIR] [(A)ll, (B)uild, (E)xtract, (I)ds, (P)ostAA, (S)tats, (U)pload, (W)ipe]"
	exit
fi

MAXPROC=$(($(grep -c ^processor /proc/cpuinfo)/2))

DATANAME=$1
REFGENOME=$2
PARASET=$3
SRCRDIR=$4


if [ "${3}" == "auto"  ]; then
	unset PARASET
elif [ "${3}" == "min" ]; then
	PARASET="-A1 -B1 -O1 -E1 -L0"
fi

RAW="${DATANAME%_rn*}"
GSTD="${RAW}_gstd"

#WORKDIR="P${3}"
WORKDIR="."
LOGSDIR="${WORKDIR}/log"
DATADIR="${WORKDIR}/final"
IDLSDIR="${WORKDIR}/ids"
IDXDIR="${WORKDIR}/index"

BWAPATH="${SRCRDIR}/bwa"
SAMPATH="${SRCRDIR}/samtools"
UTILDIR="${SRCRDIR}/libs/map_proc/utils"

function do_mkdir {
	mkdir -p "${LOGSDIR}" &> /dev/null
	mkdir -p "${DATADIR}" &> /dev/null
	mkdir -p "${IDLSDIR}" &> /dev/null
	mkdir -p "${IDXDIR}" &> /dev/null
}

function do_clean {
	echo "[clean]"
	do_mkdir
	rm -f "${LOGSDIR}/${RAW}"_* &> /dev/null
	rm -f "${DATADIR}/${RAW}"_* &> /dev/null
	rm -f "${IDLSDIR}/${RAW}"_* &> /dev/null
	rm -f "${IDXSDIR}/${RAW}"_* &> /dev/null
}

function do_wipe {
	echo "[wipe]"
	rm -rf ${WORKDIR}/ &> /dev/null
}

function do_build {
	echo "[bwa] bwt"
	${BWAPATH}/bwa index ${REFGENOME}.fasta -p index/${REFGENOME} 2>&1 | tee ${LOGSDIR}/${DATANAME}_bwt

	echo "[bwa] mem"
	${BWAPATH}/bwa mem ${PARASET} -t ${MAXPROC} -a  index/${REFGENOME} "${DATANAME}.fastq" \
	> "${WORKDIR}/${DATANAME}_all.sam" 2> >(tee "${LOGSDIR}/${DATANAME}_mem" >&2)
}

function do_extract {
	echo "[fastq and ans]"
	# filter out flags of 256 (not primary alignment) and 2048 (supplementary alignment)
	# it shall equal to the whole data
	${SAMPATH}/samtools view -S -h -F 2304 "${WORKDIR}/${DATANAME}_all.sam" 2> /dev/null | \
	grep -v "^@" | \
	awk -v FNAME="${DATADIR}/${DATANAME}.fastq" -v ANAME="${DATADIR}/${GSTD}.ans" '{
		print "@"$1"\n"$10"\n+\n"$11 > FNAME

		# filter out flag of 0x4 (read unmapped), and the mapping quality shall not be 0
		if(and($2,0x4)!=0x4 && int($5)>=1) {
			# the gold standard shall only have CIGAR flag M
			if($6 !~ /I|D|N|S|H|P|=|X/) {
				# find out the answer and export it
				for(i=2; i<=NF; i++) {
					if($i ~ /(MD:Z:)/) {
						print $1"\t"substr($i, 6) > ANAME
					}
				}
			}
		}
	}'

	echo "[convert to sfq]"
	bash ${UTILDIR}/FastqToSfq.sh "${DATADIR}/${DATANAME}.fastq" > "${DATADIR}/${DATANAME}.sfq"
}

function do_ids {
	echo "[idlist]"
	echo "${WORKDIR}"
	SYMLINK=`readlink -f "${DATANAME}.fastq"`
	RAWFILE="${SYMLINK%${SYMLINK##*/}}${RAW}.fastq"
	if [ ! -f "${RAWFILE}" ]; then
		echo "NO RAWFILE"
		exit 1
	fi

	printf "%s" "clean"
	mkdir -p "${IDLSDIR}" &> /dev/null
	rm -f "${IDLSDIR}/${RAW}"_* &> /dev/null
	printf "\n"

	LNDESC="\rextract (%s/3): %s processed"

	printf "%s" ""
	cat ${RAWFILE} | \
	awk 'NR%4==1 {printf "%s\t", substr($0, 2)} NR%4==2 {printf "%s\n", $0}' | \
	awk -v IDSET1="${IDLSDIR}/${RAW}_7_contain_N.ids" \
		-v IDSETA="${IDLSDIR}/${RAW}_0_reads.info.tmp1" \
		-v LNDESC="${LNDESC}" \
	'BEGIN {
		printf(LNDESC, "1", NR)
	}{
		if($2 ~ /N/) {
			print $1 > IDSET1
			print $1"\tN\t0" > IDSETA
		}

		# line count indicator
		if(NR%100000==0) {
			printf(LNDESC, "1", NR)
		}
	} END {
		printf(LNDESC, "1", NR)
		printf(" and %s completed", NR)
	}'
	printf "\n"

	touch "${IDLSDIR}/${RAW}_6_unmappable.ids" 
	touch "${IDLSDIR}/${RAW}_5_mappable_multi.ids" 
	touch "${IDLSDIR}/${RAW}_1_mappable_unique_noerror.ids" 
	touch "${IDLSDIR}/${RAW}_2_mappable_unique_subonly.ids" 
	touch "${IDLSDIR}/${RAW}_3_mappable_unique_clips.ids" 
	touch "${IDLSDIR}/${RAW}_4_mappable_unique_others.ids" 
	
	printf "%s" ""
	cat "${WORKDIR}/${DATANAME}_all.sam" | grep -v "^@" | \
	awk -v IDSET2="${IDLSDIR}/${RAW}_6_unmappable.ids" \
		-v IDSET3="${IDLSDIR}/${RAW}_5_mappable_multi.ids" \
		-v IDSET4="${IDLSDIR}/${RAW}_1_mappable_unique_noerror.ids" \
		-v IDSET5="${IDLSDIR}/${RAW}_2_mappable_unique_subonly.ids" \
		-v IDSET6="${IDLSDIR}/${RAW}_3_mappable_unique_clips.ids" \
		-v IDSET7="${IDLSDIR}/${RAW}_4_mappable_unique_others.ids" \
		-v IDSETA="${IDLSDIR}/${RAW}_0_reads.info.tmp2" \
		-v REPEAT="${IDLSDIR}/${RAW}_0_repeats.stock" \
		-v LNDESC="${LNDESC}" \
	'BEGIN {
		printf(LNDESC, "2", NR)

		LastNID=-1
	}{
		# the reads containing Ns have been classified in previous step
		# update the variable at the first time
		if($10 ~ /N/) {
			LastNID=$1
		}

		if($1!=LastNID) {
			# filter out flag of 0x4 (read unmapped)
			if(and($2,0x4)==0x4) {
				print $1 > IDSET2
				print $1"\tF\t0" > IDSETA
			}
			# filter out flags of 0x100 (not primary alignment) and 0x800 (supplementary alignment)
			else if(and($2,0x100)!=0x100 && and($2,0x800)!=0x800) {
				# uniquely mapped reads having mapping qulaity higher than 0
				# multiply mapped reads will be processed later
				if(int($5)>=1) {
					# extract reads having only M CIGAR flag; they
					# are perfect match and sub-only mismatch reads
					if($6 !~ /I|D|N|S|H|P|=|X/) {
						# perfect matched reads having no misatches
						for(i=12; i<=NF; i++) {
							if($i ~ /(NM:i:)/) {
								if(int(substr($i,6))==0) {
									print $1 > IDSET4
									print $1"\tP\t1" > IDSETA
								} else {
									print $1 > IDSET5
									print $1"\tS\t1" > IDSETA
								}
							}
						}
					}
					#reads containing clips
					else if($6 ~ /S|H/) {
						print $1 > IDSET6
						print $1"\tC\t1" > IDSETA
					}
					# reads of mixing types of mismatches
					else {
						print $1 > IDSET7
						print $1"\tO\t1" > IDSETA	
					}
				}
			}

			# counting for the number of read occurences
			if(and($2,0x4)!=0x4) {
				a[$1]+=$5
				b[$1]+=1
			}
		}

		# line count indicator
		if(NR%100000==0) {
			printf(LNDESC, "2", NR)
		}
	} END {
		printf(LNDESC, "2", NR)
		printf(" and %s completed", "0")

		# sort array numerically first
		n = asorti(a,sa,"@ind_num_asc")

		# export multiply mapped reads, where the sum
		# of mapping qualities of a read shall be 0
		for(i=1;i<=n;i++) {
			if(int(a[sa[i]])==0) {
				print sa[i] > IDSET3
				print sa[i]"\t"b[sa[i]] > REPEAT
				print sa[i]"\tM\t"b[sa[i]] > IDSETA
			}

			if(i%100000==0) {
				printf(LNDESC, "2", NR)
				printf(" and %s completed", i)
			}
		}

		printf(LNDESC, "2", NR)
		printf(" and %s completed", NR)
	}'
	printf "\n"

	printf "%s" "extract (3/3):"

	printf "%s" " list"
	sort -k1 -n "${IDLSDIR}/${RAW}_0_reads.info.tmp1" "${IDLSDIR}/${RAW}_0_reads.info.tmp2" 2> /dev/null > "${IDLSDIR}/${RAW}_0_reads.info"
	rm "${IDLSDIR}/${RAW}_0_reads.info.tmp1" "${IDLSDIR}/${RAW}_0_reads.info.tmp2" 2> /dev/null

	printf "%s" " count"
	wc -l `ls ${IDLSDIR}/${RAW}_*.ids` > "${IDLSDIR}/${RAW}_0_reads.cnt"

	printf "%s" " repeat"
	awk '{print $2}' "${IDLSDIR}/${RAW}_0_repeats.stock" | sort -n | uniq -c | \
	awk '{print $2"\t"$1}' > "${IDLSDIR}/${RAW}_0_repeats.dist"
	printf "\n"
}

function do_stats {
	echo "[stats]"

	printf "%-15s:\t" "${GSTD}_ans"
	java -cp ${UTILDIR} goldstdStats "${DATADIR}/${GSTD}.ans" | tee "${LOGSDIR}/${GSTD}_ans_stats"

	for d in {${GSTD},${DATANAME}}; do
		printf "%-15s:\t" "${d}"
		cat "${DATADIR}/${d}.sfq" 2> /duv/null | \
		awk \
		'{
			rl=length($2); total+=rl; max=max>rl?max:rl; min=min>rl||!min?rl:min; gc+=gsub("[GC]", "", $2)
		} END {
			printf("%d\t%d\t%d\t%d\t%d\t%d\t%.2f%%\n", NR, total, (total/NR), max, min, gc, (gc/total)*100)
		}' | tee "${LOGSDIR}/${d}_stats"
	done

	printf "\n"
}

function do_upload {
	echo "[upload data]"
	nd="${DATANAME/_ecv_/_${WORKDIR}_ecv_}"

	printf "%s" "${nd}"
	hadoop fs -rm "data/${nd}.sfq" &> /dev/null
	hadoop fs -put "${DATADIR}/${DATANAME}.sfq" "data/${nd}.sfq"
	printf "\n"

	printf "\n"

	echo "[upload ans]"
	ng="${GSTD/_ecv_/_${WORKDIR}_ecv_}"
	printf "%s" "${ng}"
	hadoop fs -rm "eval/ans/${ng}.ans" &> /dev/null
	hadoop fs -put "${DATADIR}/${GSTD}.ans" "eval/ans/${ng}.ans"
	printf "\n"
}

if [ "$5" == "A" ]; then
	do_clean
	do_build
	do_extract
	do_ids
elif [ "$5" == "B" ]; then
	do_clean
	do_build
elif [ "$5" == "E" ]; then
	do_mkdir
	do_extract
elif [ "$5" == "I" ]; then
	do_ids
elif [ "$5" == "P" ]; then
	do_clean
	rm -f "${DATADIR}/${RAW}"_* &> /dev/null
	do_build
	do_ids
elif [ "$5" == "S" ]; then
	do_stats
elif [ "$5" == "U" ]; then
	do_upload
elif [ "$5" == "W" ]; then
	do_wipe
fi
