# GenBenchmark version: 4.1
# Based on version 4 with backporting logics from version 7

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
	unset EDITDIST
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

	echo "[bwa] aln"
	${BWAPATH}/bwa aln ${EDITDIST} -t ${MAXPROC} "index/${REFGENOME}" "${DATANAME}.fastq" -f "${WORKDIR}/${DATANAME}_all.sai" 2>&1 | tee "${LOGSDIR}/${DATANAME}_aln"

	echo "[bwa] samse"
	${BWAPATH}/bwa samse "index/${REFGENOME}" "${WORKDIR}/${DATANAME}_all.sai" "${DATANAME}.fastq" -f "${WORKDIR}/${DATANAME}_all.sam" 2>&1 | tee "${LOGSDIR}/${DATANAME}_samse"
}


function do_extract {
	echo "[fastq ans ans]"

	LNDESC="\rextract (%s/2): %s processed"

	printf "%s" ""
	# filter out flags of 256 (not primary alignment) and 2048 (supplementary alignment)
	# it shall equal to the whole data
	${SAMPATH}/samtools view -S -h -F 2304 "${WORKDIR}/${DATANAME}_all.sam" 2> /dev/null | \
	grep -v "^@" | \
	awk -v FNAME="${DATADIR}/${DATANAME}.fastq" \
		-v RNFNAME="${DATADIR}/${DATANAME}_rn.fastq" \
		-v RNANAME="${DATADIR}/${GSTD%_gstd}_rn_gstd.ans" \
		-v LNDESC="${LNDESC}" \
	'BEGIN {
		printf(LNDESC, "1", NR)
	}{
		print "@"$1"\n"$10"\n+\n"$11 > FNAME

		if($10 !~ /N/) {
			print "@"$1"\n"$10"\n+\n"$11 > RNFNAME

			# filter out flag of 0x4 (read unmapped)
			if(and($2,0x4)!=0x4) {
				# we only want the "XT:A:U" tag
				for(i=12; i<=NF; i++) {
					if($i == "XT:A:U") {
						hasNoXAtag = 1

						# check if this record has "XA:Z:" tag
						for(j=NF; j>=12; j--) {
							# we do not allow record with this tag
							if($j ~ /(XA:Z:)/) {
								hasNoXAtag = 0

								break
							}
						}

						if(hasNoXAtag == 1) {
							# the gold standard shall only have CIGAR flag M
							if($6 !~ /I|D|N|S|H|P|=|X/) {
								# find out the answer and export it
								for(j=12; j<=NF; j++) {
									if($j ~ /(MD:Z:)/) {
										print $1"\t"substr($j, 6) > RNANAME

										break
									}
								}
							}
						}

						break
					}
				}
			}
		}

		# line count indicator
		if(NR%100000==0) {
			printf(LNDESC, "1", NR)
		}
	} END {
		printf(LNDESC, "1", NR)
		printf(" and %s completed", NR)
	}'
	printf '\n'

	printf "%s" "extract (2/2):"

	printf "%s" " convert"
	bash ${UTILDIR}/FastqToSfq.sh "${DATADIR}/${DATANAME}.fastq" > "${DATADIR}/${DATANAME}.sfq"
	printf '\n'
}

function do_ids {
	echo "[idlist]"

	printf "%s" "clean"
	mkdir -p "${IDLSDIR}" &> /dev/null
	rm -f "${IDLSDIR}/${RAW}"_* &> /dev/null
	printf "\n"

	LNDESC="\rextract (%s/2): %s processed"
	
	touch "${IDLSDIR}/${RAW}_2_unmappable.ids" 
	touch "${IDLSDIR}/${RAW}_3_mappable_multi.ids" 
	touch "${IDLSDIR}/${RAW}_4_mappable_unique_noerror.ids" 
	touch "${IDLSDIR}/${RAW}_5_mappable_unique_subonly.ids" 
	touch "${IDLSDIR}/${RAW}_6_mappable_unique_clips.ids" 
	touch "${IDLSDIR}/${RAW}_7_mappable_unique_others.ids" 
	#touch "${IDLSDIR}/${RAW}_8_mappable_unique_altsite.ids" 
	touch "${IDLSDIR}/${RAW}_8_mappable_special.ids" 
	touch "${IDLSDIR}/${RAW}_9_mappable_ref_contain_N.ids"

	printf "%s" ""
	cat "${WORKDIR}/${DATANAME}_all.sam" | grep -v "^@" | \
	awk -v IDSET1="${IDLSDIR}/${RAW}_1_contain_N.ids" \
		-v IDSET2="${IDLSDIR}/${RAW}_2_unmappable.ids" \
		-v IDSET3="${IDLSDIR}/${RAW}_3_mappable_multi.ids" \
		-v IDSET4="${IDLSDIR}/${RAW}_4_mappable_unique_noerror.ids" \
		-v IDSET5="${IDLSDIR}/${RAW}_5_mappable_unique_subonly.ids" \
		-v IDSET6="${IDLSDIR}/${RAW}_6_mappable_unique_clips.ids" \
		-v IDSET7="${IDLSDIR}/${RAW}_7_mappable_unique_others.ids" \
		-v IDSET8="${IDLSDIR}/${RAW}_8_mappable_special.ids" \
		-v IDSET9="${IDLSDIR}/${RAW}_9_mappable_ref_contain_N.ids" \
		-v IDSETA="${IDLSDIR}/${RAW}_0_reads.info.tmp" \
		-v REPEAT="${IDLSDIR}/${RAW}_0_repeats.stock" \
		-v LNDESC="${LNDESC}" \
	'BEGIN {
		printf(LNDESC, "1", NR)

		LastNID=-1
	}{
		# filter out reads having Ns
		if($10 ~ /N/) {
			# export the first record
			if($1!=LastNID) {
				print $1 > IDSET1
				print $1"\tN\t0" > IDSETA

				LastNID=$1
			}
		}

		# filter out records of the reads having Ns
		if($1!=LastNID) {
			# filter out flag of 0x4 (read unmapped)
			if(and($2,0x4)==0x4) {
				print $1 > IDSET2
				print $1"\tF\t0" > IDSETA
			}
			# filter out flags of 0x100 (not primary alignment) and 0x800 (supplementary alignment)
			else if(and($2,0x100)!=0x100 && and($2,0x800)!=0x800) {
				for(i=12; i<=NF; i++) {
					if($i == "XT:A:R") {
						for(j=12; j<=NF; j++) {
							if($j ~ /(X0:i:)/) {
								print $1 > IDSET3
								print $1"\t"substr($j,6) > REPEAT
								print $1"\tM\t"substr($j,6) > IDSETA

								break
							}
						}

						break
					} else if($i == "XT:A:N") {
						print $1 > IDSET9
						print $1"\tR\t1" > IDSETA

						break
					} else if($i == "XT:A:U") {
						hasNoXAtag = 1

						# check if this record has "XA:Z:" tag
						for(j=NF; j>=12; j--) {
							if($j ~ /(XA:Z:)/) {
								hasNoXAtag = 0

								break
							}
						}

						if(hasNoXAtag == 1) {
							# extract reads having only M CIGAR flag; they
							# are perfect match and sub-only mismatch reads
							if($6 !~ /I|D|N|S|H|P|=|X/) {
								# perfect matched reads having no misatches
								for(j=12; j<=NF; j++) {
									if($j ~ /(NM:i:)/) {
										if(int(substr($j,6))==0) {
											print $1 > IDSET4
											print $1"\tP\t1" > IDSETA
										} else {
											print $1 > IDSET5
											print $1"\tS\t1" > IDSETA
										}

										break
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
						} else {
							print $1 > IDSET8
							print $1"\tA\t1" > IDSETA
						}

						break
					}
				}
			}
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

	printf "%s" "extract (2/2):"

	printf "%s" " list"
	sort -k1 -n "${IDLSDIR}/${RAW}_0_reads.info.tmp" 2> /dev/null > "${IDLSDIR}/${RAW}_0_reads.info"
	rm "${IDLSDIR}/${RAW}_0_reads.info.tmp" 2> /dev/null

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
	do_build
	do_ids
elif [ "$5" == "S" ]; then
	do_stats
elif [ "$5" == "U" ]; then
	do_upload
elif [ "$5" == "W" ]; then
	do_wipe
fi

