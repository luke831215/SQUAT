#!/bin/bash

if [ "$#" -lt "1" ]; then
	echo "Using:" $0 "[input file]"
	exit
elif [ ! -e $1 ]; then
	echo "File $1 does not exist."
	exit
fi

#while read line1; do
#	read line2
#	read line3
#	read line4
#
#	echo -e "${line1:1}\t$line2\t$line4"
#done < $1

cat $1 | awk 'NR%4==1 {printf "%s\t", substr($0, 2)} NR%4==2 {printf "%s\t", $0} NR%4==0 {printf "%s\n", $0}'
