#!/bin/bash

INSTDIR="$( cd "$(dirname "$0")" ; pwd)"

cd ${INSTDIR}

#install bowtie2
printf "install bowtie2-2.3.4"
wget https://downloads.sourceforge.net/project/bowtie-bio/bowtie2/2.3.4/bowtie2-2.3.4-linux-x86_64.zip 
unzip bowtie2-2.3.4-linux-x86_64.zip
move bowtie2-2.3.4-linux-x86_64 bowtie2

#install bwa
printf "install bwa-0.7.15"
mkdir -p bwa; cd bwa
wget https://downloads.sourceforge.net/project/bio-bwa/bwa-0.7.15.tar.bz2
tar jxvf bwa-0.7.15.tar.bz2
mv bwa-0.7.15 bwa;
cd bwa; make