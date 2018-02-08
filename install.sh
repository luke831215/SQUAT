#!/bin/bash
set -e

INSTDIR="$( cd "$(dirname "$0")" ; pwd)"

cd ${INSTDIR}
if [ -d bowtie2 ]; then
	rm -r bowtie2
fi
if [ -d bwa ]; then
	rm -r bwa
fi

#install bowtie2
echo "install bowtie2-2.3.4"
wget https://downloads.sourceforge.net/project/bowtie-bio/bowtie2/2.3.4/bowtie2-2.3.4-linux-x86_64.zip 
unzip bowtie2-2.3.4-linux-x86_64.zip
rm bowtie2-2.3.4-linux-x86_64.zip
mv bowtie2-2.3.4-linux-x86_64 bowtie2

#install bwa
echo "install bwa-0.7.15"
wget https://downloads.sourceforge.net/project/bio-bwa/bwa-0.7.15.tar.bz2
tar jxvf bwa-0.7.15.tar.bz2
rm bwa-0.7.15.tar.bz2
mv bwa-0.7.15 bwa
cd bwa; make

#install KMC3
cd ${INSTDIR}/libs/gen_kmer
git clone https://github.com/refresh-bio/KMC.git
cd KMC; make DISABLE_ASMLIB=true

#compile kcnLandscape cpp
cd ${INSTDIR}/libs/kcnLandscape
g++ histAnalyze.cpp -o histAnalyze

#install python packages
echo "install python packages"
pip install fpdf