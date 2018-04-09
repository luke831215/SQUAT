#!/bin/bash
set -e

INSTDIR="$( cd "$(dirname "$0")" ; pwd)"

cd ${INSTDIR}
chmod 764 squat.sh

if [ -d bwa ]; then
	rm -r bwa
fi

if [ -d quast ]; then
	rm -r quast
fi

#install bwa
echo "install bwa-0.7.15"
wget https://downloads.sourceforge.net/project/bio-bwa/bwa-0.7.15.tar.bz2
tar jxvf bwa-0.7.15.tar.bz2
rm bwa-0.7.15.tar.bz2*
mv bwa-0.7.15 bwa
cd bwa; make

#install quast
echo "install quast"
cd ${INSTDIR}
wget https://downloads.sourceforge.net/project/quast/quast-4.6.3.tar.gz
tar -xzf quast-4.6.3.tar.gz
mv quast-4.6.3 quast;
rm quast-4.6.3.tar.gz*
cp ${INSTDIR}/library/quast_code/* ${INSTDIR}/quast/quast_libs

#build readQdist
echo "build files for generating pre-assembly quality report"
cd ${INSTDIR}/library/preQ
bash ${INSTDIR}/library/preQ/_build.sh

#install python packages
echo "install python packages"
pip install beautifulsoup4
pip install lxml