#!/bin/bash

cd data/

GENOME=papAnu2

# Download repeat-masked genome and TRF file from USCS genome browser

wget 'ftp://hgdownload.cse.ucsc.edu/goldenPath/'${GENOME}'/bigZips/'${GENOME}'.fa.out.gz' -O ${GENOME}.fa.out.gz
wget 'ftp://hgdownload.cse.ucsc.edu/goldenPath/'${GENOME}'/bigZips/'${GENOME}'.trf.bed.gz' -O ${GENOME}.trf.bed.gz

gunzip ${GENOME}.fa.out.gz
gunzip ${GENOME}.trf.bed.gz

cd ..

exit