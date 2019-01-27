#!/bin/bash

FILE_PREFIX=$1

vcftools --vcf data/${FILE_PREFIX}.vcf --plink
mv out.log reports/${FILE_PREFIX}.vcf_to_ped.log
mv out.map data/${FILE_PREFIX}.map
mv out.ped data/${FILE_PREFIX}.ped

exit