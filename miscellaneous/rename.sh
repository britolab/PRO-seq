#!/bin/bash

WRK=/workdir/users/acv46/stool_PROSeq2/fastq/raw
newnames=/workdir/users/acv46/stool_PROSeq2/fastq/raw/rename.txt

cd $WRK
while read line; do

	cp $(echo ${line} | awk '{print $2}') ./$(echo ${line} | awk '{print $1}')

done < $newnames
