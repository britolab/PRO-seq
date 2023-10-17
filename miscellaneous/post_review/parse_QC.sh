#!/bin/bash

cd /workdir/users/acv46/stool_PROSeq3/transcriptomes

for file in US*/temp/*/*.QC.log; do

	base=$(echo $file | awk -F$'/' '{print $4}')
	samp=$(echo $base | awk -F_ '{print $1}')
	repn=$(echo $base | awk -F_ '{print $2}')
	libt=$(echo $base | awk -F_ '{print $3}' | sed "s/.QC.log//g")

	inpu=$(sed '3q;d' $file)
	trm1=$(sed '8q;d' $file)
	trm2=$(sed '10q;d' $file)
	ddup=$(sed '12q;d' $file)

	if [ "$trm1" -gt "$trm2" ]; then
		trim=${trm2}
	else
		trim=${trm1}
	fi

	echo -e "${samp}\t${repn}\t${libt}\t${inpu}\t${trim}\t${ddup}"

done >> QC_log_14Dec2022.txt
