## Check all the bioinformatics tools can be called from current working directory.
for tool in cutadapt fastx_trimmer seqtk prinseq-lite.pl bwa samtools bedtools bedGraphToBigWig sort-bed
do command -v ${tool} >/dev/null 2>&1 || { echo >&2 ${tool}" is required. Please make sure you can call the bioinformatics tools from your current working directoryb.  Aborting."; exit 1; }
done

#exec 1>test.log 2>&1
exec > >(tee ${OUTPUT}/proseq2.0_Run_${tmp}.log)
exec 2>&1

## Write out the bigWigs.
echo " "
echo "Writing bigWigs:"
for f in ${TMPDIR}/*.sort.bam
 do
   j=`echo $f | awk -F"/" '{print $NF}' | rev | cut -d \. -f 3- |rev `
   echo $j > ${OUTPUT}/${j}.align.log
   if [ "${RNA5}" == "R1_5prime" ] ; then
     if [ "${OPP}" == "FALSE" ] ; then
       if [ "${MAP5}" == "TRUE" ] ; then  ## report The 5' end of the RNA. Danko lab leChRO-Seq protocol is on the 5' of _R1 readl, same strand of R1 ($9) 
          bedtools bamtobed -bedpe -mate1 -i $f 2> ${TMPDIR}/kill.warnings | awk 'BEGIN{OFS="\t"} ($9 == "+") {print $1,$2,$2+1,$7,$8,$9}; ($9 == "-") {print $1,$3-1,$3,$7,$8,$9}' | gzip > ${TMPDIR}/$j.bed.gz
       else ## report The 3' end of the RNA.  Danko lab leChRO-Seq protocol is on the 5 prime of _R2 read, opposite strand of R2 (R2 strand $10, R1 strand $9)
          bedtools bamtobed -bedpe -mate1 -i $f 2> ${TMPDIR}/kill.warnings | awk 'BEGIN{OFS="\t"} ($10 == "-") {print $1,$6-1,$6,$7,$8,$9}; ($10 == "+") {print $1,$5,$5+1,$7,$8,$9}' | gzip > ${TMPDIR}/$j.bed.gz
       fi
     elif [ "${OPP}" == "TRUE" ] ; then
       if [ "${MAP5}" == "TRUE" ] ; then  ## report The 5' end of the RNA. 
          bedtools bamtobed -bedpe -mate1 -i $f 2> ${TMPDIR}/kill.warnings | awk 'BEGIN{OFS="\t"} ($9 == "+") {print $1,$2,$2+1,$7,$8,$10}; ($9 == "-") {print $1,$3-1,$3,$7,$8,$10}' | gzip > ${TMPDIR}/$j.bed.gz
       else ## report The 3' end of the RNA.
          bedtools bamtobed -bedpe -mate1 -i $f 2> ${TMPDIR}/kill.warnings | awk 'BEGIN{OFS="\t"} ($10 == "-") {print $1,$6-1,$6,$7,$8,$10}; ($10 == "+") {print $1,$5,$5+1,$7,$8,$10}' | gzip > ${TMPDIR}/$j.bed.gz
       fi
     fi
   elif [ "${RNA5}" == "R2_5prime" ] ; then
     if [ "${OPP}" == "FALSE" ] ; then
       if [ "${MAP5}" == "TRUE" ] ; then #report the 5 prime end of RNA, in Engreitz data is 5 prime end of R2, same strand
          bedtools bamtobed -bedpe -mate1 -i $f 2> ${TMPDIR}/${j}.kill.warnings | awk 'BEGIN{OFS="\t"} ($10 == "+") {print $1,$5,$5+1,$7,$8,$10}; ($10 == "-") {print $1,$6-1,$6,$7,$8,$10}'|gzip > ${TMPDIR}/${j}.bed.gz      
       else  ## report the 3-prime end of the RNA, in Engreitz data is the 5' end of R1 read, but opposite strand
          bedtools bamtobed -bedpe -mate1 -i $f 2> ${TMPDIR}/${j}.kill.warnings | awk 'BEGIN{OFS="\t"} ($9 == "+") {print $1,$2,$2+1,$7,$8,$10}; ($9 == "-") {print $1,$3-1,$3,$7,$8,$10}' |gzip  > ${TMPDIR}/${j}.bed.gz
       fi
     elif [ "${OPP}" == "TRUE" ] ; then
       if [ "${MAP5}" == "TRUE" ] ; then #report the 5 prime end of RNA, in Engreitz data is 5 prime end of R2, same strand
          bedtools bamtobed -bedpe -mate1 -i $f 2> ${TMPDIR}/${j}.kill.warnings | awk 'BEGIN{OFS="\t"} ($10 == "+") {print $1,$5,$5+1,$7,$8,$9}; ($10 == "-") {print $1,$6-1,$6,$7,$8,$9}'|gzip > ${TMPDIR}/${j}.bed.gz      
       else  ## report the 3-prime end of the RNA, in Engreitz data is the 5' end of R1 read, but opposite strand
          bedtools bamtobed -bedpe -mate1 -i $f 2> ${TMPDIR}/${j}.kill.warnings | awk 'BEGIN{OFS="\t"} ($9 == "+") {print $1,$2,$2+1,$7,$8,$9}; ($9 == "-") {print $1,$3-1,$3,$7,$8,$9}' |gzip  > ${TMPDIR}/${j}.bed.gz
       fi
     fi
   fi
   
   echo 'Number of mappable reads:' >> ${OUTPUT}/${j}.align.log
   readCount=`zcat ${TMPDIR}/$j.bed.gz | grep "" -c`
   echo ${readCount} >> ${OUTPUT}/${j}.align.log
   
   ## Remove rRNA and reverse the strand (PRO-seq).
   zcat ${TMPDIR}/$j.bed.gz | grep "rRNA\|chrM" -v | grep "_" -v | sort-bed - | gzip > ${TMPDIR}/$j.nr.rs.bed.gz
   echo 'Number of mappable reads (excluding rRNA):' >> ${OUTPUT}/${j}.align.log
   echo `zcat ${TMPDIR}/$j.nr.rs.bed.gz | grep "" -c` >> ${OUTPUT}/${j}.align.log
   
   ## Convert to bedGraph ... Can't gzip these, unfortunately.
   bedtools genomecov -bg -i ${TMPDIR}/$j.nr.rs.bed.gz -g ${CHINFO} -strand + > ${TMPDIR}/$j\_plus.bedGraph
   bedtools genomecov -bg -i ${TMPDIR}/$j.nr.rs.bed.gz -g ${CHINFO} -strand - > ${TMPDIR}/$j\_minus.noinv.bedGraph
   
   ## Invert minus strand.
   cat ${TMPDIR}/$j\_minus.noinv.bedGraph | awk 'BEGIN{OFS="\t"} {print $1,$2,$3,-1*$4}' > ${TMPDIR}/$j\_minus.bedGraph ## Invert read counts on the minus strand.
   
   ## normalized by RPM
   cat ${TMPDIR}/$j\_plus.bedGraph | awk 'BEGIN{OFS="\t"} {print $1,$2,$3,$4*1000*1000/'$readCount'/1}' > ${TMPDIR}/$j\_plus.rpm.bedGraph  &
   cat ${TMPDIR}/$j\_minus.bedGraph | awk 'BEGIN{OFS="\t"} {print $1,$2,$3,$4*1000*1000/'$readCount'/1}' > ${TMPDIR}/$j\_minus.rpm.bedGraph  &
   wait
   ## Then to bigWig (nomalized and non-nomrmalized ones)
   bedGraphToBigWig ${TMPDIR}/$j\_plus.rpm.bedGraph ${CHINFO} ${OUTPUT}/$j\_plus.rpm.bw 
   bedGraphToBigWig ${TMPDIR}/$j\_minus.rpm.bedGraph ${CHINFO} ${OUTPUT}/$j\_minus.rpm.bw 
   bedGraphToBigWig ${TMPDIR}/$j\_plus.bedGraph ${CHINFO} ${OUTPUT}/$j\_plus.bw 
   bedGraphToBigWig ${TMPDIR}/$j\_minus.bedGraph ${CHINFO} ${OUTPUT}/$j\_minus.bw 
wait

   rm ${TMPDIR}/$j.nr.rs.bed.gz ${TMPDIR}/$j.bed.gz ${TMPDIR}/$j*.bedGraph
 done

fi #*/



#############################################
## CLEANUP!
#rm -Rf ${TMPDIR}

