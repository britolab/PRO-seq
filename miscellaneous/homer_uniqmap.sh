#$ -S /bin/bash
#$ -N uniqmap
#$ -V
#$ -o /workdir/users/acv46/stool_PROSeq2/log/uniqmap_$JOB_ID.out
#$ -e /workdir/users/acv46/stool_PROSeq2/log/uniqmap_$JOB_ID.err
#$ -wd /workdir/users/acv46/stool_PROSeq2
#$ -l h_vmem=120G
#$ -t 1-2
#$ -q long.q@cbsubrito

# creates uniqmap directory for homer
# for the metaspades contigs tested, uses >100 Gb

WRK=/workdir/users/acv46/stool_PROSeq2
UMAP=$WRK/homer
HOMER=/programs/homer/bin

list=$UMAP/samples.txt
List=${list}
        DESIGN=$(sed -n "${SGE_TASK_ID}p" $List)
        samp=`basename "$DESIGN"`

INFA=$WRK/out3/${samp}_ref/${samp}_bins_plasmids.fa
OUT=$UMAP/${samp}
mkdir -p $OUT

cd $OUT

#echo "${samp} -- splitting contigs.fasta"

##cat $INFA | sed "s/_cov.*//g" | awk '{
##        if (substr($0, 1, 1)==">") {filename=(substr($0,2) ".fa")}
##        print $0 > filename
##}'

##echo "${samp} -- starting getMappableRegions"
##$HOMER/getMappableRegions 1000000000 50 *.fa > mappableregions.50nt.txt
echo "${samp} -- starting uniqmap"
$HOMER/homerTools special uniqmap uniqmapDirectory/ mappableregions.50nt.txt

#echo "${samp} -- removing temp fastas"
#rm $OUT/*0.fa
#rm $OUT/*1.fa
#rm $OUT/*2.fa
#rm $OUT/*3.fa
#rm $OUT/*4.fa
#rm $OUT/*5.fa
#rm $OUT/*6.fa
#rm $OUT/*7.fa
#rm $OUT/*8.fa
#rm $OUT/*9.fa
