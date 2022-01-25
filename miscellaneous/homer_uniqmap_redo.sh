#$ -S /bin/bash
#$ -N uniqmap
#$ -V
#$ -o /workdir/users/acv46/mgmAssembly/log/uniqmap_$JOB_ID.out
#$ -e /workdir/users/acv46/mgmAssembly/log/uniqmap_$JOB_ID.err
#$ -wd /workdir/users/acv46/mgmAssembly
#$ -l h_vmem=200G
#$ -t 2
#$ -q long.q@cbsubrito

# creates uniqmap directory for homer
# allocate as much mem as possible

WRK=/workdir/users/acv46/mgmAssembly
UMAP=/workdir/users/acv46/stool_PROSeq2/homer
HOMER=/programs/homer/bin

list=$WRK/samples.txt
List=${list}
        DESIGN=$(sed -n "${SGE_TASK_ID}p" $List)
        samp=`basename "$DESIGN"`

INFA=$WRK/${samp}_CAB/metaSPAdes/contigs.fasta
OUT=$UMAP/${samp}
mkdir -p $OUT

cd $OUT

#echo "${samp} -- splitting contigs.fasta"

##cat $INFA | sed "s/_cov.*//g" | awk '{
##        if (substr($0, 1, 1)==">") {filename=(substr($0,2) ".fa")}
##        print $0 > filename
##}'

echo "${samp} -- starting getMappableRegions"
$HOMER/getMappableRegions 1000000000 50 *.fa > mappableregions.50nt.txt
##echo "${samp} -- starting uniqmap"
##$HOMER/homerTools special uniqmap uniqmapDirectory/ mappableregions.50nt.txt

echo "${samp} -- removing temp fastas"
rm $OUT/*0.fa &
rm $OUT/*1.fa &
rm $OUT/*2.fa &
rm $OUT/*3.fa &
rm $OUT/*4.fa &
rm $OUT/*5.fa &
rm $OUT/*6.fa &
rm $OUT/*7.fa &
rm $OUT/*8.fa &
rm $OUT/*9.fa
