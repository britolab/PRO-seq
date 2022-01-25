#! /bin/bash

# usage: sh Get_Lengths.sh dir/name.fasta outdir/lengths.txt
# gets lengths of sequences in fasta file

cat $1 | awk '$0 ~ ">" {print c; c=0;printf substr($0,2,100) "\t"; } $0 !~ ">" {c+=length($0);} END { print c; }' > $2
