#$ -S /bin/bash
#$ -N test
#$ -V
#$ -e /workdir/users/acv46/stool_PROSeq2/log/test_$JOB_ID.err
#$ -o /workdir/users/acv46/stool_PROSeq2/log/test_$JOB_ID.out
#$ -wd /workdir/users/acv46/stool_PROSeq2
#$ -l h_vmem=1M
#$ -q long.q@cbsubrito
#$ -t 1-4
#$ -v TOTAL_TASKS=4

# this doesn't guarantee that tasks will finish in the order submitted
# better to ask each task the number of output files

WRK=/workdir/users/acv46/stool_PROSeq2/scripts

echo "running job $JOB_ID, task $SGE_TASK_ID"

if [[ "$SGE_TASK_ID" != "$TOTAL_TASKS" ]]; then

	echo "finished task $SGE_TASK_ID"

else

	echo "finished task $SGE_TASK_ID, the last task!"

fi
