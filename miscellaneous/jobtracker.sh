#$ -S /bin/bash
#$ -V
#$ -N DONE_pileup
#$ -e /dev/null
#$ -o /dev/null
#$ -l h_vmem=1M
#$ -q short.q@cbsubrito
#$ -t 1
#$ -m e
#$ -M acv46@cornell.edu
#$ -hold_jid 468482

# This is a dummy job that doesn't run until a different job is done
# Used to send an email at the end of an array job, to prevent an email for every task in a large job array
# can also be launched by running separately with '-hold_jid $JOB_ID' parameter
# include -N in qsub call to give details of job in email (e.g. ${NAME}, step, fail)
