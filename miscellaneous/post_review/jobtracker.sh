#$ -S /bin/bash
#$ -V
#$ -e /dev/null
#$ -o /dev/null
#$ -l h_vmem=1M
#$ -q long.q@cbsubrito
#$ -t 1
#$ -m e

# This is a dummy job that doesn't run until a different job is done
# Used to send an email at the end of an array job, to prevent an email for every task in a large job array
# can also be launched independently with '-hold_jid $JOB_ID' parameter
# include -N in qsub call to give details of job in email (e.g. ${NAME}, step, fail)
