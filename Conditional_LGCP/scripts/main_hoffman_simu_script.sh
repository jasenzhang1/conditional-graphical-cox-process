#### submit_job.sh START ####
#!/bin/bash
#$ -cwd
# error = Merged with joblog
#$ -o script_outputs/joblog.$JOB_ID
#$ -j y
## Edit the line below as needed:
#$ -l h_rt=9:00:00,h_data=4G
## Modify the parallel environment
## and the number of cores as needed:
#$ -pe shared 6
# Email address to notify
#$ -M $USER@mail #don't change this line, finds your email in the system 
# Notify when
#$ -m bea

# echo job info on joblog:
echo "Job $JOB_ID started on:   " `hostname -s`
echo "Job $JOB_ID started on:   " `date `
echo " "

# load the job environment:
. /u/local/Modules/default/init/modules.sh
## Edit the line below as needed:

module load apptainer
module load R


## substitute the command to run your code
## in the two lines below:
echo 'running conditional LGCP'
Rscript main_simulation.R > simu_logs/output.$JOB_ID 2>&1

# echo job info on joblog:
echo "Job $JOB_ID ended on:   " `hostname -s`
echo "Job $JOB_ID ended on:   " `date `
echo " "
#### submit_job.sh STOP ####