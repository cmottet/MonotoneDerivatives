#!/bin/sh -l
# Submit job with the command: qsub Rjob
#
# To view the status of the job: qstat -u username
#
# Set the runtime limit (default 12 hours):
#$ -l h_rt=200:00:00
#
# Send email when the job is done (default: no email is sent)
#$ -m e
#
# Give the job a name (default: script name)
#$ -N runMinMaxTransformedHill_Horror_Pareto_GLP
#
# Tell the chore I need a minimum number of Memory
#$ -l mem_total=100G
#
## end of qsub options

# Load the version of R needed for computations
module load R/3.2.0

# Initialize the project structure
. /project/simulate/cmottet/Regularized2/programs/bash/init.sh

#temp=$(wc -l < $input_dir""ProbCovSyndataLogNorm.txt)
#file_lines=$(echo $temp | awk '{print $1}')

# Take off the headerline from the line count
N=256 #Number of rows in all the possible combination

for ((i = 1; i<=$N ; i++))
  do
echo $i                            
R -q --save --args $i < $prog_dir""R/runMinMaxSyntheticTransformed_Hill_Horror_Pareto_GLP.R
done