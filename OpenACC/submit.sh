#! /bin/bash -l
# The -l specifies that we are loading modules
#
## Walltime limit
#$ -l h_rt=12:30:00
#
## Give the job a name.
#$ -N OpenACC
#
## Redirect error output to standard output
#$ -j y
#
## What project to use. "paralg" is the project for the class
#$ -P paralg
#
## Ask for nodes with 4 cores, 4 cores total (so 1 node)
#$ -pe omp 4
#$ -l gpus=1.0
## -l gpu_c=6.0

# Want more flags? Look here:
# http://www.bu.edu/tech/support/research/system-usage/running-jobs/submitting-jobs/

# Load the correct modules
module load pgi
module load gcc


pgc++ -std=c++11 -Minfo=accel -acc -ta=tesla OpenACC.cpp -o OPENACC
./OPENACC 32768


exit

