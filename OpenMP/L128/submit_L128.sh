#! /bin/bash -l
# The -l specifies that we are loading modules
#
## Walltime limit
#$ -l h_rt=12:00:00
#
## Give the job a name.
#$ -N OMP_MD_vector_L128
#
## Redirect error output to standard output
#$ -j y
#
## What project to use. "paralg" is the project for the class
#$ -P paralg
#
## Ask for 4 cores
#$ -pe omp 16
#
## Have the system send you mail when your job is aborted or ends
#$ -m bae

# Want more flags? Look here:
# http://www.bu.edu/tech/support/research/system-usage/running-jobs/submitting-jobs/

# Load the correct modules
# module load gcc/5.3.0  # compiler
# module load mpich/3.2  # consistent mpi compile

# Immediately form fused output/error file, besides the one with the default name.
exec >  ${SGE_O_WORKDIR}/${JOB_NAME}-${JOB_ID}.scc.out 2>&1

echo "Hello from process ${HOSTNAME}"

for ((x = 128; x<2000; x = x*2));do
    for ((y = 1;y <20; y = y*2));do
        echo "T=$y np= $x"
        OMP_NUM_THREADS=$y ../OMP_MD_vector_v2 128 $x &> v2_t${y}_L128_np${x}.txt
    done
done

exit

