Compile and execute the code in this directory connect to SCC node.  

First, SSH into the SCC and get into queue and reserve a GPU in an interactive session  
```qrsh -l h_rt=01:45:00 -pe omp 1 -P paralg -l gpus=1.0 -l gpu_c=6.0```

Change to the project repository directory and get in /OpenACC directory  
```cd ./OpenACC```

Compile PGI code, ``np`` is the number of particles in the grid
```asm
module load pgi
module load gcc
pgc++ -std=c++11 -Minfo=accel -acc -ta=tesla OpenACC.cpp -o OPENACC
./OPENACC np
```

Request a gpu node to run if the number of particle is larger than 30,000, change the number after```./OPENACC``` in submit.sh to the number of particles and do

```qsub submit.sh```

