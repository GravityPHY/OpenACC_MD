#! /bin/bash
for ((x = 1024; x<6000; x = x*2));do
    for ((y = 1;y <20; y = y * 2));do
        echo "T=$y np= $x"
        OMP_NUM_THREADS=$y ./OMP_MD_vector_v2 1024 $x &> v2_t${y}_L1024_np${x}.txt
    done
done