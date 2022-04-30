#! /bin/bash

rm -f local_basic_*.txt local_vector_*.txt

for ((x = 64; x<300; x = x*2));do
    for ((y = 1;y <20; y = y * 2));do
        echo "T=$y np= $x"
        OMP_NUM_THREADS=$y ../OMP_MD_vector_v2 64 $x >> local_vector_2.txt
        OMP_NUM_THREADS=$y ../OMP_MD_vector_v3 64 $x >> local_vector_3.txt
        OMP_NUM_THREADS=$y ../OMP_MD_vector_v3 64 $x >> local_vector_4.txt

    done
done