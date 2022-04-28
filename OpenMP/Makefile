
CC = g++
CFLAGS = -std=c++11 -fopenmp


v0: OMP_MD_vector_v0
v1: OMP_MD_vector_v1
v2: OMP_MD_vector_v2
v3: OMP_MD_vector_v3

OMP_MD_vector_v0: OMP_MD_vector_v0.cpp 
	$(CC) $(CFLAGS) OMP_MD_vector_v0.cpp  -o $@

OMP_MD_vector_v1: OMP_MD_vector_v1.cpp 
	$(CC) $(CFLAGS) OMP_MD_vector_v1.cpp  -o $@

OMP_MD_vector_v2: OMP_MD_vector_v2.cpp 
	$(CC) $(CFLAGS) OMP_MD_vector_v2.cpp  -o $@

OMP_MD_vector_v3: OMP_MD_vector_v3.cpp 
	$(CC) $(CFLAGS) OMP_MD_vector_v3.cpp  -o $@

.PHONY: clean 

clean:
	rm -f OMP_MD_vector_v0 OMP_MD_vector_v1 OMP_MD_vector_v2 OMP_MD_vector_v3 OMP_MD_vector_v*.o
