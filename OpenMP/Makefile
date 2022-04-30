
CC = g++
CFLAGS = -std=c++11 -fopenmp


vector: OMP_MD_vector_v0 OMP_MD_vector_v1 OMP_MD_vector_v2 OMP_MD_vector_v3 OMP_MD_vector_v4

basic: OMP_MD_basic_v0 OMP_MD_basic_v1 OMP_MD_basic_v2

OMP_MD_vector_v0: OMP_MD_vector_v0.cpp 
	$(CC) $(CFLAGS) OMP_MD_vector_v0.cpp  -o $@

OMP_MD_vector_v1: OMP_MD_vector_v1.cpp 
	$(CC) $(CFLAGS) OMP_MD_vector_v1.cpp  -o $@

OMP_MD_vector_v2: OMP_MD_vector_v2.cpp 
	$(CC) $(CFLAGS) OMP_MD_vector_v2.cpp  -o $@

OMP_MD_vector_v3: OMP_MD_vector_v3.cpp 
	$(CC) $(CFLAGS) OMP_MD_vector_v3.cpp  -o $@

OMP_MD_vector_v4: OMP_MD_vector_v4.cpp 
	$(CC) $(CFLAGS) OMP_MD_vector_v4.cpp  -o $@

OMP_MD_basic_v0: OMP_MD_basic_v0.cpp 
	$(CC) $(CFLAGS) OMP_MD_basic_v0.cpp  -o $@

OMP_MD_basic_v1: OMP_MD_basic_v1.cpp 
	$(CC) $(CFLAGS) OMP_MD_basic_v1.cpp  -o $@

OMP_MD_basic_v2: OMP_MD_basic_v2.cpp 
	$(CC) $(CFLAGS) OMP_MD_basic_v2.cpp  -o $@

.PHONY: clean 

clean:
	rm -f OMP_MD_vector_v0 OMP_MD_vector_v1 OMP_MD_vector_v2 OMP_MD_vector_v3 OMP_MD_vector_v4 OMP_MD_vector_v*.o \
		OMP_MD_basic_v0 OMP_MD_basic_v1 OMP_MD_basic_v2
