INCLUDE_DIRS = -I/opt/intel/compilers_and_libraries_2020.0.166/linux/mpi/intel64/include/
LIB_DIRS = -L/opt/intel/compilers_and_libraries_2020.0.166/linux/mpi/intel64/lib/debug -L/opt/intel/compilers_and_libraries_2020.0.166/linux/mpi/intel64/lib
MPICC = mpicc
MPICXX = mpicxx
CC = gcc
CPP = g++ -std=c++11

MPIFLAGS= -g -Wall -O0 $(INCLUDE_DIRS) $(CDEFS)
CFLAGS= -Wall -pedantic
COMP= -O0 -fopenmp $(CDEFS)

PRODUCT= piseriesreduce omp_piseriesreduce serial_piseriesreduce

CFILES= piseriesreduce.c omp_piseriesreduce.c serial_piseriesreduce.c

SRCS= ${CFILES}
OBJS= ${CFILES:.cpp=.o} ${CFILES:.c=.o}

all: ${PRODUCT}

piseriesreduce: piseriesreduce.c
	$(MPICC) $(MPIFLAGS) -o $@ piseriesreduce.c $(LIB_DIRS)

omp_piseriesreduce: omp_piseriesreduce.c
	$(CC) $(COMP) -o $@ omp_piseriesreduce.c

serial_piseriesreduce: serial_piseriesreduce.c
	$(CC) $(CFLAGS) -o $@ serial_piseriesreduce.c

clean:
	-rm -f *.o *.NEW *~
	-rm -f ${PRODUCT} ${DERIVED} ${GARBAGE}
