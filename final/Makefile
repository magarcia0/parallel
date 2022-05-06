INCLUDE_DIRS = -I/opt/intel/compilers_and_libraries_2020.0.166/linux/mpi/intel64/include/
LIB_DIRS = -L/opt/intel/compilers_and_libraries_2020.0.166/linux/mpi/intel64/lib/debug -L/opt/intel/compilers_and_libraries_2020.0.166/linux/mpi/intel64/lib
MPICC = mpicc
CC=gcc

CDEFS=
CFLAGS= -O0 -fopenmp $(INCLUDE_DIRS) $(CDEFS)
MPIFLAGS= -g -Wall -O2 $(INCLUDE_DIRS) $(CDEFS)
XFLAGS= -O0 -g -pg $(INCLUDE_DIRS) $(CDEFS)
LIBS= 

HFILES= 
#CFILES= serial_matrixmult.c omp_matrixmult.c serial_gsit.c gsit.c
CFILES= serial_gsit.c gsit.c omp_gsit.c mpi_mult.c

SRCS= ${HFILES} ${CFILES}
OBJS= ${CFILES:.c=.o}

#all:	serial_matrixmult omp_matrixmult serial_gsit gsit
all:	serial_gsit gsit omp_gsit mpi_mult

clean:
	-rm -f *.o *.d
	#-rm -f serial_matrixmult omp_matrixmult serial_git gsit
	-rm -f serial_gsit gsit omp_gsit mpi_mult

distclean:
	-rm -f *.o *.d
	#-rm -f serial_matrixmult omp_matrixmult serial_gsit gsit
	-rm -f serial_gsit gsit omp_gsit mpi_mult

serial_gsit: serial_gsit.o
	$(CC) $(LDFLAGS) $(CFLAGS) -o $@ $@.o

gsit: gsit.o
	$(CC) $(LDFLAGS) $(CFLAGS) -o $@ $@.o

omp_gsit: omp_gsit.o
	$(CC) $(LDFLAGS) -fopenmp $(CFLAGS) -o $@ $@.o

mpi_gsit: mpi_gsit.o
	$(CC) $(LDFLAGS) -fopenmp $(CFLAGS) -o $@ $@.o

mpi_mult:	mpi_mult.o
	$(MPICC) $(MPIFLAGS) -o $@ mpi_mult.c $(LIB_DIRS)

#serial_matrixmult: serial_matrixmult.o
#	$(CC) $(LDFLAGS)  $(XFLAGS) -o $@ $@.o

#vmult: vmult.o
#	$(CC) $(LDFLAGS) -fopenmp $(CFLAGS) -o $@ $@.o

depend:

.c.o:
	$(CC) $(CFLAGS) -c $<
