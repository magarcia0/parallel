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
CFILES= serial_gsit.c gsit.c omp_gsit.c mpi_mult.c mpi_gsit.c

SRCS= ${HFILES} ${CFILES}
OBJS= ${CFILES:.c=.o}

all:	serial_gsit gsit omp_gsit mpi_mult mpi_gsit

clean:
	-rm -f *.o *.d
	-rm -f serial_gsit gsit omp_gsit mpi_mult mpi_gsit

distclean:
	-rm -f *.o *.d
	-rm -f serial_gsit gsit omp_gsit mpi_mult mpi_gsit

serial_gsit: serial_gsit.o
	$(CC) $(LDFLAGS) $(CFLAGS) -o $@ $@.o

gsit: gsit.o
	$(CC) $(LDFLAGS) $(CFLAGS) -o $@ $@.o

omp_gsit: omp_gsit.o
	$(CC) $(LDFLAGS) -fopenmp $(CFLAGS) -o $@ $@.o

mpi_gsit: mpi_gsit.o
	$(MPICC) $(MPIFLAGS) -o $@ mpi_gsit.c $(LIB_DIRS)

mpi_mult:	mpi_mult.o
	$(MPICC) $(MPIFLAGS) -o $@ mpi_mult.c $(LIB_DIRS)

depend:

.c.o:
	$(CC) $(CFLAGS) -c $<
