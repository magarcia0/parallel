INCLUDE_DIRS = -I/opt/intel/compilers_and_libraries_2020.0.166/linux/mpi/intel64/include/
LIB_DIRS = -L/opt/intel/compilers_and_libraries_2020.0.166/linux/mpi/intel64/lib/debug -L/opt/intel/compilers_and_libraries_2020.0.166/linux/mpi/intel64/lib
MPICC = mpicc
MPICXX = mpicxx
CC = icc
C++ = g++

CDEFS=
CFLAGS= -g -Wall $(INCLUDE_DIRS) $(CDEFS)
MPIFLAGS= -g -Wall -O0 $(INCLUDE_DIRS) $(CDEFS)
LIBS=

PRODUCT= mpi_double_train mpi_float_train

HFILES= 
CFILES= mpi_double_train.cpp mpi_float_train.cpp

SRCS= ${HFILES} ${CFILES}
OBJS= ${CFILES:.c=.o} ${CFILES:.cpp=.o}

all: ${PRODUCT}	


#mpi_trap3:	mpi_trap3.c
#	$(MPICC) $(MPIFLAGS) -o $@ mpi_trap3.c $(LIB_DIRS) -lm

#mpi_lrsum:	mpi_lrsum.cpp
#	$(MPICXX) $(MPIFLAGS) -o $@ mpi_lrsum.cpp $(LIB_DIRS)

mpi_double_train:	mpi_double_train.cpp
	$(MPICXX) $(MPIFLAGS) -o $@ mpi_double_train.cpp $(LIB_DIRS)

mpi_float_train:	mpi_float_train.cpp
	$(MPICXX) $(MPIFLAGS) -o $@ mpi_float_train.cpp $(LIB_DIRS)

#accel:	accel.cpp
#	$(C++) $(CFLAGS) -o $@ accel.cpp -lm

#lrsum:	lrsum.cpp
#	$(C++) $(CFLAGS) -o $@ lrsum.cpp

#rrsum:	rrsum.cpp
#	$(C++) $(CFLAGS) -o $@ rrsum.cpp

#trap3:	trap3.cpp
	#$(C++) $(CFLAGS) -o $@ trap3.cpp

#sum_train:	sum_train.cpp
#	$(C++) $(CFLAGS) -o $@ sum_train.cpp

#train:	train.cpp
#	$(C++) $(CFLAGS) -o $@ train.cpp

#timeinterp: timeinterp.o
#	$(CC) $(LDFLAGS) $(CFLAGS) -o $@ $@.o $(LIBS) -lm

#csvtostatic: csvtostatic.o
#	$(CC) $(LDFLAGS) $(CFLAGS) -o $@ $@.o $(LIBS) -lm

#timeprofiles: timeprofiles.o
#	$(CC) $(LDFLAGS) $(CFLAGS) -o $@ $@.o $(LIBS)

#omp_timeprofiles: omp_timeprofiles.c
#	$(CC) $(LDFLAGS) $(CFLAGS) -o $@ $@.c $(LIBS)
clean:
	-rm -f *.o *.NEW *~
	-rm -f ${PRODUCT} ${DERIVED} ${GARBAGE}
