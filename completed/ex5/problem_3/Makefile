INCLUDE_DIRS = 
LIB_DIRS = 
CC=gcc

CDEFS=
CFLAGS= -O0 -fopenmp $(INCLUDE_DIRS) $(CDEFS)
XFLAGS= -O0 -g -pg $(INCLUDE_DIRS) $(CDEFS)
LIBS= 

HFILES= 
CFILES= serial_matrixmult.c omp_matrixmult.c
#CFILES= gewpp.c gewppomp.c serial_matrixmult.c vmult.c omp_matrixmult.c

SRCS= ${HFILES} ${CFILES}
OBJS= ${CFILES:.c=.o}

all:	serial_matrixmult omp_matrixmult
#all:	gewpp gewppomp serial_matrixmult vmult omp_matrixmult

clean:
	-rm -f *.o *.d
#	-rm -f gewpp gewppomp serial_matrixmult vmult omp_matrixmult
	-rm -f serial_matrixmult omp_matrixmult

distclean:
	-rm -f *.o *.d
	#-rm -f gewpp gewppomp serial_matrixmult vmult omp_matrixmult
	-rm -f serial_matrixmult omp_matrixmult

omp_matrixmult: omp_matrixmult.o
	$(CC) $(LDFLAGS) -fopenmp $(CFLAGS) -o $@ $@.o

serial_matrixmult: serial_matrixmult.o
	$(CC) $(LDFLAGS)  $(XFLAGS) -o $@ $@.o

#gewppomp: gewppomp.o
#	$(CC) $(LDFLAGS) -fopenmp $(CFLAGS) -o $@ $@.o

#gewpp: gewpp.o
#	$(CC) $(LDFLAGS) $(CFLAGS) -o $@ $@.o

#vmult: vmult.o
#	$(CC) $(LDFLAGS) -fopenmp $(CFLAGS) -o $@ $@.o

#gsit: gsit.o
#	$(CC) $(LDFLAGS) $(CFLAGS) -o $@ $@.o

#gaussjord: gaussjord.o
#	$(CC) $(LDFLAGS) $(CFLAGS) -o $@ $@.o

#lufact: lufact.o
#	$(CC) $(LDFLAGS) $(CFLAGS) -o $@ $@.o

#testgen: testgen.o
#	$(CC) $(LDFLAGS) $(CFLAGS) -o $@ $@.o

depend:

.c.o:
	$(CC) $(CFLAGS) -c $<
