INCLUDE_DIRS = 
LIB_DIRS = 
CC=gcc

CDEFS=
CFLAGS= -O0 -fopenmp $(INCLUDE_DIRS) $(CDEFS)
XFLAGS= -O0 -g -pg $(INCLUDE_DIRS) $(CDEFS)
LIBS= 

HFILES= 
#CFILES= serial_matrixmult.c omp_matrixmult.c serial_gsit.c gsit.c
CFILES= serial_gsit.c gsit.c omp_gsit.c

SRCS= ${HFILES} ${CFILES}
OBJS= ${CFILES:.c=.o}

#all:	serial_matrixmult omp_matrixmult serial_gsit gsit
all:	serial_gsit gsit omp_gsit

clean:
	-rm -f *.o *.d
	#-rm -f serial_matrixmult omp_matrixmult serial_git gsit
	-rm -f serial_gsit gsit omp_gsit

distclean:
	-rm -f *.o *.d
	#-rm -f serial_matrixmult omp_matrixmult serial_gsit gsit
	-rm -f serial_gsit gsit omp_gsit

serial_gsit: serial_gsit.o
	$(CC) $(LDFLAGS) $(CFLAGS) -o $@ $@.o

gsit: gsit.o
	$(CC) $(LDFLAGS) $(CFLAGS) -o $@ $@.o

omp_gsit: omp_gsit.o
	$(CC) $(LDFLAGS) -fopenmp $(CFLAGS) -o $@ $@.o

#serial_matrixmult: serial_matrixmult.o
#	$(CC) $(LDFLAGS)  $(XFLAGS) -o $@ $@.o

#gewppomp: gewppomp.o
#	$(CC) $(LDFLAGS) -fopenmp $(CFLAGS) -o $@ $@.o

#gewpp: gewpp.o
#	$(CC) $(LDFLAGS) $(CFLAGS) -o $@ $@.o

#vmult: vmult.o
#	$(CC) $(LDFLAGS) -fopenmp $(CFLAGS) -o $@ $@.o


#gaussjord: gaussjord.o
#	$(CC) $(LDFLAGS) $(CFLAGS) -o $@ $@.o

#lufact: lufact.o
#	$(CC) $(LDFLAGS) $(CFLAGS) -o $@ $@.o

#testgen: testgen.o
#	$(CC) $(LDFLAGS) $(CFLAGS) -o $@ $@.o

depend:

.c.o:
	$(CC) $(CFLAGS) -c $<
