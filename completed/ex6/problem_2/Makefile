INCLUDE_DIRS = 
LIB_DIRS = 
CC = gcc
CPP = g++

CDEFS=
CFLAGS= -g -Wall -fopenmp $(INCLUDE_DIRS) $(CDEFS)
LIBS=

PRODUCT= all_newton root_finder

HFILES= 
CFILES= all_newton.c root_finder.c

SRCS= ${HFILES} ${CFILES}
OBJS= ${CFILES:.c=.o}

all:	${PRODUCT}

clean:
	-rm -f *.o *.NEW *~
	-rm -f ${PRODUCT} ${DERIVED} ${GARBAGE}

all_newton:	all_newton.c
	$(CC) $(CFLAGS) -o $@ all_newton.c -lm

root_finder: root_finder.c
	$(CC) $(CFLAGS) -o $@ root_finder.c -lm
