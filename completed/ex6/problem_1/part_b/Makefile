INCLUDE_DIRS = 
LIB_DIRS = 
CC = gcc
CPP = g++

CDEFS=
CFLAGS= -g -Wall -fopenmp $(INCLUDE_DIRS) $(CDEFS)
LIBS=

PRODUCT= all_newton

HFILES= 
CFILES= all_newton.c

SRCS= ${HFILES} ${CFILES}
OBJS= ${CFILES:.c=.o}

all:	${PRODUCT}

clean:
	-rm -f *.o *.NEW *~
	-rm -f ${PRODUCT} ${DERIVED} ${GARBAGE}

all_newton:	all_newton.c
	$(CC) $(CFLAGS) -o $@ all_newton.c -lm
