INCLUDE_DIRS = 
LIB_DIRS = 
CC = gcc
CPP = g++

CDEFS=
CFLAGS= -g -Wall -fopenmp $(INCLUDE_DIRS) $(CDEFS)
LIBS=

#PRODUCT= bisection regulafalsi newton bruteroot
PRODUCT= regulafalsi allregula

HFILES= 
#CFILES= regulafalsi.c newton.c bruteroot
CFILES= regulafalsi.c allregula.c
CPPFILES= bisection.cpp

SRCS= ${HFILES} ${CFILES}
OBJS= ${CFILES:.c=.o}

all:	${PRODUCT}

clean:
	-rm -f *.o *.NEW *~
	-rm -f ${PRODUCT} ${DERIVED} ${GARBAGE}

regulafalsi: regulafalsi.c
	$(CC) $(CFLAGS) -o $@ regulafalsi.c -lm

allregula: allregula.c
	$(CC) $(CFLAGS) -o $@ allregula.c -lm
