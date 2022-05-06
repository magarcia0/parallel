INCLUDE_DIRS = 
LIB_DIRS = 
#CC = icc
CC = gcc
CPP = g++

CDEFS=
CFLAGS= -g -Wall -fopenmp $(INCLUDE_DIRS) $(CDEFS)
LIBS=

#PRODUCT= bisection regulafalsi newton bruteroot
PRODUCT= new_newton

HFILES= 
#CFILES= regulafalsi.c newton.c bruteroot.c
CFILES= new_newton.c 
#CPPFILES= bisection.cpp

SRCS= ${HFILES} ${CFILES}
OBJS= ${CFILES:.c=.o}

all:	${PRODUCT}

clean:
	-rm -f *.o *.NEW *~
	-rm -f ${PRODUCT} ${DERIVED} ${GARBAGE}

#bruteroot:	bruteroot.c
#	$(CC) $(CFLAGS) -o $@ bruteroot.c -lm

#bisection:	bisection.cpp
#	$(CPP) $(CFLAGS) -o $@ bisection.cpp -lm

#regulafalsi:	regulafalsi.c
#	$(CC) $(CFLAGS) -o $@ regulafalsi.c -lm

new_newton:	new_newton.c
	$(CC) $(CFLAGS) -o $@ new_newton.c -lm

