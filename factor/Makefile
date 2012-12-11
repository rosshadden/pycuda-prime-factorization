CC =		gcc
CFLAGS =	-O
LDFLAGS =	-s
BINDIR =	/usr/local/bin
MANDIR =	/usr/local/man/man1

all:		factor

factor:		factor.c low_primes.h
	${CC} ${CFLAGS} factor.c ${LDFLAGS} -o factor

install:	all
	rm -f ${BINDIR}/factor
	cp factor ${BINDIR}
	rm -f ${MANDIR}/factor.1
	cp factor.1 ${MANDIR}

clean:
	rm -f factor *.o core
