# Intel compiler
CC =  icc
CFLAGS = -O3 -fno-alias -xhost

MPICC = mpicc

all: heat

heat : heat.o input.o misc.o timing.o relax_gauss.o relax_jacobi.o
	$(CC) $(CFLAGS) -o $@ $+ -lm $(PAPI_LIB)

%.o : %.c heat.h timing.h input.h
	$(CC) $(CFLAGS) $(PAPI_INC) -c -o $@ $< 

clean:
	rm -f *.o heat *~ *.ppm *.annot

remake : clean all
