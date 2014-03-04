CC =  gcc 
MPICC = mpicc
CFLAGS = -O3 -Wall -std=c99

INCLE   = -I${EXTRAE_HOME}/include/ -I.
LIBSE   = -L$(EXTRAE_HOME)/lib -lmpitrace -lxml2 

ALL	= heatmpi 
all: $(ALL)

misc.o: misc.c
	$(CC) -c $(CFLAGS) $< -o $@

heatmpi: heat-mpi.c solver-mpi.c misc.o
	$(MPICC) $(CFLAGS) -o $@ $+ -lm

clean:
	rm -fr $(ALL) *.o *.mpi *~ *.ppm *.ps *.txt *.err *.out *.prv *.pcf *.row TRACE.mpits set-0 *.sym
