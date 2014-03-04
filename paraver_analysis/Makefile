CC     =  gcc 
TAREADOR_CFLAGS = -Wall -std=c99 -O0 -g

OPT3 	= -O3
OPTG0 	= -g -O0
CFLAGS  = -Wall -std=c99
LFLAGS  = -lm 
CINCL    = -I. -I${TAREADOR_HOME}/codes/valgrind_module/tareador/lib_valgrind_annotations
CLIBS    = -L. ${TAREADOR_HOME}/codes/valgrind_module/tareador/lib_valgrind_annotations/tareador_hooks.o -lm

ALL	= heat heat-tareador 
all: $(ALL)

misc.o: misc.c
	$(CC) -c $(OPT3) $(CFLAGS) $< -o $@

heat: heat.c solver.c misc.o
	$(CC) $(OPT3) $(CFLAGS) $+ $(LFLAGS) -o $@

heat-tareador: heat-tareador.c solver-tareador.c misc.o 
	$(CC) $(OPTG0) $(CFLAGS) $(CINCL) $+ -o $@ $(CLIBS) $(LFLAGS)

clean:
	rm -fr $(ALL) .tareador.* *.o *.xdot *.ppm *.out *.err *.prv *.pcf *.row 

