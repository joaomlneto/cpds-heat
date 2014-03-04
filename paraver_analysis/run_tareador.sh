#!/bin/bash

USAGE="\n USAGE: run_tareador.sh PROG inputfile\n
        PROG       -> program name\n
	inputfile  -> input file with run arguments\n"

if (test $# -lt 1 || test $# -gt 2)
then
        echo -e $USAGE
        exit 0
fi

rm -rf .tareador.$1
mkdir .tareador.$1

echo -e "START TASK DEPENDENCY GRAPH GENERATION\n"

echo -e "Dynamic instrumentation and execution: starting (it may take a while ... please be pacient :)\n"

echo -e "valgrind --tool=tareador --trace-folder=.tareador.$1 --task-usage-verbosity=task-instance-usage ./$1 $2 >run_tareador.out 2>run_tareador.err\n"
valgrind --tool=tareador --trace-folder=.tareador.$1 --objects-perspective-report=yes --task-usage-verbosity=no-reports ./$1 $2 >run_tareador.out 2>run_tareador.err

echo -e "Dynamic instrumentation and execution: finalized\n"

echo -e "Building the task dependency graph: starting\n"

cd .tareador.$1
generate_dependency_graph.py trf_superscalar0
mv dep_graph_trf_superscalar0 ../dep_graph_$1.xdot 
cd - >& /dev/null

echo -e "Building the task dependency graph: finalized\n"

echo -e "Visualizing the task dependency graph (xdot.py dep_graph_$1.xdot)\n"
xdot.py dep_graph_$1.xdot & >& /dev/null

unset LD_PRELOAD
echo -e "END TASK DEPENDENCY GRAPH GENERATION\n"
echo -e "Note: Look for and read run_tareador.out and  " 
echo -e "      run_tareador.err files if the graph has "
echo -e "      not been generated\n"
