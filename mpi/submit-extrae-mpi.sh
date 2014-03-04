#!/bin/bash
	# @ job_name		= heatmpi.extrae
	# @ partition		= debug
	# @ initialdir		= .
	# @ output		= heatmpi.extrae.%j.out
	# @ error		= heatmpi.extrae.%j.err
	# @ total_tasks		= 4
	# @ cpus_per_task	= 1
	# @ tasks_per_node	= 4
	# @ wall_clock_limit	= 00:02:00

prog=heatmpi
procs=4

export EXTRAE_CONFIG_FILE=extrae.xml
source ${EXTRAE_HOME}/etc/extrae.sh
export LD_PRELOAD=${EXTRAE_HOME}/lib/libmpitrace.so

srun -n $procs ./$prog test.dat
srun ${EXTRAE_HOME}/bin/extrae -v ./$prog test.dat

