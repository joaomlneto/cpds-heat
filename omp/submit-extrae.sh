#!/bin/bash
	# @ job_name		= heatomp.extrae
	# @ partition		= debug
	# @ initialdir		= .
	# @ output		= heatomp.extrae.%j.out
	# @ error		= heatomp.extrae.%j.err
	# @ total_tasks		= 1
	# @ cpus_per_task	= 12
	# @ wall_clock_limit	= 00:02:00

executable=heatomp

export EXTRAE_CONFIG_FILE=extrae.xml
source ${EXTRAE_HOME}/etc/extrae.sh
export LD_PRELOAD=${EXTRAE_HOME}/lib/libmpitrace.so

n_threads=8
export OMP_NUM_THREADS=$n_threads
srun ${EXTRAE_HOME}/bin/extrae -v ./heatomp test.dat
