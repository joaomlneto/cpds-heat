#!/bin/bash
	# @ job_name		= heatomp
	# @ partition		= debug
	# @ initialdir		= .
	# @ output		= heatomp.%j.out
	# @ error		= heatomp.%j.err
	# @ total_tasks		= 1
	# @ cpus_per_task	= 12
	# @ wall_clock_limit	= 00:02:00

executable=heatomp
n_threads=1
MAX_THREADS=12
while (test $n_threads -le $MAX_THREADS)
  do
	export OMP_NUM_THREADS=$n_threads
	srun ./$executable test.dat
        n_threads=`expr $n_threads + 1`
  done
