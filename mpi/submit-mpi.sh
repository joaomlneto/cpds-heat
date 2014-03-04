#!/bin/bash
	# @ job_name		= heatmpi
	# @ partition		= debug
	# @ initialdir		= .
	# @ output		= heatmpi.%j.out
	# @ error		= heatmpi.%j.err
	# @ total_tasks		= 4
	# @ cpus_per_task	= 1
	# @ tasks_per_node	= 4
	# @ wall_clock_limit	= 00:02:00

prog=heatmpi
procs=4

srun -n $procs ./$prog test.dat
