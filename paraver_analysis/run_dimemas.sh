#!/bin/bash

USAGE="\n USAGE: run_dimemas.sh PROG n_proc \n
        PROG       -> program name\n
        n_proc  -> number of physical cpus (1,2,4,8,16,32 or 64)\n"

if (test $# -lt 2 || test $# -gt 3)
then
        echo -e $USAGE
        exit 0
fi

echo -e "START DIMEMAS SIMULATION\n"

cd .tareador.$1

echo -e "Trace merging and translation\n"
# merge the traces
merge_trf_files.py -n 1 -t superscalar
# translate to the new format
trf2trf trf_superscalar.trf trf_superscalar_new.trf

echo -e "Dimemas simulation: starting ... (it may take a while ... please be patient :)\n"
echo -e "dimemas_simulation_ncores.py trf_superscalar_new.trf $2 -t pcf_superscalar.pcf \n"
dimemas_simulation_ncores.py trf_superscalar_new.trf $2  -t pcf_superscalar.pcf 
echo -e "Dimemas simulation: finished\n"

cd ..

cp .tareador.$1/simulation_results/prv_$2cores.prv $1_sim$2cpu.prv
cp .tareador.$1/simulation_results/prv_$2cores.pcf $1_sim$2cpu.pcf
echo -e "END DIMEMAS SIMULATION\n"

echo -e "VISUALIZE THE TRACE GENERATED USING PARAVER \n"
echo -e "wxparaver $1_sim$2cpu.prv tareador.cfg \n"
