#!/usr/bin/env bash

module load openmpi

n=${1:-10}

if [ ! -f "./out/main-pcss.mpi" ]; then
  bash "compilePCSS.sh"
fi

salloc -N ${n} mpirun -n ${n} ./out/main-pcss.mpi
