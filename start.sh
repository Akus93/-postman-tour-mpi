#!/usr/bin/env bash

bash ./compile.sh

n=${1:-10}

if [ ! -f "pr_graphs" ]; then
  bash "compile.sh"
fi

mpirun -N ${n} pr_graphs
