#!/usr/bin/env bash

module load openmpi

if [ ! -d "out" ]; then
  mkdir "out"
fi

mpiCC -std=c++11 ./src/main.cpp -o ./out/main-pcss.mpi

echo "Compiled."