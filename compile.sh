#!/usr/bin/env bash

cmake .

if [ ! -d "cmake-build-debug" ]; then
  mkdir "cmake-build-debug"
fi

cmake --build ./cmake-build-debug/  --target pr_graphs -- -j 8

make

rm -R "cmake-build-debug"

if [ -d "CMakeFiles" ]; then
  rm -R "CMakeFiles"
fi

if [ -f "CMakeCache.txt" ]; then
  rm "CMakeCache.txt"
fi

if [ -f "cmake_install.cmake" ]; then
  rm "cmake_install.cmake"
fi

if [ -f "Makefile" ]; then
  rm "Makefile"
fi
