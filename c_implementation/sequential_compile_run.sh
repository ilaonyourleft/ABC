#!/bin/bash

cd src/
gcc -o ABC_sequential_lib.o -c ABC_sequential_lib.c
gcc -o ABC_sequential ABC_sequential.c ABC_sequential_lib.o
./ABC_sequential