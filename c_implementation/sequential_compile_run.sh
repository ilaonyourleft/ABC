#!/bin/bash

cd src/
gcc -o ABC_lib.o -c ABC_lib.c
gcc -o ABC_sequential ABC_sequential.c ABC_lib.o
./ABC_sequential