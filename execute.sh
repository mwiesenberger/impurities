#!/bin/bash

# Script to run a simulation

echo "Compiling the source code ... "
make impurities_hpc
echo "... Done"

FILE="${2%.*}" # outfilename without extension
echo "./impurities_hpc $1 $2"
./impurities_hpc $1 $2 &> ${FILE}.out
