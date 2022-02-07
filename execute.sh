#!/bin/bash

# Script to run a simulation

echo "Compiling the source code ... "
make impurities_hpc
echo "... Done"

echo "./impurities_hpc $1 $2"
./impurities_hpc $1 $2
