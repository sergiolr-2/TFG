#!/bin/bash

# Sample execution: ./P2-runwargs.sh RUBIC inputdata.txt rubicoutput.txt 1 1 1 1 1
# In RUBIC.c code file:
# executable = argv[0]
# input_file = argv[1]
# ouptut_file = argv[2]
# mnr = argv[3]
# mnc = argv[4]
# thr = argv[5] a cada valor de la matriz original si es menor que thr la matriz binaria sera 0 si es mayor o igual sera 1 (pasamos de valores continuos a discretos/ binarios)
# rk = argv[6]
# ck = argv[7]

# Check if at least 6 arguments are provided (program name + 5 arguments)
if [ "$#" -lt 6 ]; then
    echo "Usage: $0 <executable_file> <input_file> <output_file> <mnr> <mnc> <thr> <rk> <ck>"
    exit 1
fi

# Extract the program name and remove it from the argument list
program_name="$1"
shift

# Check if the program exists
if [ ! -x "$program_name" ]; then
    echo "Error: The program '$program_name' does not exist or is not executable."
    exit 1
fi

# Execute the program with the remaining arguments (4 processes)
mpirun -n 1 "$program_name" "$@"