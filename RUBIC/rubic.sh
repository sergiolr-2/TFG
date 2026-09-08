
#!/bin/bash
#
# SBATCH --job-name = RUBIC
# SBATCH --output =%x_%j.out
# SBATCH -N 1
# SBATCH -c 32

module load gnu12/12.2.1
module load openmpi4/4.1.4

nprocs="$1"

# Compile the C program and create an executable file
#mpicc "RUBIC_pararell.c" -o "RUBIC_pararell" -lm

# Check if compilation was successful
# if [ $? -eq 0 ]; then
#     echo "Compilation successful."
# else
#     echo "Compilation failed."
# fi

mpirun -n "$nprocs" ./RUBIC_pararell input.txt output.txt 2 2 1 1 1

# Sample execution: mpirun -n <nprocs> RUBIC inputdata.txt rubicoutput.txt 2 2 1 1 1
# In RUBIC.c code file:
# executable = argv[0]
# input_file = argv[1]
# ouptut_file = argv[2]
# mnr = argv[3]
# mnc = argv[4]
# thr = argv[5] a cada valor de la matriz original si es menor que thr la matriz binaria sera 0 si es mayor o igual sera 1 (pasamos de valores continuos a discretos/ binarios)
# rk = argv[6]
# ck = argv[7]