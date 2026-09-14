import sys
import random

# MAIN
# argv[1] -> number of rows (instances)
# argv[2] -> number of columns (attributes)
# argv[3] -> proportion of 1s in the matrix (as a float)
# argv[4] -> output file

# Checking of arguments
if len(sys.argv) != 5:
    print("matrix_generator.py <instances> <attributes> <proportion of 1s> <output file>")
    sys.exit(1)

# Arguments
n_instances = int(sys.argv[1])
n_attributes = int(sys.argv[2])
proportion = float(sys.argv[3])
output_file = sys.argv[4]

# Create binary matrix
matrix = []
for _ in range(n_instances):
    instance = ""
    for _ in range(n_attributes):
        bit = 1 if random.random() < proportion else 0
        instance = f"{instance}{bit},"
    instance = f"{instance.removesuffix(",")}\n" # remove last comma and add newline
    instance = instance.removesuffix(",")
    matrix.append(instance)

# Write binary matrix to output file
with open(output_file, "w") as output:
    output.writelines(matrix)

print("Matrix generated successfully at file ", output_file)