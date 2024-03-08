import sys
import re
import numpy as np
import subprocess


# Specify the path to the input file, and the name of the output file
input_path = 'input_psi4.tmp'
output_path= 'output_psi4.tmp'

# Set the number of threads to be used
nthreads= int(sys.argv[1])

# Set the psi4 executable path
psi4_path=sys.argv[2]

# Dictionary mapping element symbols to atomic numbers
symbol_to_atomic_number = {
    'H': 1, 'He': 2, 'Li': 3, 'Be': 4, 'B': 5, 'C': 6, 'N': 7, 'O': 8,
    'F': 9, 'Ne': 10, 'Na': 11, 'Mg': 12, 'Al': 13, 'Si': 14, 'P': 15,
    'S': 16, 'Cl': 17, 'Ar': 18,
}

# Open the file and read its content
with open(input_path, 'r') as file:
    input_text = file.read()

# Use regular expressions to find the content between "mol{" and "}"
match = re.search(r'molecule mol {(.*?)}', input_text, re.DOTALL)
if match:
    mol_content = match.group(1).strip()

    # Split the content into lines
    lines = mol_content.split('\n')

    # Separate symbols and coordinates into two arrays
    species = []
    coordinates = []
    for line in lines:
        parts = line.split()
        if len(parts) >= 4:
            symbol = parts[0]
            atomic_number = symbol_to_atomic_number.get(symbol, -1)  # -1 if not found
            species.append(atomic_number)
            coordinates.append([float(x) for x in parts[1:]])
else:
    print("Molecule content not found.")

# Run psi4 as a Bash command
bash_command = [psi4_path+" -i "+input_path+" -o "+output_path+" -n "+str(nthreads)]

# Run the Bash command
try:
    subprocess.run(bash_command, shell=True, check=True)
except subprocess.CalledProcessError as e:
    print(f"Error executing the command: {e}")

# Get natoms for reading
natoms = len(coordinates)

# Open the file and read its content
with open(output_path, 'r') as file:
    content = file.readlines()

# Search MBIS charges
pattern = "MBIS Charges: (a.u.)"
# Find the line number where the pattern is located
start_index = next((i for i, line in enumerate(content) if pattern in line), None)
start_index = start_index+1 # To skip the line containing the columns description

# Read the natoms lines following the pattern (and skiping columns descriptions)
charges=[]
if start_index is not None:
    lines_after_pattern = content[start_index + 1 : start_index + 1 + natoms]
    # Print the lines
    for line in lines_after_pattern:
        splited_line=line.split()
        charges.append(float(splited_line[-1]))
else:
    print(f"Pattern '{pattern}' not found in the file.")
# Save to numpy array
charges_PSI4=np.array(charges)

# Search energy
pattern = "Total Energy ="
# Find the line number where the pattern is located
start_index = next((i for i, line in enumerate(content) if pattern in line), None)
if start_index is not None:
    line_with_pattern = content[start_index]
    splited_line=line_with_pattern.split()
    energy_PSI4=float(splited_line[-1])
else:
    print(f"Pattern '{pattern}' not found in the file.")


# Save results to files
np.savetxt("q_psi4.tmp",charges_PSI4)
np.savetxt("e_psi4.tmp",[energy_PSI4])
