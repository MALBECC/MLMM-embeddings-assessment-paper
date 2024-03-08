import re
import torch
import torchani
import numpy as np
import warnings
warnings.filterwarnings("ignore")

# Specify the path to the input file
input_path = 'input_ml.tmp'

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

device = torch.device('cuda' if torch.cuda.is_available() else 'cpu')
coordinates = torch.tensor([coordinates],requires_grad=True, device=device)
species = torch.tensor([species], device=device)

model = torchani.models.ANI2x(periodic_table_index=True).to(device)
energy = model((species, coordinates)).energies
derivative = torch.autograd.grad(energy.sum(), coordinates)[0]
forces = -derivative

coordinates=coordinates.detach().numpy()
coordinates=coordinates[0]
forces=forces.detach().numpy()
forces_ANI=forces[0]*627.503
energy_ANI=energy.item()

# Save results to files
np.savetxt("f_ml.tmp",forces_ANI)
np.savetxt("e_ml.tmp",[energy_ANI])
