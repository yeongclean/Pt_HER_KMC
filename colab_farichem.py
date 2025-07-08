from ase.io import read, write
from ase.build import add_adsorbate, molecule
from ase.optimize import LBFGS
from fairchem.core import pretrained_mlip, FAIRChemCalculator
import numpy as np

# Load slab from POSCAR
slab = read("/content/drive/MyDrive/input_POSCAR/POSCAR_111")  #update the path

# Create adsorbate
adsorbate = molecule("H")

# Find the highest atom
z_max = max(atom.position[2] for atom in slab)
top_atoms = [atom for atom in slab if abs(atom.position[2] - z_max) < 0.5]
top_atom = top_atoms[0] 

xy = top_atom.position[:2]

# Add H on top site, adjust height as needed
add_adsorbate(slab, adsorbate,1.1,position=xy)

# Assign FAIRChemCalculator
predictor = pretrained_mlip.get_predict_unit("uma-s-1p1", device="cuda")
calc = FAIRChemCalculator(predictor, task_name="oc20")
slab.calc = calc

# Optimize
opt = LBFGS(slab)
opt.run(fmax=0.05, steps=100)

# Write final optimized structure to POSCAR
output_path = "/content/drive/MyDrive/output_POSCAR/POSCAR_111_out"
write(output_path, slab, format="vasp")

output_path
