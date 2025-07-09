from ase.io import read, write
from ase import Atom
import numpy as np
from ase.constraints import FixAtoms

# read the POSCAR file
atoms = read("POSCAR_111_H", format="vasp")

# seperate slab and adsorbates
n_pt = 36 # number of Pt atoms in the slab
slab = atoms[:n_pt]
adsorbates = atoms[n_pt:]

# fix the bottom of the slab
z_positions = slab.positions[:, 2]
z_threshold = np.median(z_positions) 
fixed_slab_indices = [i for i, z in enumerate(z_positions) if z <= z_threshold]

# calcuate the center of mass of the adsorbates
z_com_ads = adsorbates.get_center_of_mass()[2]
z_top_slab = np.max(slab.positions[:, 2])
z_bottom_slab = np.min(slab.positions[:, 2])
dz = z_com_ads - z_top_slab

# list to hold mirrored adsorbates
mirrored_adsorbates = []

# generate mirrored adsorbates
for atom in adsorbates:
    x, y, z = atom.position
    z_shifted = z_bottom_slab - dz + (z - z_com_ads)
    mirrored_atom = Atom(symbol=atom.symbol, position=(x, y, z_shifted))
    mirrored_adsorbates.append(mirrored_atom)

# add mirrored adsorbates to the original atoms list
for atom in mirrored_adsorbates:
    atoms.append(atom)

# indeces for the mirrored adsorbates
n_total = len(atoms)
n_mirrored = len(mirrored_adsorbates)
mirrored_indices = list(range(n_total - n_mirrored, n_total))

# atom fixation
fix_indices = fixed_slab_indices + mirrored_indices
atoms.set_constraint(FixAtoms(indices=fix_indices))

# save the new structure to a POSCAR file
write("111H_POSCAR_sym", atoms, format="vasp", direct=False)
