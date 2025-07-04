from ase import Atoms
from ase.build import fcc111, add_adsorbate
from ase.visualize import view
from ase.constraints import FixAtoms
from ase.io import write
from ase.visualize import view

# Pt(111), 3x2, 4 layers
slab = fcc111('Pt', size=(3,2,4), vacuum=20.0)

# z-axis centering
slab.center(axis=2)

# fix bottom
mask = [atom.tag > 2 for atom in slab]
slab.set_constraint(FixAtoms(mask=mask))

# wrtie POSCAR
write("POSCAR", slab, direct=True)
