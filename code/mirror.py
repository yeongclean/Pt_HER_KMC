from ase.io import read, write
from ase import Atoms
import numpy as np
import copy

def make_symmetric_adsorption(poscar_path, output_path, tol=2.5):
    atoms = read(poscar_path, format='vasp')

    # Separate slab atoms and adsorbates based on z-coordinates
    z_coords = atoms.positions[:, 2]
    z_mid = np.mean([min(z_coords), max(z_coords)])

    top_atoms = [atom for atom in atoms if atom.position[2] > z_mid + tol]
    bottom_atoms = [atom for atom in atoms if atom.position[2] < z_mid - tol]
    slab_atoms = [atom for atom in atoms if z_mid - tol <= atom.position[2] <= z_mid + tol]

    # Mirror top adsorbates to bottom
    copied_atoms = []
    for atom in top_atoms:
        mirrored_pos = atom.position.copy()
        mirrored_pos[2] = 2 * z_mid - atom.position[2]
        new_atom = copy.deepcopy(atom)
        new_atom.position = mirrored_pos
        copied_atoms.append(new_atom)

    all_atoms = Atoms(
        symbols=[atom.symbol for atom in slab_atoms + top_atoms + copied_atoms],
        positions=np.array([atom.position for atom in slab_atoms + top_atoms + copied_atoms]),
        cell=atoms.cell,
        pbc=atoms.pbc
    )

    write(output_path, all_atoms, format='vasp')
    print(f"{output_path}")
    

make_symmetric_adsorption("POSCAR_100_H", "POSCAR_symmetric", tol=2.5)
