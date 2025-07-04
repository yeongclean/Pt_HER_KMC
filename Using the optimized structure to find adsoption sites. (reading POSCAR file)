from ase.io import read, write
from ase.constraints import FixAtoms

from gg.modifiers import Add
from gg.predefined_sites import FlexibleSites

# 1) read optimized slab
slab = read("POSCAR_slab", format="vasp")

# 2) fix bottom of the slab
z_positions = slab.get_positions()[:,2]
z_threshold = sorted(z_positions)[int(len(z_positions)/2)]
fixed_indices = [i for i,z in enumerate(z_positions) if z < z_threshold]

slab.set_constraint(FixAtoms(indices=fixed_indices))

# 3) FlexibleSites + Add modifier : OH adsorption
FS = FlexibleSites(constraints=True,max_bond_ratio=1.2)
add_OH = Add(
    FS,
    ads="OH",
    surf_coord=[1,2,3],
    ads_id=["O"],
    surf_sym=["Pt"],
    print_movie=True,
    unique=True
)

structures = add_OH.get_modified_atoms(slab)

# 4) write POSCAR files for each structure
for i, struct in enumerate(structures):
    write(f"POSCAR_{i}", struct, format="vasp", direct=True)
