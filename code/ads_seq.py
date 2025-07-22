import os
from ase.io import read, write
from ase.constraints import FixAtoms
from gg.predefined_sites import SurfaceSites
from gg.modifiers import Add, ModifierAdder

# 1) read optimized slab
slab = read("POSCAR_111", format="vasp")

# 2) Define surface site 
max_coord = {"Pt": 12, "O": 4, "H": 2}
SS = SurfaceSites(
    max_coord=max_coord,
    max_bond_ratio=1.2,  
    com=0.5               #center of mass
)

# 3) H Modifier
add_H1 = Add(
    SS,
    ads="H",            # species
    surf_coord=[1,2,3],     
    ads_id=["H"],       # ads tag
    surf_sym=["Pt"],    
    print_movie=True    
)

# 4) OH Modifier
add_H2 = Add(
    SS,
    ads="H",
    surf_coord=[1,2,3],
    ads_id=["H"],
    surf_sym=["Pt"],
    print_movie=True
)

# 5) Sequentially modified for unique site
add_2H = ModifierAdder(
    [add_H1, add_H2],
    print_movie=True,
    unique=True
)

modified_atoms_list = add_2H.get_modified_atoms(slab)

# 6) output directory
out_dir = "111-2H"
os.makedirs(out_dir, exist_ok=True)


for i, atoms_i in enumerate(modified_atoms_list, start=1):
    filename = os.path.join(out_dir, f"POSCAR_mod_{i}")
    write(filename, atoms_i, format="vasp")
    print(f"Saved: {filename}")
    
    
# 8) Save trajectory file
traj_path = os.path.join(out_dir, "all_structures.traj")
write(traj_path, modified_atoms_list) 
print(f"Saved trajectory: {traj_path}")
