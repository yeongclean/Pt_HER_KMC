import os
from ase.io import read, write
from ase.constraints import FixAtoms
from gg.predefined_sites import SurfaceSites
from gg.modifiers import Add, ModifierAdder
from ase.neighborlist import NeighborList, natural_cutoffs


# 1) read optimized slab
slab = read("POSCAR_111", format="vasp")

# 2) Define surface site 
max_coord = {"Pt": 12, "O": 4, "H": 2}
SS = SurfaceSites(
    max_coord=max_coord,
    max_bond_ratio=1.2,  
    com=0.5               #center of mass
)

# 3) H Modifier (for atop, bridge, and hollow sites)
add_H = Add(
    SS,
    ads="H",            # species
    surf_coord=[1,2,3],     
    ads_id=["H"],       # ads tag
    surf_sym=["Pt"],    
    print_movie=True    
)

# Generate modified structures
modified_atoms_list = add_H.get_modified_atoms(slab)

# Create ouput directory
out_dir = "ads_111H"
os.makedirs(out_dir, exist_ok=True)


for i, atoms_i in enumerate(modified_atoms_list, start=1):
    #(a) Save POSCAR file
    filename = os.path.join(out_dir, f"POSCAR_{i}")
    write(filename, atoms_i, format="vasp")
    print(f"Saved: {filename}")
    
    # b) Calculate bond counts
    nl = NeighborList(natural_cutoffs(atoms_i), self_interaction=False, bothways=True)
    nl.update(atoms_i)
    
    for idx, sym in enumerate(atoms_i.get_chemical_symbols()):
        if sym == "H":
            neighbors, offsets = nl.get_neighbors(idx)
            pt_bonds = sum(1 for n in neighbors if atoms_i[n].symbol == "Pt")
            print(f"Structure {i}: H atom at index {idx} has {pt_bonds} Pt bond(s)")
    
