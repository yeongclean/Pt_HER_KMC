import os
from ase.io import read, write
from fairchem.core import pretrained_mlip, FAIRChemCalculator
from ase.optimize import LBFGS

# set the path
input_dir = "/home/users/yeong/work/25_s/VASPsol++/oc20/00-111/input"  #directory
output_dir = "/home/users/yeong/work/25_s/VASPsol++/oc20/00-111/output"
os.makedirs(output_dir, exist_ok=True)

# setting calculator
predictor = pretrained_mlip.get_predict_unit("uma-s-1p1", device="cuda")

# list for saving the results
results = []

# for loop
for fname in os.listdir(input_dir):
    if not fname.startswith("POSCAR"):
        continue
    in_path = os.path.join(input_dir, fname)
    atoms = read(in_path, format="vasp")

    # optimizing
    calc = FAIRChemCalculator(predictor, task_name="oc20")
    atoms.calc = calc

    opt = LBFGS(atoms, logfile=None)
    opt.run(fmax=0.05, steps=100)

    energy = atoms.get_potential_energy()
    out_path = os.path.join(output_dir, f"{fname}_opt.vasp")
    write(out_path, atoms, format="vasp")

    results.append({
        "file": out_path,
        "energy": energy,
        "name": fname
    })
    print(f"[{fname}] E = {energy:.4f} eV")

# select the best structure
best = sorted(results, key=lambda x: x["energy"])[0]
best_path = os.path.join(output_dir, "best_POSCAR")
write(best_path, read(best["file"], format="vasp"), format="vasp")

print(f"\ best structure : {best['name']}, E = {best['energy']:.4f} eV")
print(f"path: {best_path}")
