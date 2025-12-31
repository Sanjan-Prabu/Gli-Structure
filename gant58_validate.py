from rdkit import Chem
from rdkit.Chem import AllChem
from vina import Vina
import numpy as np
import subprocess
import os

# --- FUNCTION: Pose Analysis ---
def analyze_docked_pose(pdbqt_file):
    coords = []
    if not os.path.exists(pdbqt_file): return None
    with open(pdbqt_file, 'r') as f:
        in_model_1 = False
        for line in f:
            if line.startswith('MODEL 1'): in_model_1 = True
            elif line.startswith('ENDMDL'): break
            if in_model_1 and line.startswith('ATOM') and not ('H' in line.split()[2]):
                try:
                    coords.append(np.array([float(line[30:38]), float(line[38:46]), float(line[46:54])]))
                except ValueError: continue
    return np.mean(coords, axis=0) if coords else None

# --------------------------------------------------------------------------
# MAIN DOCKING WORKFLOW
# --------------------------------------------------------------------------

print("="*70)
print("GANT58 - GLI1 ZF1-ZF3 DOCKING (CORRECTED SITE)")
print("="*70)

receptor_file = '2gli_receptor.pdbqt'
if not os.path.exists(receptor_file):
    print(f"❌ Error: {receptor_file} not found.")
    exit()

# GANT58 structure
gant58_smiles = "c1cnccc1-c1c(-c2ccncc2)[s]c(-c2ccncc2)c1-c1ccncc1"

print("\n[1/3] Preparing GANT58 Structure...")
mol = Chem.MolFromSmiles(gant58_smiles)
mol = Chem.AddHs(mol)
AllChem.EmbedMolecule(mol, randomSeed=42)
AllChem.MMFFOptimizeMolecule(mol)

Chem.MolToPDBFile(mol, 'temp_ligand.pdb')
subprocess.run("obabel temp_ligand.pdb -O temp_ligand.pdbqt", shell=True, capture_output=True)

# [2/3] Docking Setup
print("[2/3] Initializing Vina (Broad ZF1-ZF3 Search)...")
v = Vina(sf_name='vina')
v.set_receptor(receptor_file)

# NEW COORDINATES: Centered to encompass the ZF1-ZF2-ZF3 junction
# This moves the box slightly "up" the protein chain compared to GANT61's narrow site
CORRECTED_CENTER = [-27.0, -5.0, 2.5] 
v.compute_vina_maps(center=CORRECTED_CENTER, box_size=[22, 22, 22]) # Larger box for larger molecule

# [3/3] Execute Docking
print("[3/3] Docking GANT58...")
v.set_ligand_from_file('temp_ligand.pdbqt')
v.dock(exhaustiveness=128, n_poses=1) # Increased exhaustiveness for more accurate search

energy = v.energies(n_poses=1)[0][0]
v.write_poses('gant58_target_docked.pdbqt', n_poses=1, overwrite=True)

print("-" * 70)
print(f"DONE. Binding Affinity: {energy:.3f} kcal/mol")
print(f"Site: ZF1-ZF3 Cluster (Center: {CORRECTED_CENTER})")

# Clean up
for f in ['temp_ligand.pdb', 'temp_ligand.pdbqt']:
    if os.path.exists(f): os.remove(f)

print("Output: gant58_target_docked.pdbqt")
print("=" * 70)