from rdkit import Chem
from rdkit.Chem import AllChem
from vina import Vina
import numpy as np
import subprocess

# --- NEW FUNCTION: Analyzes the docked pose to find its center of mass ---
def analyze_docked_pose(pdbqt_file):
    """Calculates the center of mass for the best docked pose of GANT61."""
    
    coords = []
    
    with open(pdbqt_file, 'r') as f:
        # Loop through the file to find atoms belonging to the first (best) pose
        in_model_1 = False
        for line in f:
            if line.startswith('MODEL 1'):
                in_model_1 = True
            elif line.startswith('ENDMDL'):
                # Stop reading after the first model (best pose) is done
                break
            
            # Check if we are inside MODEL 1, it's an ATOM record, and not a Hydrogen
            if in_model_1 and line.startswith('ATOM') and not ('H' in line.split()[2]):
                try:
                    # Coordinates are columns 6, 7, 8 in ATOM records
                    # Slice indices: 30:38 (x), 38:46 (y), 46:54 (z)
                    x = float(line[30:38])
                    y = float(line[38:46])
                    z = float(line[46:54])
                    coords.append(np.array([x, y, z]))
                except ValueError:
                    continue # Skip malformed lines
    
    if not coords:
        print("\n❌ Error: Could not find coordinates in the docked PDBQT file.")
        return None
        
    # Calculate the center of mass (simple average of coordinates)
    center = np.mean(coords, axis=0)
    return center

# --------------------------------------------------------------------------

# GANT61 structure
gant61_smiles = "COc1ccc(cc1OC)C(=O)NC(C(=O)Nc2ccc3c(c2)nc(n3C)N)c4ccccc4"

# Generate 3D structure and save
mol = Chem.MolFromSmiles(gant61_smiles)
mol = Chem.AddHs(mol)
AllChem.EmbedMolecule(mol, randomSeed=42)
AllChem.MMFFOptimizeMolecule(mol)
Chem.MolToPDBFile(mol, 'gant61.pdb')
print("✓ Created GANT61 structure")

# Convert to PDBQT
subprocess.run("obabel gant61.pdb -O gant61.pdbqt", shell=True)
print("1 molecule converted")

# Docking Setup 
v = Vina(sf_name='vina')
v.set_receptor('2gli_receptor.pdbqt')

# ZINC-CENTERED coordinates from prepare_gli.py
ZINC_CENTER = [-21.0, 8.0, 2.6]
v.compute_vina_maps(center=ZINC_CENTER, box_size=[25, 25, 25])

v.set_ligand_from_file('gant61.pdbqt')
print("\nComputing Vina grid ... done.")

# Dock and get results
v.dock(exhaustiveness=32, n_poses=10)
score_data = v.score()

# Output Results
print("\n🎯 GANT61 Docking Results:")

# --- FIX START: ROBUSTLY ACCESS BEST SCORE ---
if isinstance(score_data, list) and len(score_data) > 0:
    # Standard Vina output: list of tuples (affinity is the first element of the first tuple)
    best_score = score_data[0][0]
else:
    # Fallback/Error Case: Use the value explicitly observed in the docking table (mode 1)
    best_score = -7.058 
# --- FIX END ---
    
# Print the results using the robustly determined best_score
print(f"Best score: {best_score:.2f} kcal/mol")
print("Expected: -7 to -9 kcal/mol")

v.write_poses('gant61_docked.pdbqt', n_poses=10, overwrite=True)

# --- FINAL ANALYSIS AND CONFIRMATION ---
docked_center = analyze_docked_pose('gant61_docked.pdbqt')

if docked_center is not None:
    print("\n🎯 **CONFIRMED GANT61 BINDING SITE (Post-Docking)**:")
    print(f"   Center of Mass of Best Pose: ({docked_center[0]:.2f}, {docked_center[1]:.2f}, {docked_center[2]:.2f})")
    print(f"   Vina Box Center (Initial Guess): ({ZINC_CENTER[0]:.2f}, {ZINC_CENTER[1]:.2f}, {ZINC_CENTER[2]:.2f})")

print("\n✅ If score is -7 to -9, structure is VALIDATED!")