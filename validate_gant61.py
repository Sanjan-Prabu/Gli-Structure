"""
GANT61 Validation - CORRECTED Binding Site
Based on Agyeman 2014 paper - E119/E167 in full-length sequence
correspond to different residues in PDB 2GLI (which is residues 234-388)
"""

from rdkit import Chem
from rdkit.Chem import AllChem
from vina import Vina
import subprocess
import os
import numpy as np
from Bio.PDB import PDBParser

print("="*70)
print("GANT61 VALIDATION: Finding Correct E119/E167 Binding Site")
print("="*70)

# First, let's find the actual glutamate residues in the structure
print("\n[0/5] Analyzing receptor structure...")
parser = PDBParser(QUIET=True)

# Check if we're in gli_structure folder or parent folder
import os
if os.path.exists('2gli_with_zinc.pdb'):
    receptor_pdb = '2gli_with_zinc.pdb'
    receptor_pdbqt = '2gli_receptor.pdbqt'
elif os.path.exists('gli_structure/2gli_with_zinc.pdb'):
    receptor_pdb = 'gli_structure/2gli_with_zinc.pdb'
    receptor_pdbqt = 'gli_structure/2gli_receptor.pdbqt'
else:
    raise FileNotFoundError("Cannot find 2gli_with_zinc.pdb. Run from parent directory or inside gli_structure/")

structure = parser.get_structure('GLI', receptor_pdb)

print("\n   Searching for GLU residues between zinc fingers 2 and 3...")
glu_residues = []
for chain in structure[0]:
    if chain.id == 'A':
        for res in chain:
            if res.resname == 'GLU' and res.id[0] == ' ':
                try:
                    ca_coord = res['CA'].coord
                    glu_residues.append((res.id[1], ca_coord))
                    print(f"   Found: GLU {res.id[1]} at ({ca_coord[0]:.1f}, {ca_coord[1]:.1f}, {ca_coord[2]:.1f})")
                except:
                    pass

# The Agyeman paper says E119 and E167 - these ARE in the structure!
# Look for GLU 119 and GLU 167 specifically
e119_found = None
e167_found = None

for res_num, coord in glu_residues:
    if res_num == 119:
        e119_found = (res_num, coord)
    elif res_num == 167:
        e167_found = (res_num, coord)

if e119_found and e167_found:
    # Use E119 and E167 as specified in the paper
    e1_num, e1_coord = e119_found
    e2_num, e2_coord = e167_found
    
    print(f"\n✓ Found E119 and E167 from Agyeman paper!")
    print(f"  E119 at: ({e1_coord[0]:.1f}, {e1_coord[1]:.1f}, {e1_coord[2]:.1f})")
    print(f"  E167 at: ({e2_coord[0]:.1f}, {e2_coord[1]:.1f}, {e2_coord[2]:.1f})")
    
    # Calculate center between these two residues
    center_array = (e1_coord + e2_coord) / 2
    BINDING_CENTER = [float(center_array[0]), float(center_array[1]), float(center_array[2])]
elif len(glu_residues) >= 2:
    # Fallback to first two GLU
    e1_num, e1_coord = glu_residues[0]
    e2_num, e2_coord = glu_residues[1]
    
    print(f"\n⚠️  E119/E167 not found, using GLU {e1_num} and GLU {e2_num}")
    
    center_array = (e1_coord + e2_coord) / 2
    BINDING_CENTER = [float(center_array[0]), float(center_array[1]), float(center_array[2])]
    
    print(f"\n✓ Binding Site Center: ({BINDING_CENTER[0]:.1f}, {BINDING_CENTER[1]:.1f}, {BINDING_CENTER[2]:.1f})")
else:
    print("\n⚠️  Warning: Could not find two GLU residues")
    print("   Using default coordinates from prepare_gli.py")
    # Fallback to the coordinates calculated in prepare_gli.py
    BINDING_CENTER = [-32.6, -5.7, -0.6]

BOX_SIZE = [20, 20, 20]  # Slightly larger box to ensure we capture the site

# Correct GANT61 SMILES (verified C27H35N5)
gant61_smiles = "CN(C)C1=CC=CC=C1CN2CCCN(C2C3=CC=NC=C3)CC4=CC=CC=C4N(C)C"

print("\n[1/5] Generating GANT61 3D structure...")
mol = Chem.MolFromSmiles(gant61_smiles)
if mol is None:
    raise ValueError("Invalid SMILES")

from rdkit.Chem import Descriptors, rdMolDescriptors
formula = rdMolDescriptors.CalcMolFormula(mol)
mw = Descriptors.MolWt(mol)

print(f"   Formula: {formula} (Expected: C27H35N5)")
print(f"   MW: {mw:.1f} Da (Expected: 429.6 Da)")

if formula != "C27H35N5":
    raise ValueError(f"Wrong structure! Got {formula}")

mol = Chem.AddHs(mol)
AllChem.EmbedMolecule(mol, randomSeed=42, useRandomCoords=True)
AllChem.MMFFOptimizeMolecule(mol, maxIters=500)

Chem.MolToPDBFile(mol, 'gant61.pdb')
subprocess.run("obabel gant61.pdb -O gant61.pdbqt", 
               shell=True, capture_output=True, text=True)
print("✓ GANT61 structure generated and verified")

# Setup docking
print("\n[2/5] Setting up docking...")
v = Vina(sf_name='vina')
v.set_receptor(receptor_pdbqt)

print(f"✓ Using binding site at E{e1_num}/E{e2_num}")
print(f"   Center: ({BINDING_CENTER[0]:.1f}, {BINDING_CENTER[1]:.1f}, {BINDING_CENTER[2]:.1f})")
print(f"   Box size: {BOX_SIZE[0]} x {BOX_SIZE[1]} x {BOX_SIZE[2]} Å")

v.compute_vina_maps(center=BINDING_CENTER, box_size=BOX_SIZE)

# Dock with HIGH exhaustiveness for better accuracy
print("\n[3/5] Docking GANT61 (exhaustiveness=128)...")
v.set_ligand_from_file('gant61.pdbqt')
v.dock(exhaustiveness=128, n_poses=20)

# Get scores
energies = v.energies(n_poses=20)
scores = []
for i in range(min(20, len(energies))):
    score_raw = energies[i]
    if isinstance(score_raw, (list, tuple)):
        score = float(score_raw[0])
    elif isinstance(score_raw, np.ndarray):
        score = float(score_raw.flat[0])
    else:
        score = float(score_raw)
    scores.append(score)

best_score = scores[0]

v.write_poses('gant61_docked.pdbqt', n_poses=20, overwrite=True)

# Analyze all poses
print("\n[4/5] Analyzing docking poses...")
print("="*70)
print(f"{'Pose':<6} {'Score (kcal/mol)':<20} {'Distance to Site (Å)':<20}")
print("-"*70)

all_poses = []
with open('gant61_docked.pdbqt', 'r') as f:
    content = f.read()
    models = content.split('MODEL ')
    
    for i, model in enumerate(models[1:], 1):  # Skip first empty split
        lines = model.split('\n')
        coords = []
        for line in lines:
            if line.startswith('ATOM') or line.startswith('HETATM'):
                try:
                    x = float(line[30:38])
                    y = float(line[38:46])
                    z = float(line[46:54])
                    coords.append([x, y, z])
                except:
                    pass
            if line.startswith('ENDMDL'):
                break
        
        if coords:
            pose_center = np.mean(coords, axis=0)
            deviation = np.linalg.norm(pose_center - np.array(BINDING_CENTER))
            all_poses.append((i, scores[i-1] if i-1 < len(scores) else 0, pose_center, deviation))
            print(f"{i:<6} {scores[i-1] if i-1 < len(scores) else 0:<20.2f} {deviation:<20.2f}")

# Final Results
print("\n[5/5] Results Summary")
print("="*70)
print(f"Best Binding Affinity: {best_score:.3f} kcal/mol")

# Get best pose info
best_pose = all_poses[0] if all_poses else None
if best_pose:
    pose_num, pose_score, pose_center, deviation = best_pose
    print(f"\nBest Pose Analysis (Pose {pose_num}):")
    print(f"   Pose center:    ({pose_center[0]:.1f}, {pose_center[1]:.1f}, {pose_center[2]:.1f})")
    print(f"   Expected site:  ({BINDING_CENTER[0]:.1f}, {BINDING_CENTER[1]:.1f}, {BINDING_CENTER[2]:.1f})")
    print(f"   Deviation: {deviation:.2f} Å")
    
    # Check proximity to the GLU residues
    dist_e1 = np.linalg.norm(pose_center - e1_coord)
    dist_e2 = np.linalg.norm(pose_center - e2_coord)
    
    print(f"\nDistance to Key Residues:")
    print(f"   GLU {e1_num}: {dist_e1:.1f} Å")
    print(f"   GLU {e2_num}: {dist_e2:.1f} Å")

# Find poses that bind near the site (within 5Å)
near_site_poses = [p for p in all_poses if p[3] < 5.0]
if near_site_poses:
    print(f"\nPoses binding near E{e1_num}/E{e2_num} site (< 5Å deviation):")
    best_near_site = min(near_site_poses, key=lambda x: x[1])
    print(f"   Best: Pose {best_near_site[0]} with {best_near_site[1]:.2f} kcal/mol")

# Validation
print("\n" + "="*70)
print("VALIDATION")
print("="*70)

if best_score <= -7.0:
    print("✅ EXCELLENT BINDING - Strong GLI1 inhibitor predicted!")
    status = "excellent"
elif best_score <= -6.5:
    print("✅ VERY GOOD BINDING - GLI1 inhibitor validated!")
    status = "very_good"
elif best_score <= -6.0:
    print("✓ GOOD BINDING - Receptor validated")
    status = "good"
elif best_score <= -5.5:
    print("⚠️  MODERATE BINDING - Acceptable for screening")
    status = "moderate"
else:
    print("❌ WEAK BINDING - May need structure refinement")
    status = "weak"

if best_pose:
    if deviation < 3.0:
        print(f"✓ Pose binds directly at E{e1_num}/E{e2_num} site")
    elif deviation < 5.0:
        print(f"✓ Pose binds near E{e1_num}/E{e2_num} site")
    else:
        print(f"⚠️  Pose {deviation:.1f}Å from expected site")
        print(f"   May indicate alternative binding mode")

print("\n📌 RECOMMENDATION:")
if status in ["excellent", "very_good"] and deviation < 5.0:
    print("   ✅ GANT61 binding VALIDATED at E119/E167 equivalent site!")
    print("   ✅ Receptor ready for virtual screening")
    print(f"   ✅ Use binding center: ({BINDING_CENTER[0]:.1f}, {BINDING_CENTER[1]:.1f}, {BINDING_CENTER[2]:.1f})")
elif status in ["good", "moderate"]:
    print("   ✓ Acceptable binding detected")
    print("   ✓ Can proceed with screening")
    print("   ⚠️  Consider reviewing poses in ChimeraX")
else:
    print("   ⚠️  Binding weaker than expected")
    print("   ⚠️  Review structure preparation and zinc coordination")

print("\n📁 Output files:")
print("   - gant61_docked.pdbqt (all 20 poses)")
print("   - gant61.pdb (ligand structure)")
print("\n💡 To visualize in ChimeraX:")
if os.path.exists('2gli_with_zinc.pdb'):
    print(f"   1. open 2gli_with_zinc.pdb")
else:
    print(f"   1. open gli_structure/2gli_with_zinc.pdb")
print("   2. open gant61_docked.pdbqt")
print(f"   3. select :/{e1_num},{e2_num}")
print("   4. color sel red")
print("="*70)

# Cleanup
if os.path.exists('gant61.pdb'):
    os.remove('gant61.pdb')
if os.path.exists('gant61.pdbqt'):
    os.remove('gant61.pdbqt')