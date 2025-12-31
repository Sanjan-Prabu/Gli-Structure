from rdkit import Chem
from rdkit.Chem import AllChem, Descriptors, rdMolDescriptors

# The SMILES you're using
test_smiles = "CN(C)C1=CC=CC=C1CN2CCCN(C2C3=CC=NC=C3)CC4=CC=CC=C4N(C)C"

print("="*70)
print("GANT61 STRUCTURE VALIDATION")
print("="*70)

# Step 1: Parse SMILES
print("\n[Step 1] Parsing SMILES...")
mol = Chem.MolFromSmiles(test_smiles)

if mol is None:
    print("❌ ERROR: Could not parse SMILES!")
    print("The SMILES string is invalid.")
    exit()
else:
    print("✓ SMILES parsed successfully")

# Step 2: Check molecular properties
print("\n[Step 2] Checking molecular properties...")
formula = rdMolDescriptors.CalcMolFormula(mol)
mw = Descriptors.MolWt(mol)
num_atoms = mol.GetNumAtoms()
num_bonds = mol.GetNumBonds()

print(f"  Molecular formula: {formula}")
print(f"  Molecular weight:  {mw:.2f} Da")
print(f"  Number of atoms:   {num_atoms}")
print(f"  Number of bonds:   {num_bonds}")

# Expected values for GANT61
expected_formula = "C27H35N5"
expected_mw = 429.6

if formula == expected_formula:
    print(f"  ✓ Formula matches expected: {expected_formula}")
else:
    print(f"  ❌ Formula MISMATCH!")
    print(f"     Expected: {expected_formula}")
    print(f"     Got:      {formula}")

if abs(mw - expected_mw) < 1.0:
    print(f"  ✓ Molecular weight matches: ~{expected_mw} Da")
else:
    print(f"  ❌ Molecular weight MISMATCH!")
    print(f"     Expected: ~{expected_mw} Da")
    print(f"     Got:      {mw:.2f} Da")

# Step 3: Check for specific atoms
print("\n[Step 3] Checking atom composition...")
atom_counts = {}
for atom in mol.GetAtoms():
    symbol = atom.GetSymbol()
    atom_counts[symbol] = atom_counts.get(symbol, 0) + 1

print(f"  Atom inventory: {atom_counts}")

# GANT61 should have: C, H, N only (NO oxygen!)
if 'O' in atom_counts:
    print(f"  ❌ CRITICAL ERROR: Found {atom_counts['O']} oxygen atoms!")
    print(f"     GANT61 should NOT contain oxygen!")
    print(f"     This indicates the wrong structure.")
else:
    print(f"  ✓ No oxygen atoms (correct)")

if atom_counts.get('N', 0) == 5:
    print(f"  ✓ Has 5 nitrogen atoms (correct)")
else:
    print(f"  ❌ Wrong number of nitrogens: {atom_counts.get('N', 0)} (expected 5)")

# Step 4: Generate 3D structure
print("\n[Step 4] Generating 3D structure...")
mol = Chem.AddHs(mol)
result = AllChem.EmbedMolecule(mol, randomSeed=42)

if result == -1:
    print("  ❌ ERROR: Could not generate 3D coordinates!")
    print("  This might indicate SMILES issues.")
else:
    print("  ✓ 3D structure generated successfully")
    
    # Optimize
    AllChem.MMFFOptimizeMolecule(mol)
    print("  ✓ Energy minimization complete")

# Step 5: Save and check PDB file
print("\n[Step 5] Saving to PDB file...")
try:
    Chem.MolToPDBFile(mol, 'diagnostic_gant61.pdb')
    print("  ✓ Saved to: diagnostic_gant61.pdb")
    
    # Read back and check for oxygen
    with open('diagnostic_gant61.pdb', 'r') as f:
        pdb_content = f.read()
        if ' O ' in pdb_content or 'ATOM      O' in pdb_content:
            print("  ❌ WARNING: PDB file contains oxygen atoms!")
        else:
            print("  ✓ PDB file has no oxygen (correct)")
            
except Exception as e:
    print(f"  ❌ ERROR saving PDB: {e}")

# Step 6: Test Open Babel conversion
print("\n[Step 6] Testing Open Babel conversion...")
import subprocess
import os

if os.path.exists('diagnostic_gant61.pdb'):
    result = subprocess.run(
        "obabel diagnostic_gant61.pdb -O diagnostic_gant61.pdbqt",
        shell=True,
        capture_output=True,
        text=True
    )
    
    if result.returncode == 0:
        print("  ✓ Open Babel conversion successful")
        
        # Check PDBQT for oxygen
        with open('diagnostic_gant61.pdbqt', 'r') as f:
            pdbqt_content = f.read()
            oxygen_lines = [line for line in pdbqt_content.split('\n') 
                          if 'ATOM' in line and ' O' in line]
            
            if oxygen_lines:
                print(f"  ❌ CRITICAL: PDBQT has {len(oxygen_lines)} oxygen atoms!")
                print("  First oxygen line:")
                print(f"     {oxygen_lines[0]}")
            else:
                print("  ✓ PDBQT has no oxygen atoms (correct)")
    else:
        print(f"  ❌ Open Babel failed: {result.stderr}")

print("\n" + "="*70)
print("DIAGNOSIS COMPLETE")
print("="*70)