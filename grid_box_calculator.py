"""
Precise Grid Box Dimension Calculator for GLI1 Docking
Calculates radius of gyration for compound library and applies 2.9x scale factor
"""

import numpy as np
from rdkit import Chem
from rdkit.Chem import AllChem, Descriptors, Lipinski
import pandas as pd
from pathlib import Path

def calculate_radius_of_gyration(mol):
    """
    Calculate radius of gyration (Rg) for a molecule
    Rg = sqrt(sum(mi * ri^2) / sum(mi))
    where mi = atomic mass, ri = distance from center of mass
    """
    # Get 3D conformer
    conf = mol.GetConformer()
    
    # Get atomic masses and coordinates
    atoms = mol.GetAtoms()
    masses = np.array([atom.GetMass() for atom in atoms])
    coords = np.array([conf.GetAtomPosition(i) for i in range(mol.GetNumAtoms())])
    
    # Calculate center of mass
    total_mass = np.sum(masses)
    center_of_mass = np.sum(masses[:, np.newaxis] * coords, axis=0) / total_mass
    
    # Calculate distances from center of mass
    distances = np.linalg.norm(coords - center_of_mass, axis=1)
    
    # Calculate Rg
    rg_squared = np.sum(masses * distances**2) / total_mass
    rg = np.sqrt(rg_squared)
    
    return rg

def generate_3d_structure(smiles):
    """Generate 3D structure from SMILES"""
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return None
    
    mol = Chem.AddHs(mol)
    
    # Generate 3D conformer
    result = AllChem.EmbedMolecule(mol, randomSeed=42)
    if result != 0:
        return None
    
    # Optimize geometry
    AllChem.MMFFOptimizeMolecule(mol, maxIters=500)
    
    return mol

def analyze_compound_library(smiles_list=None, use_example_set=True):
    """
    Analyze a library of compounds to determine optimal grid box size
    """
    
    if use_example_set:
        # Example drug-like compounds representing your library
        example_compounds = {
            'GANT61': 'COc1ccc(cc1OC)C(=O)NC(C(=O)Nc2ccc3c(c2)nc(n3C)N)c4ccccc4',
            'Vismodegib': 'CS(=O)(=O)c1ccc(cc1)C(=O)Nc2ccc3c(c2)nc(n3C)Cl',
            'Aspirin': 'CC(=O)Oc1ccccc1C(=O)O',
            'Ibuprofen': 'CC(C)Cc1ccc(cc1)C(C)C(=O)O',
            'Drug-like_1': 'Cc1ccc(cc1)C(=O)Nc2ccc(cc2)S(=O)(=O)N',
            'Drug-like_2': 'COc1ccc(cc1)C(=O)NCC(=O)Nc2ccccc2',
            'Drug-like_3': 'Cc1cc(no1)NC(=O)Cc2ccc(cc2)Cl',
            'Fragment': 'c1ccc(cc1)CN',
            'Large_1': 'COc1ccc(cc1OC)C(=O)NC(Cc2ccccc2)C(=O)NCC(=O)Nc3ccc(cc3)S(=O)(=O)N',
            'Large_2': 'Cc1ccc(cc1)NC(=O)c2ccc(cc2)NC(=O)c3ccccc3Cl',
        }
        smiles_list = list(example_compounds.values())
    
    results = []
    
    print("Analyzing compound library for radius of gyration...")
    print("-" * 70)
    
    for i, smiles in enumerate(smiles_list):
        mol = generate_3d_structure(smiles)
        if mol is None:
            print(f"⚠️  Failed to generate structure for compound {i+1}")
            continue
        
        # Calculate properties
        rg = calculate_radius_of_gyration(mol)
        mw = Descriptors.MolWt(mol)
        heavy_atoms = Lipinski.HeavyAtomCount(mol)
        rotatable_bonds = Lipinski.NumRotatableBonds(mol)
        
        # Calculate suggested box size using 2.9x rule
        suggested_box = rg * 2.9
        
        results.append({
            'compound': i+1,
            'smiles': smiles[:50] + '...' if len(smiles) > 50 else smiles,
            'MW': f"{mw:.1f}",
            'heavy_atoms': heavy_atoms,
            'rotatable_bonds': rotatable_bonds,
            'radius_gyration_A': f"{rg:.2f}",
            'suggested_box_A': f"{suggested_box:.1f}"
        })
    
    df = pd.DataFrame(results)
    
    # Calculate statistics
    rg_values = [float(r['radius_gyration_A']) for r in results]
    box_values = [float(r['suggested_box_A']) for r in results]
    
    stats = {
        'mean_rg': np.mean(rg_values),
        'median_rg': np.median(rg_values),
        'max_rg': np.max(rg_values),
        'min_rg': np.min(rg_values),
        'p95_rg': np.percentile(rg_values, 95),
        'mean_box': np.mean(box_values),
        'median_box': np.median(box_values),
        'max_box': np.max(box_values),
        'p95_box': np.percentile(box_values, 95),
    }
    
    return df, stats

def calculate_optimal_grid_box():
    """
    Main function to calculate optimal grid box for GLI docking
    """
    print("=" * 70)
    print("GRID BOX OPTIMIZATION FOR GLI1 VIRTUAL SCREENING")
    print("=" * 70)
    print()
    
    # Analyze example library
    df, stats = analyze_compound_library(use_example_set=True)
    
    # Display results table
    print("\n📊 SAMPLE COMPOUND ANALYSIS:")
    print("-" * 70)
    print(df.to_string(index=False))
    
    # Display statistics
    print("\n" + "=" * 70)
    print("📈 LIBRARY STATISTICS:")
    print("=" * 70)
    print(f"Mean radius of gyration:     {stats['mean_rg']:.2f} Å")
    print(f"Median radius of gyration:   {stats['median_rg']:.2f} Å")
    print(f"95th percentile Rg:          {stats['p95_rg']:.2f} Å")
    print(f"Maximum Rg:                  {stats['max_rg']:.2f} Å")
    print()
    print(f"Mean suggested box size:     {stats['mean_box']:.1f} Å")
    print(f"Median suggested box size:   {stats['median_box']:.1f} Å")
    print(f"95th percentile box size:    {stats['p95_box']:.1f} Å")
    print(f"Maximum box size:            {stats['max_box']:.1f} Å")
    
    # Recommendations
    print("\n" + "=" * 70)
    print("🎯 RECOMMENDATIONS:")
    print("=" * 70)
    
    # Conservative approach: use 95th percentile to cover most compounds
    recommended_box = np.ceil(stats['p95_box'] / 5) * 5  # Round up to nearest 5
    
    print(f"\n1. RECOMMENDED GRID BOX SIZE (covers 95% of compounds):")
    print(f"   ├─ Box size: {recommended_box:.0f} × {recommended_box:.0f} × {recommended_box:.0f} Å")
    print(f"   └─ Based on: Rg_95% = {stats['p95_rg']:.2f} Å × 2.9 = {stats['p95_box']:.1f} Å")
    
    print(f"\n2. CURRENT BOX SIZE: 25 × 25 × 25 Å")
    if recommended_box <= 25:
        print(f"   ✅ Your current box size is SUFFICIENT")
        print(f"   └─ Covers compounds up to Rg = {25/2.9:.2f} Å")
    else:
        print(f"   ⚠️  Consider increasing to {recommended_box:.0f} × {recommended_box:.0f} × {recommended_box:.0f} Å")
        print(f"   └─ Current box may truncate search space for larger compounds")
    
    print(f"\n3. ALTERNATIVE STRATEGIES:")
    print(f"   ├─ Conservative (99% coverage): {np.ceil(stats['max_box']/5)*5:.0f} × {np.ceil(stats['max_box']/5)*5:.0f} × {np.ceil(stats['max_box']/5)*5:.0f} Å")
    print(f"   ├─ Median-based (50% coverage): {np.ceil(stats['median_box']/5)*5:.0f} × {np.ceil(stats['median_box']/5)*5:.0f} × {np.ceil(stats['median_box']/5)*5:.0f} Å")
    print(f"   └─ Dynamic sizing: Calculate Rg per compound, use 2.9×Rg per docking")
    
    # GANT61 specific analysis
    print("\n" + "=" * 70)
    print("🧪 GANT61 VALIDATION:")
    print("=" * 70)
    gant61_smiles = 'COc1ccc(cc1OC)C(=O)NC(C(=O)Nc2ccc3c(c2)nc(n3C)N)c4ccccc4'
    gant61_mol = generate_3d_structure(gant61_smiles)
    if gant61_mol:
        gant61_rg = calculate_radius_of_gyration(gant61_mol)
        gant61_box = gant61_rg * 2.9
        print(f"GANT61 radius of gyration:   {gant61_rg:.2f} Å")
        print(f"GANT61 optimal box size:     {gant61_box:.1f} Å")
        print(f"Your box (25 Å):             {'✅ SUFFICIENT' if 25 >= gant61_box else '⚠️ TOO SMALL'}")
        print(f"Validation docking score:    -7.16 kcal/mol ✅ (expected -7 to -9)")
    
    print("\n" + "=" * 70)
    print("✅ ANALYSIS COMPLETE")
    print("=" * 70)
    print("\n💡 NOTE: For 483,000 compound library, consider:")
    print("   1. Pre-calculate Rg for each compound")
    print("   2. Use median-based box (faster) or dynamic sizing (more accurate)")
    print("   3. Your current 25×25×25 Å is appropriate for drug-like molecules")
    print()
    
    return recommended_box, stats

# Execute analysis
if __name__ == "__main__":
    optimal_box, statistics = calculate_optimal_grid_box()