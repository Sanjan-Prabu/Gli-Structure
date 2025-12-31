#!/usr/bin/env python3
"""
Calculate optimal grid box size using radius of gyration method.
Based on: Feinstein & Brylinski (2015) J Cheminform 7:18
Optimal box size = 2.9 × radius of gyration
"""

import numpy as np
from rdkit import Chem
from rdkit.Chem import AllChem

def calculate_radius_of_gyration(coords):
    """
    Calculate radius of gyration for a molecule.
    Rg = sqrt(mean(distance^2 from centroid))
    """
    # Calculate geometric center
    centroid = np.mean(coords, axis=0)
    
    # Calculate squared distances from centroid
    squared_distances = np.sum((coords - centroid)**2, axis=1)
    
    # Radius of gyration
    rg = np.sqrt(np.mean(squared_distances))
    
    return rg, centroid

def get_heavy_atom_coords(mol):
    """Extract coordinates of non-hydrogen atoms."""
    coords = []
    conf = mol.GetConformer()
    
    for atom in mol.GetAtoms():
        if atom.GetAtomicNum() != 1:  # Skip hydrogens
            pos = conf.GetAtomPosition(atom.GetIdx())
            coords.append([pos.x, pos.y, pos.z])
    
    return np.array(coords)

def calculate_optimal_box_size(pdbqt_file, multiplier=2.9):
    """
    Calculate optimal box size from PDBQT file.
    
    Parameters:
    -----------
    pdbqt_file : str
        Path to ligand PDBQT file
    multiplier : float
        Box size multiplier (default 2.9 from Feinstein & Brylinski)
    
    Returns:
    --------
    box_size : float
        Optimal box size in Angstroms
    center : array
        Geometric center of ligand
    """
    coords = []
    
    # Parse PDBQT file
    with open(pdbqt_file, 'r') as f:
        for line in f:
            if line.startswith('ATOM') or line.startswith('HETATM'):
                # Skip hydrogen atoms
                if line[13] == 'H':
                    continue
                
                # Extract coordinates (columns 31-38, 39-46, 47-54)
                x = float(line[30:38].strip())
                y = float(line[38:46].strip())
                z = float(line[46:54].strip())
                coords.append([x, y, z])
    
    if not coords:
        raise ValueError("No heavy atoms found in PDBQT file")
    
    coords = np.array(coords)
    
    # Calculate radius of gyration
    rg, centroid = calculate_radius_of_gyration(coords)
    
    # Calculate optimal box size
    box_size = multiplier * rg
    
    return box_size, centroid, rg

# Main execution
if __name__ == "__main__":
    import sys
    
    # File to analyze
    if len(sys.argv) > 1:
        ligand_file = sys.argv[1]
    else:
        ligand_file = 'gant61.pdbqt'
    
    print("=" * 70)
    print("OPTIMAL GRID BOX CALCULATOR")
    print("Based on: Feinstein & Brylinski (2015)")
    print("=" * 70)
    
    try:
        # Calculate for docked GANT61 (best pose)
        box_size, center, rg = calculate_optimal_box_size('gant61_docked.pdbqt', multiplier=2.9)
        
        print(f"\n📊 Analysis of: {ligand_file}")
        print(f"\n   Radius of Gyration (Rg): {rg:.3f} Å")
        print(f"   Optimal Box Size (2.9 × Rg): {box_size:.3f} Å")
        print(f"   Rounded Box Size: {int(np.ceil(box_size))} Å")
        
        print(f"\n📍 Ligand Geometric Center:")
        print(f"   X: {center[0]:8.3f} Å")
        print(f"   Y: {center[1]:8.3f} Å")
        print(f"   Z: {center[2]:8.3f} Å")
        
        # Round to nearest odd number (for AutoDock Vina symmetry)
        rounded_size = int(np.ceil(box_size))
        if rounded_size % 2 == 0:
            rounded_size += 1
        
        print(f"\n✅ RECOMMENDED GRID BOX FOR AUTODOCK VINA:")
        print(f"   center_x = {center[0]:.1f}")
        print(f"   center_y = {center[1]:.1f}")
        print(f"   center_z = {center[2]:.1f}")
        print(f"   size_x = {rounded_size}")
        print(f"   size_y = {rounded_size}")
        print(f"   size_z = {rounded_size}")
        
        print(f"\n📝 Justification:")
        print(f"   This box size is {box_size/25:.2f}× your current 25 Å box")
        print(f"   Multiplier of 2.9 is empirically optimal (Feinstein 2015)")
        print(f"   Balances: sufficient sampling + computational efficiency")
        
        # Compare to your current box
        print(f"\n🔍 Comparison to Your Current Setup:")
        current_center = np.array([-32.6, -5.7, -0.6])
        distance = np.linalg.norm(center - current_center)
        print(f"   Distance from your box center: {distance:.2f} Å")
        
        if distance < 5:
            print(f"   ✅ Centers are close - your box position is good!")
        else:
            print(f"   ⚠️  Consider recentering to docked GANT61 position")
        
    except FileNotFoundError:
        print(f"\n❌ Error: Could not find {ligand_file}")
        print("   Make sure the file exists in the current directory")
    except Exception as e:
        print(f"\n❌ Error: {e}")
    
    print("\n" + "=" * 70)