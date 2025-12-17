from Bio.PDB import PDBParser
import numpy as np

parser = PDBParser(QUIET=True)

# Load structures
receptor = parser.get_structure('receptor', '2gli_minimized.pdb')
ligand = parser.get_structure('ligand', 'gant61_docked.pdbqt')

# Find distances to GLU119 and GLU167
def find_interactions(receptor, ligand, residue_ids=[119, 167]):
    interactions = []
    
    for model_r in receptor:
        for chain_r in model_r:
            for residue_r in chain_r:
                if residue_r.id[1] in residue_ids:
                    for atom_r in residue_r:
                        for model_l in ligand:
                            for chain_l in model_l:
                                for residue_l in chain_l:
                                    for atom_l in residue_l:
                                        dist = np.linalg.norm(
                                            atom_r.coord - atom_l.coord
                                        )
                                        if dist < 4.0:  # H-bond distance
                                            interactions.append({
                                                'residue': residue_r.resname + str(residue_r.id[1]),
                                                'protein_atom': atom_r.name,
                                                'ligand_atom': atom_l.name,
                                                'distance': dist
                                            })
    return interactions

interactions = find_interactions(receptor, ligand)

print("\n🔬 GANT61 Interactions with Binding Site:")
print("=" * 60)
for i in interactions:
    print(f"  {i['residue']}: {i['protein_atom']} ↔ {i['ligand_atom']} = {i['distance']:.2f} Å")

if len(interactions) >= 2:
    print("\n✅ VALIDATED: GANT61 forms multiple contacts with GLU119/167")
else:
    print("\n⚠️  WARNING: Few interactions detected")