"""
URGENT: Determine WHERE GANT61 actually docked
Check if it's at GLN119/GLU167 or an allosteric pocket
"""

import numpy as np
from Bio.PDB import PDBParser, NeighborSearch, Selection
import warnings
warnings.filterwarnings('ignore')

def analyze_gant61_binding_site():
    """
    Analyze the docked GANT61 pose to find:
    1. Which residues it actually contacts
    2. If it's near GLN119 (wrong) or an allosteric pocket
    3. Distance to zinc ions (ZF2-3 interface)
    """
    
    print("=" * 80)
    print("URGENT: GANT61 BINDING SITE ANALYSIS")
    print("=" * 80)
    
    # Load protein structure
    parser = PDBParser(QUIET=True)
    
    try:
        protein = parser.get_structure('GLI', 'gli_structure/2gli_receptor.pdbqt')
    except:
        try:
            protein = parser.get_structure('GLI', '2gli_receptor.pdbqt')
        except:
            print("❌ Cannot find 2gli_receptor.pdbqt")
            return
    
    # Load GANT61 docked pose
    try:
        ligand = parser.get_structure('GANT61', 'gli_structure/gant61_docked.pdbqt')
    except:
        try:
            ligand = parser.get_structure('GANT61', 'gant61_docked.pdbqt')
        except:
            print("❌ Cannot find gant61_docked.pdbqt")
            print("   Run docking first to generate this file")
            return
    
    print("\n✓ Loaded structures")
    
    # Get all protein atoms
    protein_atoms = list(protein.get_atoms())
    ligand_atoms = list(ligand.get_atoms())
    
    # Calculate ligand center of mass
    ligand_coords = np.array([atom.coord for atom in ligand_atoms])
    ligand_center = np.mean(ligand_coords, axis=0)
    
    print(f"\nGANT61 binding center: ({ligand_center[0]:.2f}, {ligand_center[1]:.2f}, {ligand_center[2]:.2f})")
    
    # Find all residues within 5 Å of GANT61
    print("\n" + "=" * 80)
    print("RESIDUES IN CONTACT WITH GANT61 (within 5 Å):")
    print("=" * 80)
    
    ns = NeighborSearch(protein_atoms)
    
    contacting_residues = set()
    for ligand_atom in ligand_atoms:
        nearby = ns.search(ligand_atom.coord, 5.0, level='R')
        for res in nearby:
            if res.id[0] == ' ':  # Standard residue
                contacting_residues.add(res)
    
    # Organize by residue number
    contact_list = []
    for res in contacting_residues:
        res_num = res.id[1]
        res_name = res.resname
        
        # Calculate closest distance
        min_dist = 999
        for res_atom in res:
            for lig_atom in ligand_atoms:
                dist = np.linalg.norm(res_atom.coord - lig_atom.coord)
                if dist < min_dist:
                    min_dist = dist
        
        contact_list.append((res_num, res_name, min_dist))
    
    # Sort by residue number
    contact_list.sort(key=lambda x: x[0])
    
    print(f"\nFound {len(contact_list)} residues in contact:\n")
    
    # Critical residue flags
    found_gln119 = False
    found_glu167 = False
    found_glu_in_contact = []
    found_arg_in_contact = []
    
    for res_num, res_name, dist in contact_list:
        # Mark critical residues
        marker = ""
        if res_num == 119:
            marker = " ⚠️  GLN119 (MUTANT - should be GLU)"
            found_gln119 = True
        elif res_num == 167:
            marker = " ✓ GLU167 (CORRECT - acidic)"
            found_glu167 = True
        elif res_name == 'GLU':
            marker = " ✓ Glutamic acid (acidic)"
            found_glu_in_contact.append(res_num)
        elif res_name == 'ARG' or res_name == 'LYS':
            marker = " + Positive residue"
            found_arg_in_contact.append(res_num)
        
        print(f"  {res_name}{res_num:3d}: {dist:4.2f} Å{marker}")
    
    # Find zinc ions and their distances
    print("\n" + "=" * 80)
    print("DISTANCE TO ZINC IONS (ZF2-3 interface check):")
    print("=" * 80)
    
    zinc_distances = []
    chain = list(protein.get_chains())[0]
    
    for res in chain:
        if res.resname == 'ZN':
            for atom in res:
                if 'ZN' in atom.name or atom.element == 'ZN':
                    zn_coord = atom.coord
                    dist = np.linalg.norm(ligand_center - zn_coord)
                    zinc_distances.append((res.id[1], dist))
                    print(f"  ZN {res.id[1]}: {dist:.2f} Å from GANT61")
    
    # Analysis
    print("\n" + "=" * 80)
    print("CRITICAL ANALYSIS:")
    print("=" * 80)
    
    print("\n🔍 BINDING SITE IDENTIFICATION:")
    
    if found_gln119:
        print("\n  ❌ PROBLEM: GANT61 is contacting GLN119 (neutral)")
        print("     Literature says GANT61 needs GLU119 (acidic)")
        print("     This is likely a FALSE POSITIVE docking")
    
    if found_glu167:
        print("\n  ✓ GLU167 (acidic) is in contact")
        if not found_gln119:
            print("    → GANT61 may be binding primarily to GLU167")
    
    if len(found_glu_in_contact) > 0:
        print(f"\n  ✓ Other GLU residues in contact: {found_glu_in_contact}")
        print("    → These could be alternative binding sites")
    
    # Check if it's in allosteric pocket
    if zinc_distances:
        closest_zn_dist = min([d[1] for d in zinc_distances])
        if closest_zn_dist < 15:
            print(f"\n  🎯 GANT61 is {closest_zn_dist:.1f} Å from nearest zinc")
            print("     → Binding NEAR zinc fingers (possible allosteric pocket)")
            
            # Check if between ZF2-3
            if len(zinc_distances) >= 3:
                zn2_dist = zinc_distances[1][1]
                zn3_dist = zinc_distances[2][1]
                if zn2_dist < 20 and zn3_dist < 20:
                    print("     → Located between ZF2-3 interface ✓")
        else:
            print(f"\n  ⚠️  GANT61 is {closest_zn_dist:.1f} Å from zinc (FAR)")
            print("     → Not in zinc finger region")
    
    # VERDICT
    print("\n" + "=" * 80)
    print("VERDICT & NEXT STEPS:")
    print("=" * 80)
    
    if found_gln119 and not found_glu167:
        print("\n❌ CRITICAL ISSUE: Docking to GLN119 (wrong residue)")
        print("\n📋 SOLUTION OPTIONS:")
        print("   1. MUTATE GLN119→GLU119 in structure (RECOMMENDED)")
        print("   2. Search for alternative binding site with GLU residues")
        print("   3. Use different GLI structure (AlphaFold full-length)")
        
    elif found_glu167 or len(found_glu_in_contact) > 0:
        print("\n✓ POTENTIALLY VALID: Binding to acidic residues")
        if found_glu167:
            print("  → GLU167 present (correct residue)")
        if found_glu_in_contact:
            print(f"  → Alternative GLU contacts: {found_glu_in_contact}")
        print("\n📋 RECOMMEND: Verify this is the real GANT61 binding mode")
        print("   - Check literature for specific residues")
        print("   - Consider mutagenesis validation")
        
    else:
        print("\n⚠️  UNCLEAR: No acidic residues in strong contact")
        print("\n📋 URGENT ACTION NEEDED:")
        print("   1. Check if binding site makes chemical sense")
        print("   2. Re-dock with corrected structure (GLU119)")
        print("   3. Use experimental data to validate")
    
    print("\n" + "=" * 80)
    
    return contact_list, found_gln119, found_glu167, found_glu_in_contact

def suggest_solutions():
    """Provide concrete next steps"""
    
    print("\n" + "=" * 80)
    print("SOLUTION 1: MUTATE GLN119 → GLU119 (FASTEST)")
    print("=" * 80)
    print("""
This is the most direct fix:
1. Edit the PDB file and change GLN119 to GLU119
2. Re-minimize the structure
3. Re-dock GANT61 to validate

Script: I'll create this for you
    """)
    
    print("\n" + "=" * 80)
    print("SOLUTION 2: SEARCH FOR ALLOSTERIC POCKET")
    print("=" * 80)
    print("""
Use pocket detection tools:
1. Run FPocket on structure to find all pockets
2. Filter for pockets with GLU/ASP residues
3. Check if any are between ZF2-3
4. Re-dock to alternative site

Tool: FPocket or SiteMap
    """)
    
    print("\n" + "=" * 80)
    print("SOLUTION 3: USE ALPHAFOLD FULL-LENGTH GLI1")
    print("=" * 80)
    print("""
Get the complete GLI1 structure:
1. Download AlphaFold GLI1 (UniProt P08151)
2. Extract zinc finger region with CORRECT sequence
3. This should have GLU119 (not GLN)
4. Prepare and re-dock

Risk: Lower confidence regions (pLDDT < 70)
    """)

if __name__ == "__main__":
    contacts, has_gln119, has_glu167, other_glu = analyze_gant61_binding_site()
    
    print("\n")
    suggest_solutions()
    
    print("\n" + "=" * 80)
    print("⚠️  DO NOT PROCEED WITH SCREENING UNTIL THIS IS RESOLVED")
    print("=" * 80)