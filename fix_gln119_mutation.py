"""
FIX GLN119 → GLU119 MUTATION
This corrects the structure to have the acidic residue GANT61 needs
"""

from Bio.PDB import PDBParser, PDBIO, Select
import numpy as np
import warnings
warnings.filterwarnings('ignore')

def mutate_gln119_to_glu119(input_pdb='gli_structure/2gli_with_zinc.pdb',
                             output_pdb='gli_structure/2gli_GLU119_corrected.pdb'):
    """
    Mutate GLN119 to GLU119 by changing atom names
    GLN has: N, CA, C, O, CB, CG, CD, OE1, NE2
    GLU has: N, CA, C, O, CB, CG, CD, OE1, OE2
    
    Change: NE2 (amine) → OE2 (carboxyl)
    """
    
    print("=" * 80)
    print("MUTATING GLN119 → GLU119")
    print("=" * 80)
    
    parser = PDBParser(QUIET=True)
    
    try:
        structure = parser.get_structure('GLI', input_pdb)
    except:
        print(f"❌ Cannot find {input_pdb}")
        return False
    
    print(f"\n✓ Loaded: {input_pdb}")
    
    # Find residue 119
    chain = list(structure.get_chains())[0]
    
    residue_119 = None
    for res in chain:
        if res.id[1] == 119 and res.id[0] == ' ':
            residue_119 = res
            break
    
    if residue_119 is None:
        print("❌ Residue 119 not found in structure")
        return False
    
    original_name = residue_119.resname
    print(f"\nFound residue 119: {original_name}")
    
    if original_name == 'GLU':
        print("✓ Already GLU - no mutation needed!")
        return True
    
    if original_name != 'GLN':
        print(f"⚠️  Residue 119 is {original_name}, not GLN")
        print("   Cannot perform GLN→GLU mutation")
        return False
    
    # Perform mutation
    print("\nPerforming mutation GLN→GLU:")
    
    # Change residue name
    residue_119.resname = 'GLU'
    print("  ✓ Changed residue name: GLN → GLU")
    
    # Find and modify NE2 → OE2
    if 'NE2' in residue_119:
        ne2_atom = residue_119['NE2']
        
        # Change atom properties
        ne2_atom.name = 'OE2'
        ne2_atom.fullname = ' OE2'
        ne2_atom.element = 'O'
        
        print("  ✓ Changed atom: NE2 (amine) → OE2 (carboxyl)")
        
        # Adjust geometry slightly (OE2 should be sp2, not sp3)
        # Move slightly to correct bond angle
        if 'CD' in residue_119 and 'OE1' in residue_119:
            cd_coord = residue_119['CD'].coord
            oe1_coord = residue_119['OE1'].coord
            
            # Calculate ideal OE2 position (opposite to OE1)
            cd_to_oe1 = oe1_coord - cd_coord
            oe2_ideal = cd_coord - cd_to_oe1 * 0.95  # Slightly shorter C=O
            
            ne2_atom.coord = oe2_ideal
            print("  ✓ Adjusted OE2 geometry")
    else:
        print("  ⚠️  NE2 atom not found - check structure")
    
    # Save corrected structure
    io = PDBIO()
    io.set_structure(structure)
    io.save(output_pdb)
    
    print(f"\n✓ Saved corrected structure: {output_pdb}")
    
    # Verify mutation
    print("\n" + "=" * 80)
    print("VERIFICATION:")
    print("=" * 80)
    
    # Re-load and check
    verify_structure = parser.get_structure('verify', output_pdb)
    verify_chain = list(verify_structure.get_chains())[0]
    
    for res in verify_chain:
        if res.id[1] == 119 and res.id[0] == ' ':
            print(f"\nResidue 119 after mutation: {res.resname}")
            
            if res.resname == 'GLU':
                print("✅ MUTATION SUCCESSFUL")
                
                # Check atoms
                print("\nAtoms in GLU119:")
                for atom in res:
                    print(f"  {atom.name:4s} ({atom.element})")
                
                # Check for OE1 and OE2
                has_oe1 = 'OE1' in res
                has_oe2 = 'OE2' in res
                
                if has_oe1 and has_oe2:
                    print("\n✓ Both carboxyl oxygens present (OE1, OE2)")
                    print("✓ GLU119 is now ACIDIC (can bind GANT61)")
                    return True
                else:
                    print("\n⚠️  Missing carboxyl oxygen atoms")
                    return False
            else:
                print("❌ MUTATION FAILED")
                return False
    
    print("❌ Could not verify residue 119")
    return False

def prepare_corrected_receptor():
    """
    After mutation, prepare the corrected receptor for docking
    """
    print("\n" + "=" * 80)
    print("PREPARING CORRECTED RECEPTOR FOR DOCKING")
    print("=" * 80)
    
    print("""
Next steps:
1. Add hydrogens to corrected structure
2. Energy minimize
3. Convert to PDBQT
4. Re-dock GANT61 to validate

Run this command:
    python prepare_gli.py  (with corrected input file)
    
Or manually:
    obabel 2gli_GLU119_corrected.pdb -O corrected.pdb -h
    prepare_receptor4.py -r corrected.pdb -o corrected.pdbqt
    """)

def main():
    print("\n" + "=" * 80)
    print("GLI1 GLN119→GLU119 CORRECTION")
    print("=" * 80 + "\n")
    
    # Try both possible file locations
    input_files = [
        'gli_structure/2gli_with_zinc.pdb',
        '2gli_with_zinc.pdb',
        'gli_structure/2gli_min.pdb',
        '2gli_min.pdb',
        'gli_structure/2gli_clean.pdb',
        '2gli_clean.pdb'
    ]
    
    success = False
    for input_file in input_files:
        try:
            if mutate_gln119_to_glu119(input_file):
                success = True
                break
        except:
            continue
    
    if not success:
        print("\n❌ Could not perform mutation")
        print("   Make sure GLI structure files are in current directory")
        return
    
    prepare_corrected_receptor()
    
    print("\n" + "=" * 80)
    print("NEXT: RE-DOCK GANT61 WITH CORRECTED STRUCTURE")
    print("=" * 80)
    print("""
After preparing the corrected receptor:
1. Re-run GANT61 validation docking
2. Check if score improves (should be < -7.16)
3. Verify GANT61 now contacts GLU119 (not GLN119)
4. If successful, proceed with 483K screening
    """)

if __name__ == "__main__":
    main()