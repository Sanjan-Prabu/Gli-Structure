"""
Download AlphaFold GLI1 and verify GLU119/GLU167 exist
Then prepare for targeted docking
"""

import urllib.request
import os
import numpy as np
from Bio.PDB import PDBParser, PDBIO, Select
import warnings
warnings.filterwarnings('ignore')

def download_alphafold_gli1():
    """Download full-length AlphaFold GLI1 structure"""
    
    print("=" * 80)
    print("STEP 1: DOWNLOADING ALPHAFOLD GLI1")
    print("=" * 80)
    
    os.makedirs('alphafold_gli', exist_ok=True)
    os.chdir('alphafold_gli')
    
    # AlphaFold GLI1 (UniProt P08151)
    print(f"\nDownloading GLI1 from AlphaFold DB...")
    print(f"UniProt: P08151 (GLI1)")
    
    # Try multiple URL formats
    urls_to_try = [
        "https://alphafold.ebi.ac.uk/files/AF-P08151-F1-model_v4.pdb",
        "https://alphafold.ebi.ac.uk/files/AF-P08151-F1-model_v3.pdb",
        "https://alphafold.ebi.ac.uk/files/AF-P08151-F1-model_v2.pdb",
    ]
    
    for url in urls_to_try:
        try:
            print(f"\nTrying: {url}")
            urllib.request.urlretrieve(url, 'GLI1_alphafold_full.pdb')
            print("✓ Downloaded: GLI1_alphafold_full.pdb")
            return 'GLI1_alphafold_full.pdb'
        except Exception as e:
            print(f"  Failed: {e}")
            continue
    
    # If automatic download fails, provide manual instructions
    print("\n" + "=" * 80)
    print("❌ AUTOMATIC DOWNLOAD FAILED")
    print("=" * 80)
    print("\n📋 MANUAL DOWNLOAD INSTRUCTIONS:")
    print("   1. Go to: https://alphafold.ebi.ac.uk/entry/P08151")
    print("   2. Click 'Download structure' button")
    print("   3. Select 'PDB format'")
    print("   4. Save as: GLI1_alphafold_full.pdb")
    print("   5. Move file to: alphafold_gli/ directory")
    print("   6. Re-run this script")
    print("\n" + "=" * 80)
    
    return None

def verify_glu119_glu167(pdb_file):
    """
    Check if GLU119 and GLU167 exist in AlphaFold structure
    """
    
    print("\n" + "=" * 80)
    print("STEP 2: VERIFYING GLU119 AND GLU167")
    print("=" * 80)
    
    parser = PDBParser(QUIET=True)
    structure = parser.get_structure('GLI1', pdb_file)
    
    chain = list(structure.get_chains())[0]
    
    # Find residues 119 and 167
    residue_119 = None
    residue_167 = None
    
    for res in chain:
        if res.id[0] == ' ':  # Standard residue
            if res.id[1] == 119:
                residue_119 = res
            elif res.id[1] == 167:
                residue_167 = res
    
    print("\n🔍 RESIDUE VERIFICATION:")
    
    # Check residue 119
    if residue_119:
        res_name = residue_119.resname
        print(f"\n  Residue 119: {res_name}")
        
        if res_name == 'GLU':
            print("  ✅ CORRECT - GLU119 (Glutamic acid - ACIDIC)")
            glu119_correct = True
        elif res_name == 'GLN':
            print("  ❌ WRONG - GLN119 (Glutamine - NEUTRAL)")
            print("     Literature requires GLU (acidic)")
            glu119_correct = False
        else:
            print(f"  ❌ UNEXPECTED - {res_name}119")
            glu119_correct = False
        
        # Show atoms
        if 'CA' in residue_119:
            coord = residue_119['CA'].coord
            print(f"  Location: ({coord[0]:.2f}, {coord[1]:.2f}, {coord[2]:.2f})")
    else:
        print("\n  ❌ Residue 119 NOT FOUND")
        glu119_correct = False
    
    # Check residue 167
    if residue_167:
        res_name = residue_167.resname
        print(f"\n  Residue 167: {res_name}")
        
        if res_name == 'GLU':
            print("  ✅ CORRECT - GLU167 (Glutamic acid - ACIDIC)")
            glu167_correct = True
        else:
            print(f"  ❌ WRONG - Expected GLU, found {res_name}")
            glu167_correct = False
        
        # Show atoms
        if 'CA' in residue_167:
            coord = residue_167['CA'].coord
            print(f"  Location: ({coord[0]:.2f}, {coord[1]:.2f}, {coord[2]:.2f})")
    else:
        print("\n  ❌ Residue 167 NOT FOUND")
        glu167_correct = False
    
    # Calculate binding site center if both present
    if residue_119 and residue_167:
        if 'CA' in residue_119 and 'CA' in residue_167:
            coord_119 = residue_119['CA'].coord
            coord_167 = residue_167['CA'].coord
            
            center = (coord_119 + coord_167) / 2
            
            print("\n" + "=" * 80)
            print("BINDING SITE CENTER (GLU119/GLU167):")
            print("=" * 80)
            print(f"\n  Center: ({center[0]:.2f}, {center[1]:.2f}, {center[2]:.2f})")
            
            distance = np.linalg.norm(coord_119 - coord_167)
            print(f"  Distance between residues: {distance:.2f} Å")
            
            if distance > 50:
                print("  ⚠️  Residues very far apart - check structure")
            elif distance < 5:
                print("  ⚠️  Residues very close - check numbering")
            else:
                print("  ✓ Reasonable separation")
            
            return glu119_correct and glu167_correct, center
    
    return False, None

def check_confidence_scores(pdb_file):
    """
    Check AlphaFold pLDDT confidence scores for residues 119 and 167
    """
    
    print("\n" + "=" * 80)
    print("STEP 3: ALPHAFOLD CONFIDENCE SCORES")
    print("=" * 80)
    
    parser = PDBParser(QUIET=True)
    structure = parser.get_structure('GLI1', pdb_file)
    
    chain = list(structure.get_chains())[0]
    
    print("\npLDDT confidence scores (0-100):")
    print("  > 90: Very high confidence")
    print("  70-90: High confidence")
    print("  50-70: Low confidence")
    print("  < 50: Very low confidence")
    
    for res in chain:
        if res.id[0] == ' ':
            res_num = res.id[1]
            
            # Check residues in zinc finger region (100-200)
            if res_num in [119, 167] or (100 <= res_num <= 200):
                # pLDDT is stored in B-factor column
                if 'CA' in res:
                    plddt = res['CA'].bfactor
                    
                    if res_num in [119, 167]:
                        confidence = "High" if plddt > 70 else "Low"
                        marker = "✓" if plddt > 70 else "⚠️"
                        print(f"\n  {marker} {res.resname}{res_num}: pLDDT = {plddt:.1f} ({confidence})")
    
    print("\n" + "=" * 80)

def extract_zinc_finger_region(pdb_file, output_file='GLI1_ZF_domain.pdb'):
    """
    Extract zinc finger domain (residues 230-380 or region around 119-167)
    """
    
    print("\n" + "=" * 80)
    print("STEP 4: EXTRACTING RELEVANT DOMAIN")
    print("=" * 80)
    
    parser = PDBParser(QUIET=True)
    structure = parser.get_structure('GLI1', pdb_file)
    
    # Determine region to extract
    # Option 1: If 119/167 are in zinc finger domain (230-380)
    # Option 2: Extract broader region around 119 and 167
    
    chain = list(structure.get_chains())[0]
    
    # Find actual residue range
    residue_numbers = []
    for res in chain:
        if res.id[0] == ' ':
            residue_numbers.append(res.id[1])
    
    print(f"\nTotal residues: {len(residue_numbers)}")
    print(f"Range: {min(residue_numbers)} to {max(residue_numbers)}")
    
    # Strategy: Extract region that includes both 119 and 167
    # Plus 50 residues buffer on each side
    
    if 119 in residue_numbers and 167 in residue_numbers:
        start = max(1, 119 - 50)
        end = 167 + 50
        
        print(f"\nExtracting residues {start} to {end}")
        print("  (GLU119/GLU167 region + buffer)")
        
        class DomainSelect(Select):
            def accept_residue(self, residue):
                if residue.id[0] != ' ':
                    return False
                return start <= residue.id[1] <= end
        
        io = PDBIO()
        io.set_structure(structure)
        io.save(output_file, DomainSelect())
        
        print(f"\n✓ Saved: {output_file}")
        return output_file
    
    else:
        print("\n⚠️  Could not find residues 119 and 167")
        print("   Extracting full structure instead")
        
        io = PDBIO()
        io.set_structure(structure)
        io.save(output_file)
        
        return output_file

def generate_chimerax_script(center, pdb_file='GLI1_ZF_domain.pdb'):
    """
    Generate ChimeraX script to visualize GLU119/GLU167 binding site
    """
    
    print("\n" + "=" * 80)
    print("STEP 5: GENERATING CHIMERAX VISUALIZATION")
    print("=" * 80)
    
    script = f"""# ChimeraX Visualization Script for GLI1 GLU119/GLU167
# Save this as: visualize_binding_site.cxc

# Open structure
open {pdb_file}

# Color by pLDDT confidence (AlphaFold)
color bfactor palette alphafold

# Show cartoon representation
cartoon
hide atoms

# Highlight GLU119 and GLU167
select :119,167
color sel red
show sel atoms
label sel residues
style sel sphere

# Show binding site center as sphere
shape sphere radius 3 center {center[0]:.2f},{center[1]:.2f},{center[2]:.2f} color blue

# Show grid box (25x25x25 Å)
shape box size 25,25,25 center {center[0]:.2f},{center[1]:.2f},{center[2]:.2f} color yellow

# Focus view
view sel

# Show nearby residues (within 10 Å)
select :119,167 zr<10
show sel atoms
color sel byhetero

# Label acidic residues near binding site
select :GLU,ASP zr<15 from :119,167
label sel residues
color sel orange

# Final view
~select
view

# Save session
save GLI1_binding_site.cxs

echo "✓ Visualization complete!"
echo "  Red spheres: GLU119 and GLU167"
echo "  Blue sphere: Binding site center"
echo "  Yellow box: 25x25x25 Å docking grid"
echo "  Orange: Nearby acidic residues"
"""
    
    with open('visualize_binding_site.cxc', 'w') as f:
        f.write(script)
    
    print("\n✓ Created: visualize_binding_site.cxc")
    print("\nTo visualize in ChimeraX:")
    print("  1. Open ChimeraX")
    print("  2. File → Open → visualize_binding_site.cxc")
    print("  3. Or run: chimerax visualize_binding_site.cxc")

def prepare_for_targeted_docking(center, output_config='targeted_vina_config.txt'):
    """
    Generate Vina config file for targeted docking at GLU119/GLU167
    """
    
    print("\n" + "=" * 80)
    print("STEP 6: TARGETED DOCKING CONFIGURATION")
    print("=" * 80)
    
    config = f"""# AutoDock Vina Configuration
# Targeted docking at GLU119/GLU167 site

# Binding site center (midpoint of GLU119 and GLU167)
center_x = {center[0]:.2f}
center_y = {center[1]:.2f}
center_z = {center[2]:.2f}

# Grid box size (covers GLU119-GLU167 region)
size_x = 25
size_y = 25
size_z = 25

# Docking parameters
exhaustiveness = 32  # Higher for validation
num_modes = 20       # More poses to analyze
energy_range = 5

# Receptor and ligand (update these paths)
receptor = GLI1_prepared.pdbqt
ligand = gant61.pdbqt

# Output
out = gant61_GLU119_167_docked.pdbqt
log = gant61_GLU119_167_docking.log
"""
    
    with open(output_config, 'w') as f:
        f.write(config)
    
    print(f"\n✓ Created: {output_config}")
    print(f"\nBinding site centered at: ({center[0]:.2f}, {center[1]:.2f}, {center[2]:.2f})")
    print("  → Midpoint between GLU119 and GLU167")

def main():
    """
    Complete pipeline to verify and prepare AlphaFold GLI1 structure
    """
    
    print("\n" + "=" * 80)
    print("ALPHAFOLD GLI1 VERIFICATION & TARGETED DOCKING SETUP")
    print("=" * 80)
    print("\nObjective: Verify GLU119/GLU167 exist and prepare for targeted docking")
    print("=" * 80 + "\n")
    
    # Step 1: Check if file already exists
    possible_files = [
        'alphafold_gli/AF-P08151-F1-v6.pdb',
        'alphafold_gli/GLI1_alphafold_full.pdb',
        'AF-P08151-F1-v6.pdb',
        'GLI1_alphafold_full.pdb',
        'alphafold_gli/GLI1_alphafold_full.pdb',
    ]
    
    pdb_file = None
    for filepath in possible_files:
        if os.path.exists(filepath):
            pdb_file = filepath
            print(f"✓ Found existing file: {filepath}")
            break
    
    if not pdb_file:
        # Try to download
        pdb_file = download_alphafold_gli1()
        if not pdb_file:
            print("\n❌ No AlphaFold structure found")
            print("   Please download manually and place in alphafold_gli/ folder")
            return
    
    # Step 2: Verify residues
    correct, center = verify_glu119_glu167(pdb_file)
    
    if not correct:
        print("\n" + "=" * 80)
        print("❌ CRITICAL ISSUE: GLU119 and/or GLU167 not correct")
        print("=" * 80)
        return
    
    # Step 3: Check confidence
    check_confidence_scores(pdb_file)
    
    # Step 4: Extract domain
    domain_file = extract_zinc_finger_region(pdb_file)
    
    # Step 5: Generate visualization
    generate_chimerax_script(center, domain_file)
    
    # Step 6: Prepare docking config
    prepare_for_targeted_docking(center)
    
    # Final summary
    print("\n" + "=" * 80)
    print("✅ SETUP COMPLETE")
    print("=" * 80)
    
    print(f"""
Next Steps:

1. VISUALIZE IN CHIMERAX:
   chimerax visualize_binding_site.cxc
   
   → Verify GLU119 and GLU167 are visible
   → Confirm they are in DNA-binding groove
   → Check pLDDT confidence (color coded)

2. PREPARE RECEPTOR:
   - Add hydrogens: obabel {domain_file} -O prepared.pdb -h
   - Add zinc ions (if needed)
   - Minimize structure
   - Convert to PDBQT

3. DOCK GANT61 (Targeted):
   - Use binding center: ({center[0]:.2f}, {center[1]:.2f}, {center[2]:.2f})
   - Grid box: 25×25×25 Å
   - Exhaustiveness: 32 (high accuracy)

4. VERIFY BINDING:
   - Check if GANT61 contacts GLU119 and GLU167
   - Score should be -7 to -9 kcal/mol
   - If successful → proceed with 483K screening

IF GANT61 DOES NOT BIND TO GLU119/GLU167:
   → Literature may be incorrect
   → Use alternative binding site found earlier
   → Validate with experimental data
    """)
    
    print("=" * 80)
    print(f"✓ AlphaFold structure: {pdb_file}")
    print(f"✓ Extracted domain: {domain_file}")
    print(f"✓ ChimeraX script: visualize_binding_site.cxc")
    print(f"✓ Vina config: targeted_vina_config.txt")
    print("=" * 80 + "\n")

if __name__ == "__main__":
    main()