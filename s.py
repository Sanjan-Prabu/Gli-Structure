import os  # <-- FIX 1: Must import os to use os.path.exists and os.rename
from Bio.PDB import PDBList, PDBParser

print("\nDownloading fresh 2GLI from RCSB...")
pdbl = PDBList()
# BioPython downloads PDB 2GLI as 'pdb2gli.ent' by default.
pdbl.retrieve_pdb_file('2GLI', pdir='.', file_format='pdb')

# Define the expected downloaded filename and the desired new filename
downloaded_filename = 'pdb2gli.ent'
new_filename = 'fresh_2gli.pdb'

# --- FIX 2: Ensure we check for the expected file name ---
if os.path.exists(downloaded_filename):
    os.rename(downloaded_filename, new_filename)
elif os.path.exists('./2gli.pdb'):
    # Check for another common BioPython output format
    os.rename('./2gli.pdb', new_filename)
elif os.path.exists('2gli.pdb'):
    # Check for another common download output
    os.rename('2gli.pdb', new_filename)
else:
    print(f"❌ Error: Could not find the downloaded PDB file ({downloaded_filename} or 2gli.pdb).")
    exit() # Stop the script if the file isn't found

parser = PDBParser(QUIET=True)

# --- FIX 3: Ensure the new filename is used ---
try:
    structure = parser.get_structure('fresh', new_filename)

    print("\n" + "="*60)
    print("FRESH 2GLI FROM RCSB")
    print("="*60)

    # Check the residues at positions 119 and 167
    for chain in structure[0]:
        if chain.id == 'A':
            for res in chain:
                # res.id[1] is the residue number
                # res.id[0] == ' ' ensures it's a standard AA residue (not HETATM)
                if res.id[1] in [119, 167] and res.id[0] == ' ':
                    print(f"Residue {res.id[1]}: {res.resname}")

except FileNotFoundError:
    print(f"❌ Error: The file {new_filename} was not found after renaming.")
except Exception as e:
    print(f"An error occurred during parsing: {e}")