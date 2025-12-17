import py3Dmol

# Load your docked GANT61
with open('gant61_docked.pdbqt', 'r') as f:
    gant61_data = f.read()

# Load GLI receptor
with open('2gli_receptor.pdbqt', 'r') as f:
    receptor_data = f.read()

# Create 3D viewer
view = py3Dmol.view(width=1800, height=1000)

# Add receptor (protein)
view.addModel(receptor_data, 'pdbqt')
view.setStyle({'model': 0}, {'cartoon': {'color': 'spectrum'}})

# Add GANT61 (ligand)
view.addModel(gant61_data, 'pdbqt')
view.setStyle({'model': 1}, {'stick': {'colorscheme': 'greenCarbon'}})

# Highlight binding site residues (GLU119, GLU167)
view.addStyle({'model': 0, 'resi': [119, 167]}, 
              {'stick': {'colorscheme': 'redCarbon', 'radius': 0.3}})

# Add labels
view.addLabel('GLU119 Binding Site', {'resi': 119}, 
              {'backgroundColor': 'white', 'fontSize': 14})
view.addLabel('GLU167 Binding Site', {'resi': 167}, 
              {'backgroundColor': 'white', 'fontSize': 14})

# Center on binding site
view.zoomTo({'resi': [119, 167]})

# Save
view.write_html('gli_gant61_complex.html')
print("✓ Open gli_gant61_complex.html in browser to see 3D structure")