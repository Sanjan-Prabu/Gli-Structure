import subprocess

# Install SMINA
# conda install -c conda-forge smina

# Dock with SMINA (different scoring)
cmd = """smina --receptor 2gli_receptor.pdbqt \
              --ligand gant61.pdbqt \
              --center_x -32.6 --center_y -5.7 --center_z -0.6 \
              --size_x 25 --size_y 25 --size_z 25 \
              --exhaustiveness 32 \
              --out gant61_smina.pdbqt"""

subprocess.run(cmd, shell=True)