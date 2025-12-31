import os

# Your grid box parameters
center = [-32.6, -5.7, -0.6]
size = [25, 25, 25]

# Calculate box corners
half_size = [s/2 for s in size]
x_min, x_max = center[0] - half_size[0], center[0] + half_size[0]
y_min, y_max = center[1] - half_size[1], center[1] + half_size[1]
z_min, z_max = center[2] - half_size[2], center[2] + half_size[2]

# Create BILD file
bild_content = f""".color red
.transparency 0.7

# Bottom face
.cylinder {x_min} {y_min} {z_min} {x_max} {y_min} {z_min} 0.3
.cylinder {x_max} {y_min} {z_min} {x_max} {y_max} {z_min} 0.3
.cylinder {x_max} {y_max} {z_min} {x_min} {y_max} {z_min} 0.3
.cylinder {x_min} {y_max} {z_min} {x_min} {y_min} {z_min} 0.3

# Top face
.cylinder {x_min} {y_min} {z_max} {x_max} {y_min} {z_max} 0.3
.cylinder {x_max} {y_min} {z_max} {x_max} {y_max} {z_max} 0.3
.cylinder {x_max} {y_max} {z_max} {x_min} {y_max} {z_max} 0.3
.cylinder {x_min} {y_max} {z_max} {x_min} {y_min} {z_max} 0.3

# Vertical edges
.cylinder {x_min} {y_min} {z_min} {x_min} {y_min} {z_max} 0.3
.cylinder {x_max} {y_min} {z_min} {x_max} {y_min} {z_max} 0.3
.cylinder {x_max} {y_max} {z_min} {x_max} {y_max} {z_max} 0.3
.cylinder {x_min} {y_max} {z_min} {x_min} {y_max} {z_max} 0.3

# Center marker
.color yellow
.sphere {center[0]} {center[1]} {center[2]} 1.5
"""

# Get current directory
current_dir = os.getcwd()
bild_file = os.path.join(current_dir, 'gridbox.bild')

# Write file
with open(bild_file, 'w') as f:
    f.write(bild_content)

print("✓ Created gridbox.bild")
print(f"📂 Location: {bild_file}")
print(f"\n📋 Copy this command for ChimeraX:")
print(f"   open {bild_file}")

# Create ChimeraX command file
chimerax_script = f"""# ChimeraX Visualization Script
# Copy and paste these commands into ChimeraX

# Open files
open {os.path.join(current_dir, '2gli_receptor.pdbqt')}
open {bild_file}
open {os.path.join(current_dir, 'gant61_docked.pdbqt')}

# Style protein
cartoon
color bychain

# Highlight E119 and E167
select /A:119,167
color sel red
show sel atoms
style sel sphere
label sel text "{{0.name}}{{0.number}}" height 2

# Style GANT61
color ligand yellow
style ligand stick

# Measure distances
distance /A:119@CA marker #2:1
distance /A:167@CA marker #2:1

# Center view
view all
"""

script_file = os.path.join(current_dir, 'chimerax_commands.cxc')
with open(script_file, 'w') as f:
    f.write(chimerax_script)

print(f"\n✓ Created ChimeraX script: {script_file}")
print(f"\n🚀 To use:")
print(f"   1. Open ChimeraX")
print(f"   2. Run: open {script_file}")