import numpy as np

# Your grid box parameters
center = [-32.6, -5.7, -0.6]
size = [25, 25, 25]

# Calculate box corners
half_size = [s/2 for s in size]
x_min, x_max = center[0] - half_size[0], center[0] + half_size[0]
y_min, y_max = center[1] - half_size[1], center[1] + half_size[1]
z_min, z_max = center[2] - half_size[2], center[2] + half_size[2]

# Create BILD format file (ChimeraX geometric objects)
with open('gridbox.bild', 'w') as f:
    f.write('.color red\n')
    f.write('.transparency 0.7\n')
    
    # Draw the 12 edges of the box
    # Bottom face (z_min)
    f.write(f'.cylinder {x_min} {y_min} {z_min} {x_max} {y_min} {z_min} 0.3\n')
    f.write(f'.cylinder {x_max} {y_min} {z_min} {x_max} {y_max} {z_min} 0.3\n')
    f.write(f'.cylinder {x_max} {y_max} {z_min} {x_min} {y_max} {z_min} 0.3\n')
    f.write(f'.cylinder {x_min} {y_max} {z_min} {x_min} {y_min} {z_min} 0.3\n')
    
    # Top face (z_max)
    f.write(f'.cylinder {x_min} {y_min} {z_max} {x_max} {y_min} {z_max} 0.3\n')
    f.write(f'.cylinder {x_max} {y_min} {z_max} {x_max} {y_max} {z_max} 0.3\n')
    f.write(f'.cylinder {x_max} {y_max} {z_max} {x_min} {y_max} {z_max} 0.3\n')
    f.write(f'.cylinder {x_min} {y_max} {z_max} {x_min} {y_min} {z_max} 0.3\n')
    
    # Vertical edges
    f.write(f'.cylinder {x_min} {y_min} {z_min} {x_min} {y_min} {z_max} 0.3\n')
    f.write(f'.cylinder {x_max} {y_min} {z_min} {x_max} {y_min} {z_max} 0.3\n')
    f.write(f'.cylinder {x_max} {y_max} {z_min} {x_max} {y_max} {z_max} 0.3\n')
    f.write(f'.cylinder {x_min} {y_max} {z_min} {x_min} {y_max} {z_max} 0.3\n')
    
    # Add center marker
    f.write('.color yellow\n')
    f.write(f'.sphere {center[0]} {center[1]} {center[2]} 1.5\n')

print("✓ Created gridbox.bild")
print(f"  Center: ({center[0]}, {center[1]}, {center[2]})")
print(f"  Size: {size[0]} × {size[1]} × {size[2]} Å")
print(f"  X range: {x_min:.1f} to {x_max:.1f} Å")
print(f"  Y range: {y_min:.1f} to {y_max:.1f} Å")
print(f"  Z range: {z_min:.1f} to {z_max:.1f} Å")