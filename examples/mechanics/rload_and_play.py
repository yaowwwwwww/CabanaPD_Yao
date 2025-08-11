from paraview.simple import *
import glob
import os

# Set your output directory
output_dir = "/home/wuwen/program/CabanaPD_Yao/output"

# === Load only the last file ===
# Get a sorted list of all files matching the pattern
files = sorted(glob.glob(os.path.join(output_dir, "particles_*.silo")))
if not files:
    raise RuntimeError(f"❌ No files found: {output_dir}/particles_*.silo")

# Select only the last file from the sorted list
last_file = files[-1]
print(f"Loading only the last file: {last_file}")

# Open the last file directly
reader = OpenDataFile(last_file)

# === Get the active view ===
renderView1 = GetActiveViewOrCreate('RenderView')

# === Set active point arrays for coloring and display ===
# We only need 'displacements' for this task, so we set that specifically.
if hasattr(reader, 'PointArrayStatus'):
    reader.PointArrayStatus = ["rank_0/displacements"]

# === Show the data ===
# The clip filter is removed, so we show the reader's output directly.
display = Show(reader, renderView1)

# === Set display properties to Point Gaussian ===
display.Representation = 'Point Gaussian'
display.GaussianRadius = 0.001  # Adjust as needed
# Color by displacement magnitude
ColorBy(display, ('POINTS', 'rank_0/displacements', 'Magnitude'))
# Ensure the color bar is visible
displacementsLUT = GetColorTransferFunction('rank_0/displacements')
displacementsLUT.ApplyPreset('Viridis (matplotlib)', True) # Example: use a different colormap

# === Render the scene ===
Render()

# === Set camera properties ===
# Adjust these based on your model's actual extent and desired view
renderView1.CameraPosition = [0, -0.05, 0.02]
renderView1.CameraFocalPoint = [0, 0, 0.005]
renderView1.CameraViewUp = [0, 0, 1]
ResetCamera() # Use ResetCamera() to automatically fit the model in the view

# === Play animation is removed, as we only loaded one timestep ===