python
import os

# Get PyMOL's current working directory (where this script is loaded from)
script_dir = os.getcwd()

# Go up one directory to get the project base (drug-target-pocket-funnel/)
base_dir = os.path.abspath(os.path.join(script_dir, ".."))

# Construct paths to input PDBs
af_path = os.path.join(base_dir, "outputs", "ftmap_af_ppi_mode", "fftmap.170801.pdb")
f42_path = os.path.join(base_dir, "outputs", "ftmap_1f42_ppi_mode", "fftmap.145217.pdb")

# Load and rename AF objects
cmd.load(af_path)
af_objects = cmd.get_names('objects')
for name in af_objects:
    cmd.set_name(name, name + '_af')

# Load and rename 1f42 objects
cmd.load(f42_path)
f42_objects = cmd.get_names('objects')
for name in f42_objects:
    if not name.endswith('_af'):
        cmd.set_name(name, name + '_1f42')

# Align 1f42 protein to AF protein
cmd.align("protein_1f42", "protein_af")

# Apply the same transformation to all other 1f42 objects
matrix = cmd.get_object_matrix("protein_1f42")
for obj in cmd.get_names('objects'):
    if obj.endswith("_1f42") and obj != "protein_1f42":
        cmd.transform_object(obj, matrix)

# Zoom to the scene
cmd.zoom("all")

python end
