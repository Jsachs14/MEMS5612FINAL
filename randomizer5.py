#!/usr/bin/env python3
import os
import sys
import random

###############################################################################
# Automated Pipeline:
#  1) Use Atomsk to:
#      - Create a BCC Fe cell (lattice = 3.3 Å).
#      - Duplicate it 60×60×60 times (~20 nm per side, i.e. 198 Å per side).
#      - Generate a 30-grain polycrystal (poly.cfg).
#
#  2) Parse poly.cfg (where all atoms are initially "Fe").
#
#  3) Randomly assign ~20% each to V, Nb, Ta, Ti, Zr.
#
#  4) Check if atom coordinates are given in fractional units by computing
#     the range. If the range is small (<2), we scale all coordinates by 198.
#
#  5) Write a LAMMPS data file (final_lammps.data) with 5 atom types:
#         type 1: V
#         type 2: Nb
#         type 3: Ta
#         type 4: Ti
#         type 5: Zr
#
#  NOTE:
#    - Ensure Atomsk is installed and on your PATH.
#    - Adjust the Atomsk commands if your version requires different syntax.
###############################################################################

# -------------
# ATOMSK PART
# -------------

print("[STEP 1] Creating base BCC cell (3.3 Å, Fe)...")
cmd_create = "atomsk --create bcc 3.3 Fe orient [100] [010] [001] base.xsf"
ret = os.system(cmd_create)
if ret != 0:
    sys.exit("[ERROR] Atomsk command failed at step 1")

print("[STEP 2] Replicating base cell 60×60×60 (~20 nm per side)...")
cmd_duplicate = "atomsk base.xsf -duplicate 60 60 60 base-large.xsf"
ret = os.system(cmd_duplicate)
if ret != 0:
    sys.exit("[ERROR] Atomsk command failed at step 2")

print("[STEP 3] Creating a 30-grain polycrystal -> poly.cfg...")
# For Atomsk version ≥ 0.14, use:
cmd_poly = "atomsk --polycrystal base-large.xsf voronoi_random.txt poly.cfg"
# If using an older version, you might need:
# cmd_poly = "atomsk base-large.xsf polycrystal 30 poly.cfg"
ret = os.system(cmd_poly)
if ret != 0:
    sys.exit("[ERROR] Atomsk command failed at step 3")

print("[INFO] poly.cfg created successfully.")

# -------------
# PYTHON PART: Parse poly.cfg, assign random species, and write LAMMPS data file.
# -------------

input_file  = "poly.cfg"
output_file = "final_lammps.data"

print(f"[STEP 4] Parsing '{input_file}' and writing LAMMPS data file '{output_file}'")

with open(input_file, 'r') as f:
    lines = f.readlines()

header_part = []
atom_part   = []

# Use the marker line "Fe" to separate header from atom lines.
found_fe_line = False
for line_num, line in enumerate(lines):
    if not found_fe_line:
        header_part.append(line)
        if "Fe" in line.strip():
            found_fe_line = True
    else:
        atom_part.append(line)

print(f"[INFO] Header lines: {len(header_part)}")
print(f"[INFO] Atom lines (raw): {len(atom_part)}")

# Parse atom lines (each should have: x y z grainID)
atoms = []
for idx, line in enumerate(atom_part):
    stripped = line.strip()
    if not stripped:
        continue
    parts = stripped.split()
    if len(parts) < 4:
        continue
    try:
        x = float(parts[0])
        y = float(parts[1])
        z = float(parts[2])
        grain_id = parts[3]
    except ValueError:
        continue
    atoms.append({
        'x': x,
        'y': y,
        'z': z,
        'grain': grain_id,
        'type': 0  # to be assigned later
    })

n_atoms = len(atoms)
print(f"[INFO] Parsed {n_atoms} atoms from poly.cfg.")

if n_atoms == 0:
    sys.exit("[ERROR] No atoms parsed. Exiting.")

# Determine if coordinates are fractional by checking the range
xs = [a['x'] for a in atoms]
ys = [a['y'] for a in atoms]
zs = [a['z'] for a in atoms]
x_range = max(xs) - min(xs)
y_range = max(ys) - min(ys)
z_range = max(zs) - min(zs)
print(f"[DEBUG] x_range: {x_range:.6f}, y_range: {y_range:.6f}, z_range: {z_range:.6f}")

# If the range is very small (e.g. <2), assume fractional coordinates and scale.
if x_range < 2 and y_range < 2 and z_range < 2:
    scale_factor = 198.0  # Expected absolute box size (in Å)
    print(f"[INFO] Detected fractional coordinates. Scaling all positions by {scale_factor}.")
    for a in atoms:
        a['x'] *= scale_factor
        a['y'] *= scale_factor
        a['z'] *= scale_factor
else:
    print("[INFO] Coordinates appear to be absolute. No scaling applied.")

# Randomly assign species: split atoms into 5 groups (≈20% each)
fe_indices = list(range(n_atoms))
random.shuffle(fe_indices)

n1 = int(round(0.2 * n_atoms))
n2 = int(round(0.2 * n_atoms))
n3 = int(round(0.2 * n_atoms))
n4 = int(round(0.2 * n_atoms))
# The remainder goes to type 5.
i0, i1 = 0, n1
i2 = i1 + n2
i3 = i2 + n3
i4 = i3 + n4

V_indices  = fe_indices[i0:i1]
Nb_indices = fe_indices[i1:i2]
Ta_indices = fe_indices[i2:i3]
Ti_indices = fe_indices[i3:i4]
Zr_indices = fe_indices[i4:]  # remainder

for idx in V_indices:
    atoms[idx]['type'] = 1  # V
for idx in Nb_indices:
    atoms[idx]['type'] = 2  # Nb
for idx in Ta_indices:
    atoms[idx]['type'] = 3  # Ta
for idx in Ti_indices:
    atoms[idx]['type'] = 4  # Ti
for idx in Zr_indices:
    atoms[idx]['type'] = 5  # Zr

print("[INFO] Random species assignment complete: 1=V, 2=Nb, 3=Ta, 4=Ti, 5=Zr.")

# Set simulation box explicitly (we know the intended dimensions)
xlo, xhi = 0.0, 198.0
ylo, yhi = 0.0, 198.0
zlo, zhi = 0.0, 198.0
print(f"[INFO] Setting simulation box: x:[{xlo}, {xhi}], y:[{ylo}, {yhi}], z:[{zlo}, {zhi}]")

# Hardcode masses (in atomic mass units) for the 5 species:
mass_dict = {
    1: 50.9415,   # V
    2: 92.90638,  # Nb
    3: 180.94788, # Ta
    4: 47.867,    # Ti
    5: 91.224     # Zr
}

with open(output_file, 'w') as f:
    f.write("LAMMPS data file generated by automated_pipeline.py\n\n")
    f.write(f"{n_atoms} atoms\n")
    f.write("5 atom types\n\n")
    
    f.write(f"{xlo:.6f} {xhi:.6f} xlo xhi\n")
    f.write(f"{ylo:.6f} {yhi:.6f} ylo yhi\n")
    f.write(f"{zlo:.6f} {zhi:.6f} zlo zhi\n\n")
    
    f.write("Masses\n\n")
    for t in range(1, 6):
        f.write(f"{t} {mass_dict[t]}\n")
    
    f.write("\nAtoms\n\n")
    # LAMMPS atoms format: "atom-ID atom-type x y z"
    for i, atom in enumerate(atoms, start=1):
        f.write(f"{i} {atom['type']} {atom['x']:.6f} {atom['y']:.6f} {atom['z']:.6f}\n")

print(f"[INFO] Wrote LAMMPS data file to '{output_file}'.")
print("[INFO] Pipeline complete. Load final_lammps.data in OVITO or LAMMPS.")
print("[INFO] Done!")