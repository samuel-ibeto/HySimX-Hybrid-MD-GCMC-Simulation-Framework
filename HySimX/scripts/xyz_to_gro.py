#!/usr/bin/env python3

# xyz_to_gro.py
# Converts last frame from multi-frame .xyz to .gro format
# Handles multiple molecule types defined in internal dictionary

import sys

# === Configuration: Molecule name ? (atoms per molecule, number of molecules in the system)
molecules = {
    "PIM" : (1902, 10), # this should match the residue name  resname used in the gromacs topology for ease of use when running the md
    "TF"  : (15,  400), 
    "CO2" : (3,   1),
}

box_size = 7.0414   # Box size in nm
input_xyz = "C-2PIM.out.xyz"
output_gro = "output.gro"

# === Step 1: Read XYZ and detect last frame ===
with open(input_xyz, "r") as f:
    lines = f.readlines()

frame_indices = [i for i, line in enumerate(lines) if line.strip().isdigit()]
if not frame_indices:
    print("No valid frames found.")
    sys.exit(1)

last_frame_start = frame_indices[-1]
n_atoms = int(lines[last_frame_start].strip())
frame_data = lines[last_frame_start + 2:last_frame_start + 2 + n_atoms]

# === Step 2: Begin writing GRO file ===
with open(output_gro, "w") as g:
    g.write("Converted from last XYZ frame\n")
    g.write(f"{n_atoms:5d}\n")

    idx = 0
    res_id = 1
    atom_id = 1

    for mol_name, (atoms_per_mol, mol_count) in molecules.items():
        for m in range(mol_count):
            for a in range(atoms_per_mol):
                if idx >= len(frame_data):
                    print(f"Error: ran out of atoms in frame for {mol_name}")
                    sys.exit(1)

                parts = frame_data[idx].split()
                atom_type = parts[0]
                x, y, z = map(float, parts[1:4])
                g.write(f"{res_id:5d}{mol_name:<5}{atom_type:>5}{atom_id:5d}{x/10+box_size/2:8.3f}{y/10+box_size/2:8.3f}{z/10+box_size/2:8.3f}\n")
                atom_id += 1
                idx += 1
            res_id += 1

    g.write(f"{box_size:10.5f}{box_size:10.5f}{box_size:10.5f}\n")

print(f" Done !!! .gro file written with {n_atoms} atoms from last frame.")
