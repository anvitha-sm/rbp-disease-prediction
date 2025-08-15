import pandas as pd
import numpy as np
import ast
import requests
import io
import os

df = pd.read_csv("Delivery/Dataframe_Chunk_8.csv")

df["Complex_PDBs"] = df["Complex_PDBs"].apply(lambda x: ast.literal_eval(x) if isinstance(x, str) else x)
df["PDB_Strand"] = df["PDB_Strand"].apply(lambda x: ast.literal_eval(x) if isinstance(x, str) else x)

from pymol import cmd

def load_local_cif(cif_path, object_name):
    cmd.delete(object_name)  # just in case
    cmd.load(cif_path, object_name)
    cmd.remove("solvent")

def compute_close_contacts(object_name, chain_id, cutoff=10.0):
    # Select protein CA atoms and RNA atoms
    cmd.select("protein_ca", "{} and polymer.protein and name CA and chain {}".format(object_name, chain_id))
    cmd.select("rna_atoms", "{} and polymer.nucleic".format(object_name))

    protein_model = cmd.get_model("protein_ca")
    rna_model = cmd.get_model("rna_atoms")

    residue_contacts = {}

    for ca_atom in protein_model.atom:
        ca_coord = ca_atom.coord
        closest_dist = None

        for rna_atom in rna_model.atom:
            rna_coord = rna_atom.coord
            dist = ((ca_coord[0] - rna_coord[0])**2 +
                    (ca_coord[1] - rna_coord[1])**2 +
                    (ca_coord[2] - rna_coord[2])**2) ** 0.5

            if closest_dist is None or dist < closest_dist:
                closest_dist = dist

        # Store if within cutoff
        if closest_dist is not None and closest_dist <= cutoff:
            res_key = (ca_atom.resi, ca_atom.resn)
            residue_contacts[res_key] = round(closest_dist, 2)

    cmd.delete("all")

    return [
        {"residue": res[0], "resname": res[1], "min_distance_to_rna": d}
        for res, d in residue_contacts.items()
    ]

def compute_close_contacts(object_name, chain_id, cutoff=10.0):
    # Select all protein atoms and RNA atoms
    cmd.select("protein_atoms", f"{object_name} and polymer.protein and chain {chain_id}")
    cmd.select("rna_atoms", f"{object_name} and polymer.nucleic")

    protein_model = cmd.get_model("protein_atoms")
    rna_model = cmd.get_model("rna_atoms")

    # Group protein atoms by residue
    residue_atoms = {}
    for atom in protein_model.atom:
        key = (atom.resi, atom.resn)
        residue_atoms.setdefault(key, []).append(atom.coord)

    # Compute minimum distance to RNA for each residue
    residue_contacts = {}
    for res_key, atom_coords in residue_atoms.items():
        min_dist = None
        for coord in atom_coords:
            for rna_atom in rna_model.atom:
                rna_coord = rna_atom.coord
                dist = ((coord[0] - rna_coord[0])**2 +
                        (coord[1] - rna_coord[1])**2 +
                        (coord[2] - rna_coord[2])**2) ** 0.5
                if min_dist is None or dist < min_dist:
                    min_dist = dist
        if min_dist is not None and min_dist <= cutoff:
            residue_contacts[res_key] = round(min_dist, 2)

    cmd.delete("all")

    return [
        {"residue": res[0], "resname": res[1], "min_distance_to_rna": d}
        for res, d in residue_contacts.items()
    ]

updated_rows = []

for idx, row in df.iterrows():
    pdb_ids = row["Complex_PDBs"]  # assumed to be a list
    distance_map = {}

    for pdb_id in pdb_ids:
        chain_id = row["PDB_Strand"].get(pdb_id)
        if not chain_id:
            continue

        cif_path = "../../CIF Files/Biological_Assemblies/{}_assembly1.cif".format(pdb_id)
        load_local_cif(cif_path, pdb_id)
        contacts = compute_close_contacts(pdb_id, chain_id, cutoff=1000.0)

        residue_distance_map = {
            entry["residue"]: entry["min_distance_to_rna"]
            for entry in contacts
        }

        distance_map[pdb_id] = residue_distance_map

    row["Closest_Distance_to_RNA_Per_Residue"] = distance_map
    updated_rows.append(row)

# Final updated DataFrame
final_df = pd.DataFrame(updated_rows)
final_df.to_csv("Closest_Output_Chunk_8.csv",index=False)
