from pymol import cmd
from collections import defaultdict
import numpy as np
from sklearn.cluster import DBSCAN

def total_sasa():
    print("Total SASA:", cmd.get_area("all"))

def chain_sasa(chain="A"):
    print(f"Chain {chain} SASA:", cmd.get_area(f"chain {chain}"))

def domain_sasa(name="domain1", resi_range="1-100"):
    cmd.select(name, f"resi {resi_range}")
    print(f"{name} SASA:", cmd.get_area(name))

def hydrophobic_surface_area():
    cmd.dss("all") 
    cmd.h_add("all") 
    cmd.set("solvent_radius", 1.4) 

    surface_residues_selection_parts = []
    sasa_threshold = 0.5
    processed_residues = set() # To store (chain, resi) tuples to avoid duplicates

    all_unique_residues = []
    cmd.iterate("all", "all_unique_residues.append((chain, resi))", space={"all_unique_residues": all_unique_residues})
    
    for chain, resi in sorted(list(set(all_unique_residues))): 
        res_selection_string = f"(chain {chain} and resi {resi})"
        sasa = cmd.get_area(res_selection_string)
        
        if sasa >=sasa_threshold:
            surface_residues_selection_parts.append(res_selection_string)

    if surface_residues_selection_parts:
        cmd.select("surface_residues", " or ".join(surface_residues_selection_parts))
        print(f"Selected {cmd.count_atoms('surface_residues', 'residues')} surface residues (SASA >= {sasa_threshold:.1f} Å²).")
        
        cmd.hide("all")
        cmd.show("cartoon", "all")
        cmd.color("gray70", "all")
        cmd.show("sticks", "surface_residues")
        cmd.color("orange", "surface_residues")
        print("Visualizing surface residues.")
  
    else:
        print("No surface residues found with the given SASA threshold.")
    hydrophobic_res_names = "ALA+VAL+LEU+ILE+MET+PHE+TRP+PRO"
    cmd.select("surface_hydrophobic", f"surface_residues and resn {hydrophobic_res_names}")
    print(f"Selected {cmd.count_atoms('surface_hydrophobic', 'residues')} hydrophobic residues on the surface.")

def hydrophobic_patch_clusters():
    model = cmd.get_model("hydrophobic")
    coords = [(a.coord[0], a.coord[1], a.coord[2]) for a in model.atom]
    if not coords:
        print("No hydrophobic surface residues found.")
        return
    X = np.array(coords)
    db = DBSCAN(eps=4, min_samples=2).fit(X)
    n_clusters = len(set(db.labels_)) - (1 if -1 in db.labels_ else 0)
    print("Hydrophobic patch clusters:", n_clusters)

def polar_surface_area():
    cmd.select("hydrophilic_d1", "domain1 and resn SER+THR+TYR+ASN+GLN+CYS and b > 25")
    print("Polar Surface Area (domain1):", cmd.get_area("hydrophilic_d1"))

def charged_surface_area():
    cmd.select("charged", "resn LYS+ARG+HIS+ASP+GLU and b > 25")
    print("Charged Surface Area:", cmd.get_area("charged"))

def disulfide_sasa():
    cmd.select("disulfide_cys", "name SG and resn CYS and bound_to resn CYS")
    print("Disulfide Bond SASA:", cmd.get_area("disulfide_cys"))

def residue_sasa():
    cmd.get_area("all", load_b=1)
    model = cmd.get_model("all")
    sasa_dict = defaultdict(float)
    for atom in model.atom:
        key = (atom.chain, atom.resn, atom.resi)
        sasa_dict[key] += atom.b
    for key, sasa in sasa_dict.items():
        print(f"Residue {key}: {sasa:.2f} Å²")

max_sasa_dict = {
    "ALA": 129, "ARG": 274, "ASN": 195, "ASP": 193, "CYS": 167,
    "GLN": 223, "GLU": 225, "GLY": 104, "HIS": 224, "ILE": 197,
    "LEU": 201, "LYS": 236, "MET": 224, "PHE": 240, "PRO": 159,
    "SER": 155, "THR": 172, "TRP": 285, "TYR": 263, "VAL": 174
}

def relative_sasa():
    cmd.get_area("all", load_b=1)
    model = cmd.get_model("all")
    residue_sasa = defaultdict(float)
    for atom in model.atom:
        key = (atom.chain, atom.resn, atom.resi)
        residue_sasa[key] += atom.b
    rsa = {}
    for (chain, resn, resi), sasa in residue_sasa.items():
        if resn in max_sasa_dict:
            rsa[(chain, resn, resi)] = sasa / max_sasa_dict[resn]
    surface_exposed = [key for key, val in rsa.items() if val > 0.25 and key[0] == "A"]
    print(f"Surface-exposed residues in chain A (RSA > 0.25): {len(surface_exposed)}")

def count_surface_hbond_atoms():
    cmd.h_add()
    cmd.select("surface_donors", "donor and b > 25")
    cmd.select("surface_acceptors", "acceptor and b > 25")
    d = cmd.count_atoms("surface_donors")
    a = cmd.count_atoms("surface_acceptors")
    print(f"Surface donors: {d}, acceptors: {a}")

def roughness():
    area = cmd.get_area("all")
    min_xyz, max_xyz = cmd.get_extent("all")
    dx, dy, dz = [max_xyz[i] - min_xyz[i] for i in range(3)]
    volume = dx * dy * dz
    print(f"SASA: {area:.2f} Å²")
    print(f"Bounding box volume: {volume:.2f} Å³")
    print(f"Surface roughness: {area/volume:.4f} Å⁻¹")

cmd.extend("total_sasa", total_sasa)
cmd.extend("chain_sasa", chain_sasa)
cmd.extend("domain_sasa", domain_sasa)
cmd.extend("hydrophobic_surface_area", hydrophobic_surface_area)
cmd.extend("hydrophobic_patch_clusters", hydrophobic_patch_clusters)
cmd.extend("polar_surface_area", polar_surface_area)
cmd.extend("charged_surface_area", charged_surface_area)
cmd.extend("disulfide_sasa", disulfide_sasa)
cmd.extend("residue_sasa", residue_sasa)
cmd.extend("relative_sasa", relative_sasa)
cmd.extend("count_surface_hbond_atoms", count_surface_hbond_atoms)
cmd.extend("roughness", roughness)
