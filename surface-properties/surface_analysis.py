from pymol import cmd, util
from collections import defaultdict
import numpy as np
import pandas as pd
import os
import sys
from sklearn.cluster import DBSCAN
import json

import __main__
__main__.pymol_argv = ['pymol','-qc']
import pymol
from pymol import cmd, util
pymol.finish_launching()

def _prepare_structure_for_sasa():
    cmd.dss("all")
    cmd.h_add("all")
    cmd.set("solvent_radius", 1.4)

def _get_per_atom_sasa(selection="all"):
    cmd.get_area(selection, load_b=1)
    model = cmd.get_model(selection)
    return {atom.id: atom.b for atom in model.atom}

def _get_surface_atoms_by_sasa_threshold(sasa_threshold=0.5):
    _prepare_structure_for_sasa()
    per_atom_sasa = _get_per_atom_sasa("all")
    
    surface_atom_ids = [atom_id for atom_id, sasa_val in per_atom_sasa.items() if sasa_val >= sasa_threshold]
    
    if not surface_atom_ids:
        cmd.delete("surface_atoms_temp")
        return "none"
        
    surface_selection_string = f"id {'+'.join(map(str, surface_atom_ids))}"
    cmd.select("surface_atoms_temp", surface_selection_string)
    
    return "surface_atoms_temp"

def get_total_surface_sasa():
    try:
        if cmd.count_atoms("all") == 0:
            return 0.0

        _prepare_structure_for_sasa()
        area = cmd.get_area("all")
        return area
    except Exception:
        return 0.0

def get_chain_surface_sasa(chain="A"):
    try:
        selection_string = f"chain {chain}"
        if cmd.count_atoms(selection_string) == 0:
            return 0.0

        _prepare_structure_for_sasa()
        area = cmd.get_area(selection_string)
        return area
    except Exception:
        return 0.0

def get_surface_residues_by_type_counts(sasa_threshold_for_residue=0.5):
    try:
        _prepare_structure_for_sasa()

        all_unique_residues_data = []
        cmd.iterate("all and polymer", "all_unique_residues_data.append({'chain':chain, 'resn':resn, 'resi':resi})", space={"all_unique_residues_data": all_unique_residues_data})
        
        hydrophobic_res_names = ["ALA", "VAL", "LEU", "ILE", "MET", "PHE", "TRP", "PRO"]
        polar_res_names = ["SER", "THR", "TYR", "ASN", "GLN", "CYS"]
        charged_res_names = ["LYS", "ARG", "HIS", "ASP", "GLU"]
        
        hydrophobic_count = 0
        polar_count = 0
        charged_count = 0
        other_count = 0

        sorted_unique_residues = sorted(list({(r['chain'], r['resn'], r['resi']) for r in all_unique_residues_data}), 
                                        key=lambda x: (x[0], int(x[2])))

        for chain, resn, resi in sorted_unique_residues:
            res_selection_string = f"(chain {chain} and resi {resi})"
            if cmd.count_atoms(res_selection_string) == 0:
                continue

            sasa = cmd.get_area(res_selection_string)
            
            if sasa >= sasa_threshold_for_residue:
                if resn in hydrophobic_res_names:
                    hydrophobic_count += 1
                elif resn in polar_res_names:
                    polar_count += 1
                elif resn in charged_res_names:
                    charged_count += 1
                else:
                    other_count += 1
        
        return hydrophobic_count, polar_count, charged_count, other_count
            
    except Exception:
        return 0, 0, 0, 0

def get_hydrophobic_patch_clusters(eps=4.0, min_samples=2):
    try:
        _prepare_structure_for_sasa()
        all_unique_residues_data = []
        cmd.iterate("all and polymer", "all_unique_residues_data.append({'chain':chain, 'resn':resn, 'resi':resi})", space={"all_unique_residues_data": all_unique_residues_data})
        
        hydrophobic_res_names = ["ALA", "VAL", "LEU", "ILE", "MET", "PHE", "TRP", "PRO"]
        hydrophobic_surface_residues = []
        
        sorted_unique_residues = sorted(list({(r['chain'], r['resn'], r['resi']) for r in all_unique_residues_data}), 
                                        key=lambda x: (x[0], int(x[2])))

        for chain, resn, resi in sorted_unique_residues:
            res_selection_string = f"(chain {chain} and resi {resi})"
            if cmd.count_atoms(res_selection_string) == 0:
                continue
            sasa = cmd.get_area(res_selection_string)
            if sasa >= 0.5 and resn in hydrophobic_res_names:
                hydrophobic_surface_residues.append((chain, resn, resi))

        if not hydrophobic_surface_residues:
            return 0

        coords = []
        for chain, resn, resi in hydrophobic_surface_residues:
            ca_atom_model = cmd.get_model(f"(chain {chain} and resi {resi} and name CA)")
            if ca_atom_model.atom:
                coords.append(ca_atom_model.atom[0].coord)
            else:
                any_atom_model = cmd.get_model(f"(chain {chain} and resi {resi} and not solvent and not H)")
                if any_atom_model.atom:
                    coords.append(any_atom_model.atom[0].coord)
        
        if not coords:
            return 0
            
        X = np.array(coords)
        
        db = DBSCAN(eps=eps, min_samples=min_samples).fit(X)
        labels = db.labels_

        n_clusters = len(set(labels)) - (1 if -1 in labels else 0)
        
        return n_clusters

    except Exception:
        return 0

def get_polar_surface_area(b_factor_threshold=25, sasa_threshold_for_atom=0.1):
    try:
        _prepare_structure_for_sasa()
        
        surface_atoms_selection_name = _get_surface_atoms_by_sasa_threshold(sasa_threshold_for_atom)
        if surface_atoms_selection_name == "none":
            return 0.0

        polar_res_names = "SER+THR+TYR+ASN+GLN+CYS"
        
        selection_name = cmd.get_unused_name("temp_polar_surface_selection")
        cmd.select(selection_name, f"({surface_atoms_selection_name}) and resn {polar_res_names} and b > {b_factor_threshold}")

        if cmd.count_atoms(selection_name) == 0:
            cmd.delete(selection_name)
            cmd.delete(surface_atoms_selection_name)
            return 0.0

        area = cmd.get_area(selection_name)
        cmd.delete(selection_name)
        cmd.delete(surface_atoms_selection_name)
        return area
    except Exception:
        cmd.delete("temp_polar_surface_selection")
        cmd.delete("surface_atoms_temp")
        return 0.0

def get_charged_surface_area(b_factor_threshold=25, sasa_threshold_for_atom=0.1):
    try:
        _prepare_structure_for_sasa()

        surface_atoms_selection_name = _get_surface_atoms_by_sasa_threshold(sasa_threshold_for_atom)
        if surface_atoms_selection_name == "none":
            return 0.0

        charged_res_names = "LYS+ARG+HIS+ASP+GLU"
        
        selection_name = cmd.get_unused_name("temp_charged_surface_selection")
        cmd.select(selection_name, f"({surface_atoms_selection_name}) and resn {charged_res_names} and b > {b_factor_threshold}")

        if cmd.count_atoms(selection_name) == 0:
            cmd.delete(selection_name)
            cmd.delete(surface_atoms_selection_name)
            return 0.0

        area = cmd.get_area(selection_name)
        cmd.delete(selection_name)
        cmd.delete(surface_atoms_selection_name)
        return area
    except Exception:
        cmd.delete("temp_charged_surface_selection")
        cmd.delete("surface_atoms_temp")
        return 0.0

def get_disulfide_surface_sasa(sasa_threshold_for_atom=0.1):
    try:
        _prepare_structure_for_sasa()

        surface_atoms_selection_name = _get_surface_atoms_by_sasa_threshold(sasa_threshold_for_atom)
        if surface_atoms_selection_name == "none":
            return 0.0

        selection_name = cmd.get_unused_name("temp_disulfide_surface_selection")
        cmd.select(selection_name, f"({surface_atoms_selection_name}) and name SG and resn CYS and bonded_to name SG and resn CYS")

        if cmd.count_atoms(selection_name) == 0:
            cmd.delete(selection_name)
            cmd.delete(surface_atoms_selection_name)
            return 0.0

        area = cmd.get_area(selection_name)
        cmd.delete(selection_name)
        cmd.delete(surface_atoms_selection_name)
        return area
    except Exception:
        cmd.delete("temp_disulfide_surface_selection")
        cmd.delete("surface_atoms_temp")
        return 0.0

def get_avg_residue_surface_sasa(sasa_threshold_for_atom=0.1):
    try:
        if cmd.count_atoms("all") == 0:
            return 0.0

        _prepare_structure_for_sasa()
        per_atom_sasa = _get_per_atom_sasa("all")

        residue_sasa_raw = defaultdict(float)
        model = cmd.get_model("all")

        for atom in model.atom:
            if atom.id in per_atom_sasa and per_atom_sasa[atom.id] >= sasa_threshold_for_atom:
                key = (atom.chain, atom.resn, atom.resi)
                residue_sasa_raw[key] += per_atom_sasa[atom.id]
        
        if not residue_sasa_raw:
            return 0.0
            
        total_sasa = sum(residue_sasa_raw.values())
        num_residues = len(residue_sasa_raw)
        
        return total_sasa / num_residues
    except Exception:
        return 0.0

def get_num_residues_with_surface_sasa(sasa_threshold_for_atom=0.1):
    try:
        if cmd.count_atoms("all") == 0:
            return 0

        _prepare_structure_for_sasa()
        per_atom_sasa = _get_per_atom_sasa("all")

        surface_residues_set = set()
        model = cmd.get_model("all")

        for atom in model.atom:
            if atom.id in per_atom_sasa and per_atom_sasa[atom.id] >= sasa_threshold_for_atom:
                key = (atom.chain, atom.resn, atom.resi)
                surface_residues_set.add(key)
        
        return len(surface_residues_set)
    except Exception:
        return 0

max_sasa_dict = {
    "ALA": 129, "ARG": 274, "ASN": 195, "ASP": 193, "CYS": 167,
    "GLN": 223, "GLU": 225, "GLY": 104, "HIS": 224, "ILE": 197,
    "LEU": 201, "LYS": 236, "MET": 224, "PHE": 240, "PRO": 159,
    "SER": 155, "THR": 172, "TRP": 285, "TYR": 263, "VAL": 174
}

def get_avg_relative_surface_sasa(rsa_threshold=0.25, target_chain="A", sasa_threshold_for_atom=0.1):
    try:
        if cmd.count_atoms("all") == 0:
            return 0.0

        _prepare_structure_for_sasa()
        per_atom_sasa = _get_per_atom_sasa("all")
        
        residue_sasa_raw = defaultdict(float)
        model = cmd.get_model("all")

        for atom in model.atom:
            if atom.id in per_atom_sasa and per_atom_sasa[atom.id] >= sasa_threshold_for_atom:
                key = (atom.chain, atom.resn, atom.resi)
                residue_sasa_raw[key] += per_atom_sasa[atom.id]
            
        rsa_values_for_avg = []
        for (chain, resn, resi), sasa in residue_sasa_raw.items():
            if resn in max_sasa_dict and chain == target_chain:
                rsa_val = sasa / max_sasa_dict[resn]
                if rsa_val > rsa_threshold:
                    rsa_values_for_avg.append(rsa_val)
        
        if not rsa_values_for_avg:
            return 0.0
        
        return np.mean(rsa_values_for_avg)
    except Exception:
        return 0.0

def get_num_highly_exposed_residues(rsa_threshold=0.25, target_chain="A", sasa_threshold_for_atom=0.1):
    try:
        if cmd.count_atoms("all") == 0:
            return 0

        _prepare_structure_for_sasa()
        per_atom_sasa = _get_per_atom_sasa("all")
        
        residue_sasa_raw = defaultdict(float)
        model = cmd.get_model("all")

        for atom in model.atom:
            if atom.id in per_atom_sasa and per_atom_sasa[atom.id] >= sasa_threshold_for_atom:
                key = (atom.chain, atom.resn, atom.resi)
                residue_sasa_raw[key] += per_atom_sasa[atom.id]
            
        count = 0
        for (chain, resn, resi), sasa in residue_sasa_raw.items():
            if resn in max_sasa_dict and chain == target_chain:
                rsa_val = sasa / max_sasa_dict[resn]
                if rsa_val > rsa_threshold:
                    count += 1
        
        return count
    except Exception:
        return 0

def get_surface_area_per_residue_type(sasa_threshold=0.1):
    try:
        _prepare_structure_for_sasa()
        per_atom_sasa = _get_per_atom_sasa("all")
        model = cmd.get_model("all")

        residue_area = defaultdict(float)
        for atom in model.atom:
            if atom.id in per_atom_sasa and per_atom_sasa[atom.id] >= sasa_threshold:
                residue_area[atom.resn] += per_atom_sasa[atom.id]
        return dict(residue_area)
    except Exception:
        return {}


def get_avg_rsa_per_chain(target_chain="A", rsa_threshold=0.25, sasa_threshold=0.1):
    try:
        _prepare_structure_for_sasa()
        per_atom_sasa = _get_per_atom_sasa("all")
        model = cmd.get_model("all")
        
        residue_sasa = defaultdict(float)
        for atom in model.atom:
            if atom.id in per_atom_sasa and per_atom_sasa[atom.id] >= sasa_threshold:
                key = (atom.chain, atom.resn, atom.resi)
                residue_sasa[key] += per_atom_sasa[atom.id]
        
        rsa_vals = []
        for (chain, resn, resi), sasa in residue_sasa.items():
            if chain == target_chain and resn in max_sasa_dict:
                rsa = sasa / max_sasa_dict[resn]
                if rsa >= rsa_threshold:
                    rsa_vals.append(rsa)
        
        return np.mean(rsa_vals) if rsa_vals else 0.0
    except Exception:
        return 0.0


def get_count_surface_hbond_atoms(b_factor_threshold=25, sasa_threshold_for_atom=0.1):
    try:
        if cmd.count_atoms("all") == 0:
            return 0, 0
            
        _prepare_structure_for_sasa()
        
        surface_atoms_selection_name = _get_surface_atoms_by_sasa_threshold(sasa_threshold_for_atom)
        if surface_atoms_selection_name == "none":
            return 0, 0

        temp_donor_sel = cmd.get_unused_name("temp_donor_sel")
        temp_acceptor_sel = cmd.get_unused_name("temp_acceptor_sel")

        cmd.select(temp_donor_sel, f"({surface_atoms_selection_name}) and donor and b > {b_factor_threshold}")
        cmd.select(temp_acceptor_sel, f"({surface_atoms_selection_name}) and acceptor and b > {b_factor_threshold}")

        d = cmd.count_atoms(temp_donor_sel)
        a = cmd.count_atoms(temp_acceptor_sel)

        cmd.delete(temp_donor_sel)
        cmd.delete(temp_acceptor_sel)
        cmd.delete(surface_atoms_selection_name)
        return d, a
    except Exception:
        cmd.delete("temp_donor_sel")
        cmd.delete("temp_acceptor_sel")
        cmd.delete("surface_atoms_temp")
        return 0, 0

residue_charge = {
    "ASP": -1, "GLU": -1,
    "LYS": +1, "ARG": +1,
    "HIS": +0.5
}

def get_net_surface_charge(sasa_threshold=0.1):
    try:
        _prepare_structure_for_sasa()
        per_atom_sasa = _get_per_atom_sasa("all")
        model = cmd.get_model("all")
        
        residue_sasa = defaultdict(float)
        for atom in model.atom:
            if atom.id in per_atom_sasa and per_atom_sasa[atom.id] >= sasa_threshold:
                key = (atom.chain, atom.resn, atom.resi)
                residue_sasa[key] += per_atom_sasa[atom.id]
        
        net_charge = 0.0
        for (chain, resn, resi), sasa in residue_sasa.items():
            if resn in residue_charge:
                net_charge += residue_charge[resn]
        
        return net_charge
    except Exception:
        return 0.0


def get_roughness(sasa_threshold_for_atom=0.1):
    try:
        if cmd.count_atoms("all") == 0:
            return 0.0

        _prepare_structure_for_sasa()
        
        surface_atoms_selection_name = _get_surface_atoms_by_sasa_threshold(sasa_threshold_for_atom)
        if surface_atoms_selection_name == "none":
            return 0.0
            
        area = cmd.get_area(surface_atoms_selection_name)
        
        min_xyz, max_xyz = cmd.get_extent("all") 

        dx = max_xyz[0] - min_xyz[0]
        dy = max_xyz[1] - min_xyz[1]
        dz = max_xyz[2] - min_xyz[2]

        volume = dx * dy * dz

        if volume <= 0:
            cmd.delete(surface_atoms_selection_name)
            return 0.0

        roughness_val = area / volume
        cmd.delete(surface_atoms_selection_name)
        return roughness_val
    except Exception:
        cmd.delete("surface_atoms_temp")
        return 0.0

def get_residue_rsa_json(sasa_threshold=0.1):
    try:
        _prepare_structure_for_sasa()
        per_atom_sasa = _get_per_atom_sasa("all")
        model = cmd.get_model("all")
        
        residue_sasa = defaultdict(float)
        rsa_records = []
        
        for atom in model.atom:
            if atom.id in per_atom_sasa and per_atom_sasa[atom.id] >= sasa_threshold:
                key = (atom.chain, atom.resn, atom.resi)
                residue_sasa[key] += per_atom_sasa[atom.id]
        
        for (chain, resn, resi), sasa in residue_sasa.items():
            rsa = sasa / max_sasa_dict[resn] if resn in max_sasa_dict else None
            rsa_records.append({
                "chain": chain,
                "resn": resn,
                "resi": resi,
                "sasa": sasa,
                "rsa": rsa
            })
        
        return json.dumps(rsa_records)
    except Exception:
        return json.dump()

def process_cif_files_to_dataframe(folder_path=".", task_id=1, total_tasks=35):
    results = []
    
    default_sasa_threshold_for_atom = 0.1
    default_sasa_threshold_for_residue = 0.5
    default_b_factor_threshold = 25
    default_rsa_threshold = 0.25
    default_dbscan_eps = 4.0
    default_dbscan_min_samples = 2

    all_cif_files = sorted([f for f in os.listdir(folder_path) if f.endswith(".cif")])
    files_to_process = [f for i, f in enumerate(all_cif_files) if (i % total_tasks) + 1 == task_id]

    print(f"Task {task_id} will process {len(files_to_process)} files.")
    sys.stdout.flush()

    for i, filename in enumerate(files_to_process, 1):
        print(f"[{i}/{len(files_to_process)}] Processing {filename}")
        sys.stdout.flush()

        file_path = os.path.join(folder_path, filename)

        cmd.reinitialize()

        try:
            cmd.load(file_path, filename.split('.')[0])

            if cmd.count_atoms("all") == 0:
                continue

            file_data = {
                "CIF_File_Name": filename,
                "Total_Surface_SASA": get_total_surface_sasa(),
                "Chain_A_Surface_SASA": get_chain_surface_sasa(chain="A"),
                "Chain_B_Surface_SASA": get_chain_surface_sasa(chain="B"),
            }

            hydro_count, polar_count, charged_count, other_count = \
                get_surface_residues_by_type_counts(sasa_threshold_for_residue=default_sasa_threshold_for_residue)
            file_data["Surface_Hydrophobic_Residues_Count"] = hydro_count
            file_data["Surface_Polar_Residues_Count"] = polar_count
            file_data["Surface_Charged_Residues_Count"] = charged_count
            file_data["Surface_Other_Residues_Count"] = other_count

            file_data["Hydrophobic_Patch_Clusters_Count"] = \
                get_hydrophobic_patch_clusters(eps=default_dbscan_eps, min_samples=default_dbscan_min_samples)

            file_data["Polar_Surface_Area_B_gt_25"] = \
                get_polar_surface_area(b_factor_threshold=default_b_factor_threshold, sasa_threshold_for_atom=default_sasa_threshold_for_atom)

            file_data["Charged_Surface_Area_B_gt_25"] = \
                get_charged_surface_area(b_factor_threshold=default_b_factor_threshold, sasa_threshold_for_atom=default_sasa_threshold_for_atom)

            file_data["Disulfide_Surface_SASA"] = \
                get_disulfide_surface_sasa(sasa_threshold_for_atom=default_sasa_threshold_for_atom)

            file_data["Avg_Residue_Surface_SASA"] = \
                get_avg_residue_surface_sasa(sasa_threshold_for_atom=default_sasa_threshold_for_atom)

            file_data["Num_Residues_with_Surface_SASA"] = \
                get_num_residues_with_surface_sasa(sasa_threshold_for_atom=default_sasa_threshold_for_atom)

            file_data["Avg_Relative_Surface_SASA_Chain_A"] = \
                get_avg_relative_surface_sasa(rsa_threshold=default_rsa_threshold, target_chain="A", sasa_threshold_for_atom=default_sasa_threshold_for_atom)

            file_data["Num_Highly_Exposed_Residues_Chain_A"] = \
                get_num_highly_exposed_residues(rsa_threshold=default_rsa_threshold, target_chain="A", sasa_threshold_for_atom=default_sasa_threshold_for_atom)

            donors, acceptors = get_count_surface_hbond_atoms(b_factor_threshold=default_b_factor_threshold, sasa_threshold_for_atom=default_sasa_threshold_for_atom)
            file_data["Surface_Hbond_Donors_Count"] = donors
            file_data["Surface_Hbond_Acceptors_Count"] = acceptors

            file_data["Roughness"] = get_roughness(sasa_threshold_for_atom=default_sasa_threshold_for_atom)
            file_data["Surface_Area_Per_Residue_Type"] = json.dumps(
                get_surface_area_per_residue_type(sasa_threshold=default_sasa_threshold_for_atom)
            )
            file_data["Net_Surface_Charge"] = get_net_surface_charge(sasa_threshold=default_sasa_threshold_for_atom)
            file_data["Avg_RSA_Chain_A"] = get_avg_rsa_per_chain(
                target_chain="A", rsa_threshold=default_rsa_threshold,
                sasa_threshold=default_sasa_threshold_for_atom
            )
            file_data["Residue_RSA_Values_JSON"] = get_residue_rsa_json(
                sasa_threshold=default_sasa_threshold_for_atom
            )
            results.append(file_data)

        except Exception:
            pass
        finally:
            cmd.delete("all")

    if results:
        df = pd.DataFrame(results)
        return df
    else:
        return pd.DataFrame()
    
cif_folder_path = "/u/project/kappel/fraza/CIF Files/AlphaFold/alphafold_pdbs"
task_id_str = os.getenv("SGE_TASK_ID", "1")
try:
    task_id = int(task_id_str)
except ValueError:
    print(f"Invalid task_id value: {task_id_str}, defaulting to 1")
    task_id = 1
analysis_df = process_cif_files_to_dataframe(cif_folder_path, task_id=task_id)
print(analysis_df.columns)
sys.stdout.flush()
print("hi")
sys.stdout.flush()
if not analysis_df.empty:
    output_csv_path = f"protein_surface_analysis_{task_id}.csv"
    analysis_df.to_csv(output_csv_path, index=False)
    
pymol.cmd.quit()
