import pandas as pd
import requests
import os
import urllib.request
import json
from collections import defaultdict
from pymol import cmd

SPREADSHEET_PATH = "/u/home/t/tchhabri/project-kappel/with_alphafold_codes.xlsx" #RBDpep Table
PDB_LIST_PATH = "/u/home/t/tchhabri/project-kappel/PDB_to_run.xlsx"
OUTPUT_PATH = "/u/home/t/tchhabri/project-kappel/pdb_analysis_results6.xlsx" #output file 
UNIPROT_MAPPING_PATH = "/u/home/t/tchhabri/project-kappel/uniprot_mappings.txt"  # Your saved file that has the mapping for each chain to what UNIPROT for all the PDBs


three_to_one = {
    'ALA':'A', 'ARG':'R', 'ASN':'N', 'ASP':'D', 'CYS':'C', 'GLN':'Q', 'GLU':'E',
    'GLY':'G', 'HIS':'H', 'ILE':'I', 'LEU':'L', 'LYS':'K', 'MET':'M', 'PHE':'F',
    'PRO':'P', 'SER':'S', 'THR':'T', 'TRP':'W', 'TYR':'Y', 'VAL':'V',
}

def get_chain_to_uniprot_mapping(pdb_id):
    pdb_id = pdb_id.lower()

    # Load the saved JSON data once (you could also do this once globally outside the function)
    with open(UNIPROT_MAPPING_PATH, "r") as f:
        data = json.load(f)

    if pdb_id not in data:
        raise ValueError(f"No data found for PDB ID {pdb_id} in local mappings file")

    mapping = {}
    # The saved JSON structure matches API response, e.g. data[pdb_id][pdb_id]['UniProt']
    uniprot_entries = data[pdb_id][pdb_id].get('UniProt', {})

    for uniprot_id, entries in uniprot_entries.items():
        for mapping_entry in entries['mappings']:
            chain = mapping_entry['chain_id']
            mapping[chain] = uniprot_id

    return mapping

def fetch_and_color_full(pdb_id):
    obj_name = pdb_id.lower()
    loaded = False

    # Try PDB format biological assembly
    try:
        cmd.fetch(pdb_id, name=obj_name, type="pdb1", async_=0)
        if cmd.count_atoms(obj_name) > 0:
            print(f"✅ Loaded PDB (bioassembly): {pdb_id}")
            loaded = True
    except Exception as e:
        print(f"⚠️ Failed to load PDB bioassembly: {e}")

    # If that fails, try CIF format biological assembly from RCSB URL
    if not loaded:
        try:
            cif_url = f"https://files.rcsb.org/download/{pdb_id}-assembly1.cif"
            local_cif = f"{pdb_id}_bioassembly.cif"
            urllib.request.urlretrieve(cif_url, local_cif)
            cmd.load(local_cif, obj_name)
            if cmd.count_atoms(obj_name) > 0:
                print(f"✅ Loaded CIF (bioassembly): {pdb_id}")
                loaded = True
            os.remove(local_cif)
        except Exception as e2:
            print(f"⚠️ Failed to load CIF bioassembly: {e2}")

    # Fallback to asymmetric unit CIF
    if not loaded:
        try:
            cmd.fetch(pdb_id, name=obj_name, type="cif", async_=0)
            if cmd.count_atoms(obj_name) > 0:
                print(f"✅ Loaded CIF (asymmetric unit): {pdb_id}")
                loaded = True
        except Exception as e3:
            print(f"❌ Failed to load asymmetric unit: {e3}")

    # Final check
    if not loaded or cmd.count_atoms(obj_name) == 0:
        raise RuntimeError(f"❌ Could not load structure for {pdb_id}")

    # Apply color and visibility
    try:
        cmd.color("grey", f"{obj_name} and polymer.protein")
        cmd.color("orange", f"{obj_name} and (resn A+C+G+U) and backbone")
        cmd.show("cartoon", f"{obj_name} and polymer.protein")
        cmd.show("cartoon", f"{obj_name} and (resn A+C+G+U)")
        cmd.hide("everything", f"{obj_name} and not (polymer.protein or resn A+C+G+U)")
        cmd.create(f"{obj_name}_RNA", f"{obj_name} and (resn A+C+G+U)")
    except Exception as viz_error:
        print(f"⚠️ Visualization failed: {viz_error}")

def create_uniprot_objects(pdb_id):
    mapping = get_chain_to_uniprot_mapping(pdb_id)
    fetch_and_color_full(pdb_id)

    uniprot_to_chains = defaultdict(list)
    for chain, uniprot in mapping.items():
        uniprot_to_chains[uniprot].append(chain)

    uniprot_objects = {}
    for uniprot_id, chains in uniprot_to_chains.items():
        chain_selector = " + ".join([f"chain {ch}" for ch in chains])
        obj_name = f"{pdb_id}_{uniprot_id}"
        cmd.create(obj_name, f"{pdb_id} and ({chain_selector})")
        uniprot_objects[uniprot_id] = obj_name
        print(f"✅ Created object: {obj_name} from chains: {chains}")
    return uniprot_objects, mapping

def get_real_seq_from_pymol(obj_name, start_resi, stop_resi):
    seq = []
    for resi in range(start_resi, stop_resi + 1):
        sele = f"{obj_name} and resi {resi} and polymer.protein"
        model = cmd.get_model(sele)
        if len(model.atom) == 0:
            seq.append('X')
            continue
        three_letter = model.atom[0].resn.upper()
        one_letter = three_to_one.get(three_letter, 'X')
        seq.append(one_letter)
    return "".join(seq)

def sequences_match(sheet_seq, real_seq):
    if len(sheet_seq) != len(real_seq):
        return False
    for s, r in zip(sheet_seq, real_seq):
        if r == 'X':
            continue
        if s != r:
            return False
    return True

def get_residues_near_rna_per_chain(obj_name, rna_obj, dist=10.0):
    model = cmd.get_model(f"{obj_name} and polymer.protein")
    chain_residues = defaultdict(set)
    for atom in model.atom:
        try:
            chain_residues[atom.chain].add(int(atom.resi))
        except Exception:
            continue

    near_rna_model = cmd.get_model(f"{obj_name} and polymer.protein within {dist} of {rna_obj}")
    near_rna_residues = defaultdict(set)
    for atom in near_rna_model.atom:
        try:
            near_rna_residues[atom.chain].add(int(atom.resi))
        except Exception:
            continue

    return chain_residues, near_rna_residues

def encode_sequence_per_object(obj_name, chain):
    color_index_to_code = {
        3: 1,  # green
        2: 2,  # blue
        5: 3,  # cyan
        4: 6,  # red
        8: 5,  # magenta
        6: 4,  # yellow
        24: 0 #grey
    }

    resi_to_color = {}

    def store_color(resi, color):
        try:
            resi_to_color[int(resi)] = color
        except:
            pass

    cmd.iterate(f"{obj_name} and chain {chain} and name CA",
                "store_color(resi, color)",
                space={"store_color": store_color})

    if not resi_to_color:
        return "", ""

    min_resi = min(resi_to_color)
    max_resi = max(resi_to_color)
    full_resis = list(range(min_resi, max_resi + 1))

    res_nums_str_parts = []
    color_codes_str_parts = []

    prev_was_gap = False

    for r in full_resis:
        if r in resi_to_color:
            res_nums_str_parts.append(str(r))
            color_codes_str_parts.append(str(color_index_to_code.get(resi_to_color[r], 0)))
            prev_was_gap = False
        else:
            if not prev_was_gap:
                res_nums_str_parts.append("...")
                color_codes_str_parts.append("...")
                prev_was_gap = True
            # skip if we're still in a gap

    res_nums_str = " ".join(res_nums_str_parts)
    color_codes_str = " ".join(color_codes_str_parts)



    return res_nums_str, color_codes_str






def process_pdb(pdb_id):
    df = pd.read_excel(SPREADSHEET_PATH)
    uniprot_objects, chain_to_uniprot = create_uniprot_objects(pdb_id)
    rna_obj = f"{pdb_id}_RNA"

    red_residues_info = defaultdict(list)

    def pdb_match(pdb_ids_str):
        return pdb_id.lower() in [p.strip().lower() for p in str(pdb_ids_str).split(';')]

    df_pdb = df[df['PDB_IDs'].apply(pdb_match)]

    ms_residues_per_protein = defaultdict(set)
    for _, row in df_pdb.iterrows():
        uniprot_id = row['ProtID']
        try:
            ms_start = int(row['Start'])
            ms_stop = int(row['Stop'])
        except Exception:
            continue
        ms_residues_per_protein[uniprot_id].update(range(ms_start, ms_stop + 1))

    dist = 10.0
    for uniprot_id, obj_name in uniprot_objects.items():
        chain_residues, near_rna_residues = get_residues_near_rna_per_chain(obj_name, rna_obj, dist=dist)
        #print(f"Near RNA residues for UniProt {uniprot_id} ({obj_name}):")
        #for chain, residues in near_rna_residues.items():
            #print(f"  Chain {chain}: {sorted(residues)} (count: {len(residues)})")
        for chain in chain_residues:
            all_res = chain_residues[chain]
            near_res = near_rna_residues.get(chain, set())

            # Convert to PyMOL selection string (e.g., "resi 79+80+81")
            if all_res:
                resi_str = '+'.join(str(r) for r in sorted(all_res))
                selection = f"({obj_name} and chain {chain} and resi {resi_str})"
                cmd.color ("green", selection)
            if near_res:
                resi_blue = '+'.join(str(r) for r in sorted(near_res))
                blue_selection = f"({obj_name} and chain {chain} and resi {resi_blue})"
                #print(f"Coloring blue: {blue_selection}")
                cmd.color("blue", blue_selection)


    for _, row in df_pdb.iterrows():
        uniprot_id = row['ProtID']
        if uniprot_id not in uniprot_objects:
            print(f"⚠️ UniProt {uniprot_id} not found, skipping row")
            continue
        obj_name = uniprot_objects[uniprot_id]

        try:
            lysc_start = int(row['fragmentStart'])
            lysc_stop = int(row['fragmentStop'])
            ms_start = int(row['Start'])
            ms_stop = int(row['Stop'])
        except Exception:
            print(f"⚠️ Invalid fragment range in row, skipping: {row}")
            continue

        lysc_seq = get_real_seq_from_pymol(obj_name, lysc_start, lysc_stop)
        ms_seq = get_real_seq_from_pymol(obj_name, ms_start, ms_stop)

        lysc_frag_seq = str(row['LysC/ArgC proteolytic fragment (RBDpep)']).upper().replace(" ", "")
        ms_frag_seq = str(row['MS-identified tryptic peptide (N-link)']).upper().replace(" ", "")

        if not sequences_match(lysc_frag_seq, lysc_seq):
            print(f"⚠️ LysC fragment seq mismatch at {lysc_start}-{lysc_stop} UniProt {uniprot_id}")
            continue

        if not sequences_match(ms_frag_seq, ms_seq):
            print(f"⚠️ MS fragment seq mismatch at {ms_start}-{ms_stop} UniProt {uniprot_id}")
            continue

        global_ms_set = ms_residues_per_protein[uniprot_id]
        lysc_range = set(range(lysc_start, lysc_stop + 1))
        exclusive_residues = lysc_range - global_ms_set

        # Color entire LysC fragment cyan
        for resi in lysc_range:
            sel = f"{obj_name} and resi '{resi}'"
            cmd.color("cyan", sel)

        # Color exclusive LysC residues magenta
        for resi in exclusive_residues:
            sel = f"{obj_name} and resi '{resi}'"
            cmd.color("magenta", sel)

        # Color LysC residues near RNA yellow (non-exclusive only)
        for resi in lysc_range:
            if resi in exclusive_residues:
                continue
            sel = f"{obj_name} and resi '{resi}'"
            nearby = cmd.select("tmp_near", f"{sel} within {dist} of {rna_obj}")
            if nearby > 0:
                cmd.color("yellow", sel)

        # Color exclusive LysC residues near RNA red (overrides magenta)
        for resi in exclusive_residues:
            sel = f"{obj_name} and resi '{resi}'"
            nearby = cmd.select("tmp_near", f"{sel} within {dist} of {rna_obj}")
            if nearby > 0:
                cmd.color("red", sel)
                model = cmd.get_model(sel)
                if model.atom:
                    chain = model.atom[0].chain
                    red_residues_info[chain].append((resi, lysc_start, lysc_stop, uniprot_id))

        cmd.deselect()

    results = []
    chains_processed = set()

    # Get all chains in the structure (regardless of uniprot)
    obj_name_all = pdb_id.lower()
    model = cmd.get_model(f"{obj_name_all} and polymer.protein")
    chains_in_structure = set(atom.chain for atom in model.atom)

    # Calculate stats for chains with UniProt objects
    for uniprot_id, obj_name in uniprot_objects.items():
        chain_residues, near_rna_residues = get_residues_near_rna_per_chain(obj_name, rna_obj, dist=dist)
        for chain in chain_residues:
            chains_processed.add(chain)
            # All residues in this chain (e.g., set of ints)
            all_res = chain_residues[chain]

            # Near RNA
            near_res = near_rna_residues.get(chain, set())

            # MS hits for this UniProt (subset of LysC)
            global_ms_set = ms_residues_per_protein.get(uniprot_id, set())

            # LysC-predicted fragment residues
            lysc_res = set()
            for _, row in df_pdb.iterrows():
                if row['ProtID'] != uniprot_id:
                    continue
                try:
                    frag_start = int(row['fragmentStart'])
                    frag_stop = int(row['fragmentStop'])
                    lysc_res.update(range(frag_start, frag_stop + 1))
                except Exception:
                    continue

            # 🔴 Red: near RNA, not in MS, in LysC
            red_res = {
                resi for resi, frag_start, frag_stop, u_id in red_residues_info.get(chain, [])
                if u_id == uniprot_id
            }

            # 🟡 Yellow: near RNA, in MS
            yellow = global_ms_set & near_res

            # 🔵 Blue: near RNA, not in LysC (i.e., not predicted or detected)
            blue_residues = near_res - yellow - red_res

            # 🟣 Magenta: in LysC, not in MS, not near RNA
            magenta_res = lysc_res - global_ms_set - red_res

            # 🔷 Cyan: in MS, not near RNA
            cyan = global_ms_set - near_res

            # 🟢 Green: not in LysC, not in MS, not near RNA
            #excluded = lysc_res.union(near_res).union(global_ms_set)
            green = all_res - red_res - magenta_res - yellow - cyan - blue_residues


            #count N-links close as pos
            tp_residues_n = red_res | yellow
            fp_residues_n = magenta_res | cyan
            fn_residues_n = blue_residues
            tn_residues_n = green

            tp_n = len(tp_residues_n)
            fn_n = len(fn_residues_n)
            fp_n = len(fp_residues_n)
            tn_n = len(tn_residues_n)


            tpr_n = (tp_n / (tp_n + fn_n) * 100) if (tp_n + fn_n) else 0.0    # Sensitivity / Recall
            tnr_n = (tn_n / (tn_n + fp_n) * 100) if (tn_n + fp_n) else 0.0    # Specificity
            fpr_n = (fp_n / (fp_n + tn_n) * 100) if (fp_n + tn_n) else 0.0   # False Positive Rate
            fnr_n = (fn_n / (tp_n + fn_n) * 100) if (tp_n + fn_n) else 0.0   # False Negative Rate

            #ignore N-links treat cyan as green & yellow as blue
            tp_residues = red_res 
            fp_residues = magenta_res 
            fn_residues = blue_residues | yellow
            tn_residues = green | cyan

            #total_close_to_rna = len(tp_residues | fn_residues)

            tp = len(tp_residues)
            fn = len(fn_residues)
            fp = len(fp_residues)
            tn = len(tn_residues)

            tpr = (tp / (tp + fn) * 100) if (tp + fn) else 0.0    # Sensitivity / Recall
            tnr = (tn / (tn + fp) * 100) if (tn + fp) else 0.0    # Specificity
            fpr = (fp / (fp + tn) * 100) if (fp + tn) else 0.0   # False Positive Rate
            fnr = (fn / (tp + fn) * 100) if (tp + fn) else 0.0   # False Negative Rate



            red_res_str = ", ".join([
                f"{resi} (frag {frag_start}-{frag_stop}, UniProt {u_id})"
                for (resi, frag_start, frag_stop, u_id) in red_residues_info.get(chain, [])
                if u_id == uniprot_id
            ]) or "None"

            entry = {
                "PDB_ID": pdb_id,
                "Chain": chain,
                "UniProt": uniprot_id,
                "All Residues": len(all_res),
                "Near Residues": len(near_res),
                "Red Residues": len(red_res),
                "Magenta Residues": len(magenta_res),
                "Yellow Residues": len(yellow),
                "Cyan Residues: ": len(cyan),
                "Blue Residues: ": len(blue_residues),
                "Green Residues: ": len(green),
                "True Positive with N-links": tp_n,
                "True Negative with N-links": tn_n,
                "False Positive with N-links": fp_n,
                "False Negative with N-links": fn_n,
                "False Negative Rate with N-links": round(fnr_n, 2),
                "True Positive Rate with N-links": round(tpr_n, 2),
                "False Positive Rate with N-links": round(fpr_n, 2),
                "True Negative Rate with N-links": round(tnr_n, 2),
                "False Negative Rate": round(fnr, 2),
                "True Positive Rate": round(tpr, 2),
                "False Positive Rate": round(fpr, 2),
                "True Negative Rate": round(tnr, 2),
                "Red Residues Positioning": red_res_str
            }

            # Add residue numbers and colors columns here
            try:
                res_nums_str, color_codes_str = encode_sequence_per_object(obj_name, chain)
            except Exception as e:
                print(f"⚠️ Failed residue/color collection for {uniprot_id} chain {chain}: {e}")
                residue_numbers, residue_colors = [], []

            
            entry["Residue_Numbers"] = res_nums_str
            entry["Residue_Colors"] = color_codes_str
            entry["Unresolved (Grey) Residues"] = color_codes_str.split().count("0")


            results.append(entry)

    # Now include chains NOT in uniprot_objects (i.e. no UniProt info) but present in structure
    chains_without_uniprot = chains_in_structure - chains_processed
    for chain in chains_without_uniprot:
        results.append({
            "PDB_ID": pdb_id,
            "Chain": chain,
            "UniProt": "Unknown",
            "False Negative Rate with N-links": 0.0,
            "True Positive Rate with N-links": 0.0,
            "False Positive Rate with N-links": 0.0,
            "True Negative Rate with N-links": 0.0,
            "False Negative Rate": 0.0,
            "True Positive Rate": 0.0,
            "False Positive Rate": 0.0,
            "True Negative Rate": 0.0,
            "Red Residues Positioning": "None",
            "Residue_Numbers": "",
            "Residue_Colors": ""
        })


    return results

def batch_process_all_pdbs():
    pdb_df = pd.read_excel(PDB_LIST_PATH)
    pdb_list = pdb_df.iloc[:, 0].dropna().unique().tolist()  # change [:1] to run more/all

    all_results = []
    for pdb_id in pdb_list:
        print(f"Processing {pdb_id} ...")
        try:
            cmd.reinitialize()
            cmd.delete("all")
            results = process_pdb(pdb_id)
            all_results.extend(results)
        except Exception as e:
            print(f"Error processing {pdb_id}: {e}")

    if all_results:
        out_df = pd.DataFrame(all_results)
        out_df.to_excel(OUTPUT_PATH, index=False)
        print(f"Saved combined results to {OUTPUT_PATH}")
    else:
        print("No results to save.")

if __name__ == "__main__":
    batch_process_all_pdbs()

