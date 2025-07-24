import pandas as pd
import requests
from collections import defaultdict
from pymol import cmd

SPREADSHEET_PATH = "with_alphafold_codes.xlsx"
PDB_LIST_PATH = "PDB_to_run.xlsx"
OUTPUT_PATH = "batch_pdb_results.xlsx"

three_to_one = {
    'ALA':'A', 'ARG':'R', 'ASN':'N', 'ASP':'D', 'CYS':'C', 'GLN':'Q', 'GLU':'E',
    'GLY':'G', 'HIS':'H', 'ILE':'I', 'LEU':'L', 'LYS':'K', 'MET':'M', 'PHE':'F',
    'PRO':'P', 'SER':'S', 'THR':'T', 'TRP':'W', 'TYR':'Y', 'VAL':'V',
}

def get_chain_to_uniprot_mapping(pdb_id):
    pdb_id = pdb_id.lower()
    url = f"https://www.ebi.ac.uk/pdbe/api/mappings/uniprot/{pdb_id}"
    r = requests.get(url)
    if not r.ok:
        raise ValueError(f"Failed to fetch UniProt mapping for {pdb_id}")
    data = r.json()
    mapping = {}
    if pdb_id not in data:
        raise ValueError(f"No data found for PDB ID {pdb_id}")
    for uniprot_id, entries in data[pdb_id]['UniProt'].items():
        for mapping_entry in entries['mappings']:
            chain = mapping_entry['chain_id']
            mapping[chain] = uniprot_id
    return mapping

def fetch_and_color_full(pdb_id):
    obj_name = pdb_id.lower()
    loaded = False
    try:
        cmd.fetch(pdb_id, name=obj_name, type="pdb1", async_=0)
        if cmd.count_atoms(obj_name) > 0:
            print(f"✅ Loaded PDB (bioassembly): {pdb_id}")
            loaded = True
    except Exception as e:
        print(f"⚠️ Failed to load PDB bioassembly: {e}")

    if not loaded:
        try:
            cmd.fetch(pdb_id, name=obj_name, type="cif", async_=0)
            if cmd.count_atoms(obj_name) > 0:
                print(f"✅ Loaded CIF (asymmetric unit): {pdb_id}")
                loaded = True
        except Exception as e2:
            print(f"❌ Failed to load CIF too: {e2}")

    if not loaded or cmd.count_atoms(obj_name) == 0:
        raise RuntimeError(f"❌ Could not load structure for {pdb_id}")

    try:
        cmd.color("green", f"{obj_name} and polymer.protein")
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

    # Color protein residues near RNA blue
    for uniprot_id, obj_name in uniprot_objects.items():
        near_rna_sel = f"{obj_name} and polymer.protein within {dist} of {rna_obj}"
        cmd.color("blue", near_rna_sel)

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
            all_res = chain_residues[chain]
            near_res = near_rna_residues.get(chain, set())
            global_ms_set = ms_residues_per_protein.get(uniprot_id, set())

            red_res = {
                resi for resi, frag_start, frag_stop, u_id in red_residues_info.get(chain, [])
                if u_id == uniprot_id
            }

            lysc_res = set()
            for _, row in df_pdb.iterrows():
                if row['ProtID'] != uniprot_id:
                    continue
                try:
                    frag_start = int(row['fragmentStart'])
                    frag_stop = int(row['fragmentStop'])
                except Exception:
                    continue
                lysc_res.update(range(frag_start, frag_stop + 1))

            xlink_only = lysc_res - global_ms_set
            magenta_res = xlink_only - near_res
            tn_res = [r for r in all_res if r not in lysc_res and r not in near_res]

            yellow = global_ms_set & near_res & lysc_res
            cyan = global_ms_set - near_res & lysc_res

            close_total = len(red_res | yellow | (near_res - global_ms_set))  # blue residues = near_res - global_ms_set

            blue_residues = near_res - global_ms_set
            fn_percent = (len(blue_residues) / close_total * 100) if close_total else 0.0
            tp_percent = (len(red_res) / close_total * 100) if close_total else 0.0
            fp_percent = (len(magenta_res) / (len(magenta_res) + len(red_res)) * 100) if (len(magenta_res) + len(red_res)) else 0.0
            tn_percent = (len(tn_res) / len(all_res) * 100) if all_res else 0.0
            nlink_close_percent = (len(yellow) / (len(yellow) + len(cyan)) * 100) if (len(yellow) + len(cyan)) else 0.0

            red_res_str = ", ".join([
                f"{resi} (frag {frag_start}-{frag_stop}, UniProt {u_id})"
                for (resi, frag_start, frag_stop, u_id) in red_residues_info.get(chain, [])
                if u_id == uniprot_id
            ]) or "None"

            results.append({
                "PDB_ID": pdb_id,
                "Chain": chain,
                "UniProt": uniprot_id,
                "False_Negative_%": round(fn_percent, 2),
                "True_Positive_%": round(tp_percent, 2),
                "False_Positive_%": round(fp_percent, 2),
                "True_Negative_%": round(tn_percent, 2),
                "N_Link_Close_%": round(nlink_close_percent, 2),
                "Red_Residues": red_res_str
            })

    # Now include chains NOT in uniprot_objects (i.e. no UniProt info) but present in structure
    chains_without_uniprot = chains_in_structure - chains_processed
    for chain in chains_without_uniprot:
        results.append({
            "PDB_ID": pdb_id,
            "Chain": chain,
            "UniProt": "Unknown",
            "False_Negative_%": 0.0,
            "True_Positive_%": 0.0,
            "False_Positive_%": 0.0,
            "True_Negative_%": 0.0,
            "N_Link_Close_%": 0.0,
            "Red_Residues": "None"
        })

    return results

def batch_process_all_pdbs():
    pdb_df = pd.read_excel(PDB_LIST_PATH)
    pdb_list = pdb_df.iloc[:, 0].dropna().unique().tolist()[:3]  # change [:3] to run more/all

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
