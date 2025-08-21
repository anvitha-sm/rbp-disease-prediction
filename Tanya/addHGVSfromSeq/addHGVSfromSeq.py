import pandas as pd
import requests

# Input / output
input_csv = "ZincFinger_Classical_RBD.csv"
output_csv = "sequences_with_np.csv"

# Load CSV with all columns as strings to avoid DtypeWarning
df = pd.read_csv(input_csv, dtype=str)

def fetch_matching_nps(uniprot_id, sequence):
    """
    Fetch all RefSeq NP_ accessions linked to a UniProt ID
    where the input sequence matches exactly.
    Returns a comma-separated list and description strings.
    """
    url = f"https://rest.uniprot.org/uniprotkb/{uniprot_id}.json"
    try:
        r = requests.get(url, timeout=10)
        if r.status_code != 200:
            return "No NP found", "No description found"
    except requests.RequestException:
        return "No NP found", "No description found"

    data = r.json()
    xrefs = data.get("uniProtKBCrossReferences", [])
    refseq_ids = [x["id"] for x in xrefs if x.get("database") == "RefSeq" and x["id"].startswith("NP_")]

    matching_nps = []
    descriptions = []

    for np_id in refseq_ids:
        np_url = f"https://www.ncbi.nlm.nih.gov/sviewer/viewer.fcgi?id={np_id}&db=protein&report=fasta&retmode=text"
        try:
            np_r = requests.get(np_url, timeout=10)
            if np_r.status_code != 200:
                continue
            lines = np_r.text.splitlines()
            np_seq = "".join([l.strip() for l in lines if not l.startswith(">")])
            # Check for exact match
            if sequence == np_seq:
                matching_nps.append(np_id)
                descriptions.append(f"{np_id} {lines[0][1:]}")
        except requests.RequestException:
            continue

    if matching_nps:
        return ",".join(matching_nps), ",".join(descriptions)
    else:
	return "No NP found", "No description found"

# Process all rows
np_ids = []
descriptions = []

for idx, row in df.iterrows():
    uniprot_id = str(row["uniprot_id"]).strip()
    sequence = str(row["sequence"]).strip()

    np_id, desc = fetch_matching_nps(uniprot_id, sequence)

    np_ids.append(np_id)
    descriptions.append(desc)

    print(f"{uniprot_id}\t{sequence}\t{np_id}")

# Insert NP column as second
df.insert(1, "ProteinHGVS", np_ids)
df.insert(2, "HGVSDescription", descriptions)

# Save CSV
df.to_csv(output_csv, index=False)

print(f"\n  Done! All results saved to {output_csv}")
