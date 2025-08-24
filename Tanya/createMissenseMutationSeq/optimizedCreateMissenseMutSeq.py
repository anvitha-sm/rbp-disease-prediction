import pandas as pd
from Bio import Entrez, SeqIO
import time
import os

# === Settings ===
input_csv = "gnomad_with_refseq_all_rows.csv"
output_file = "clinvar_with_seqs_all.csv"
Entrez.email = "tchhabria@gmail.com"  # required by NCBI

# Amino acid mappings
one_to_three = {
    'A':'Ala','R':'Arg','N':'Asn','D':'Asp','C':'Cys','E':'Glu','Q':'Gln','G':'Gly',
    'H':'His','I':'Ile','L':'Leu','K':'Lys','M':'Met','F':'Phe','P':'Pro','S':'Ser',
    'T':'Thr','W':'Trp','Y':'Tyr','V':'Val'
}
three_to_one = {v:k for k,v in one_to_three.items()}

# Cache for fetched sequences
seq_cache = {}

def fetch_sequence(refseq_id):
    """Fetch protein sequence from NCBI and cache results."""
    if refseq_id in seq_cache:
        return seq_cache[refseq_id]
    try:
        handle = Entrez.efetch(db="protein", id=refseq_id, rettype="fasta", retmode="text")
        record = SeqIO.read(handle, "fasta")
        seq = str(record.seq)
        desc = record.description
        seq_cache[refseq_id] = (seq, desc)
        return seq, desc
    except Exception as e:
        print(f"Error fetching {refseq_id}: {e}")
        seq_cache[refseq_id] = (None, None)
        return None, None

def mutate_sequence(seq, position, aa_from_3, aa_to_3):
    """Perform mutation on the given sequence."""
    aa_from = three_to_one.get(aa_from_3)
    aa_to = three_to_one.get(aa_to_3)
    if not aa_from or not aa_to:
        return None, f"invalid_aa_code: {aa_from_3}->{aa_to_3}"
    idx = int(position) - 1
    if idx < 0 or idx >= len(seq):
        return None, f"out_of_bounds: position {position}, length {len(seq)}"
    if seq[idx] != aa_from:
        return None, f"mismatch: expected {aa_from_3}, found {one_to_three.get(seq[idx],'X')} at {position}"
    mutated_seq = seq[:idx] + aa_to + seq[idx+1:]
    return mutated_seq, None

# === MAIN ===
print("Loading data...")
df = pd.read_csv(input_csv)

# --- Step 1: Collect all unique isoforms ---
all_isoforms = set()
for isoforms in df["ProteinHGVS"].dropna().astype(str):
    first_isoform = isoforms.split(";")[0].strip()
    all_isoforms.add(first_isoform)

print(f"Found {len(all_isoforms)} unique RefSeq IDs. Fetching from NCBI...")

# --- Step 2: Fetch all sequences once ---
for i, refseq_id in enumerate(all_isoforms, 1):
    fetch_sequence(refseq_id)
    if i % 100 == 0:
        print(f"Fetched {i}/{len(all_isoforms)} sequences...")
    time.sleep(0.3)  # throttle NCBI

print("✅ All sequences cached.")

# --- Step 3: Process rows using the cache ---
df["GeneIsoform"] = None
df["UnmutatedSeq"] = None
df["MutatedSeq"] = None
df["Error"] = None

for idx, row in df.iterrows():
    isoforms = str(row["ProteinHGVS"]).split(";")
    if not isoforms:
        continue
    first_isoform = isoforms[0].strip()

    seq, desc = seq_cache.get(first_isoform, (None, None))
    if seq:
        mutated_seq, err = mutate_sequence(seq, row["AA Position"], row["Previous AA"], row["New AA"])
        df.at[idx, "GeneIsoform"] = desc
        df.at[idx, "UnmutatedSeq"] = seq
        df.at[idx, "MutatedSeq"] = mutated_seq
        df.at[idx, "Error"] = err
    else:
        df.at[idx, "GeneIsoform"] = first_isoform
        df.at[idx, "Error"] = f"fetch_failed: {first_isoform}"

print("Saving results...")
df.to_csv(output_file, index=False)
print(f"✅ Done. Results saved to {output_file}")
