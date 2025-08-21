import pandas as pd
from Bio import Entrez, SeqIO
import time
import os

# === Settings ===
input_xlsx = "full_rbdpep_uniprot_missense_mutation_fixed_v4.xlsx"
input_csv = "full_rbdpep_uniprot_missense_mutation_fixed_v4.csv"
output_file = "clinvar_with_seqs_all.csv"
chunk_size = 500      # number of rows to process at once
cache_limit = 1000    # clear sequence cache after this many entries
Entrez.email = "tchhabria@gmail.com"  # required by NCBI

# Convert Excel to CSV (one-time)
if not os.path.exists(input_csv):
    print("Converting Excel to CSV...")
    df_full = pd.read_excel(input_xlsx)
    df_full.to_csv(input_csv, index=False)
    print("✅ Conversion complete.")

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
    if refseq_id in seq_cache:
        return seq_cache[refseq_id]
    try:
        handle = Entrez.efetch(db="protein", id=refseq_id, rettype="fasta", retmode="text")
        record = SeqIO.read(handle, "fasta")
        seq = str(record.seq)
        desc = record.description  # full description e.g., NP_001092104.1 RNA-binding protein 47 isoform a [Homo sapiens]
        seq_cache[refseq_id] = (seq, desc)
        return seq, desc
    except Exception as e:
        print(f"Error fetching {refseq_id}: {e}")
        return None, None

def mutate_sequence(seq, position, aa_from_3, aa_to_3):
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

# Remove output file if it exists
if os.path.exists(output_file):
    os.remove(output_file)

# Process CSV in chunks
for chunk in pd.read_csv(input_csv, chunksize=chunk_size):
    chunk["GeneIsoform"] = None
    chunk["UnmutatedSeq"] = None
    chunk["MutatedSeq"] = None
    chunk["Error"] = None

    for idx, row in chunk.iterrows():
        isoforms = str(row["ProteinHGVSCode"]).split(";")
        if not isoforms:
            continue
        first_isoform = isoforms[0].strip()

        seq, desc = fetch_sequence(first_isoform)
        if seq:
            mutated_seq, err = mutate_sequence(seq, row["ProteinPosition"], row["MutatedFrom"], row["MutatedTo"])
            chunk.at[idx, "GeneIsoform"] = desc
            chunk.at[idx, "UnmutatedSeq"] = seq
            chunk.at[idx, "MutatedSeq"] = mutated_seq
            chunk.at[idx, "Error"] = err
        else:
            chunk.at[idx, "GeneIsoform"] = first_isoform
            chunk.at[idx, "Error"] = f"fetch_failed: {first_isoform}"

        # Small delay to avoid NCBI throttling
        time.sleep(0.3)

    # Append processed chunk to CSV
    if not os.path.exists(output_file):
        chunk.to_csv(output_file, index=False)
    else:
        chunk.to_csv(output_file, index=False, mode='a', header=False)

    # Clear cache periodically to save memory
    if len(seq_cache) > cache_limit:
        seq_cache.clear()

print(f"✅ Done. All entries processed and saved to {output_file}")

