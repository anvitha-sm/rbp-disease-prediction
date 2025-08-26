#!/usr/bin/env python3
import os
import time
import argparse
import pandas as pd
import requests
from io import StringIO
from Bio import SeqIO

# -------------------
# Amino acid maps
# -------------------
one_to_three = {
    'A':'Ala','R':'Arg','N':'Asn','D':'Asp','C':'Cys','E':'Glu','Q':'Gln','G':'Gly',
    'H':'His','I':'Ile','L':'Leu','K':'Lys','M':'Met','F':'Phe','P':'Pro','S':'Ser',
    'T':'Thr','W':'Trp','Y':'Tyr','V':'Val'
}
three_to_one = {v:k for k,v in one_to_three.items()}

# -------------------
# Argument parser
# -------------------
parser = argparse.ArgumentParser(description="Mutate sequences from CSV")
parser.add_argument("--input", required=True, help="Input CSV")
parser.add_argument("--start", type=int, required=True, help="Start row")
parser.add_argument("--end", type=int, required=True, help="End row")
args = parser.parse_args()

INPUT_CSV = args.input
START = args.start
END = args.end
OUTPUT_FILE = f"gnomad_mutated_{START}_{END}.csv"
NP_CACHE = {}

# -------------------
# Load UniProt FASTA (for fallback)
# -------------------
UNIPROT_FASTA = "uniprots_sequences.fasta"
uniprot_seqs = {}
if os.path.exists(UNIPROT_FASTA):
    for record in SeqIO.parse(UNIPROT_FASTA, "fasta"):
        try:
            uid = record.id.split('|')[1]  # middle field is accession
            uniprot_seqs[uid] = (str(record.seq), record.description)
        except IndexError:
            continue

# -------------------
# Functions
# -------------------
def fetch_np_sequence(refseq_id):
    """Fetch NP sequence from NCBI via direct FASTA URL, cache it"""
    if refseq_id in NP_CACHE:
        return NP_CACHE[refseq_id]

    url = f"https://www.ncbi.nlm.nih.gov/sviewer/viewer.fcgi?id={refseq_id}&db=protein&report=fasta"
    try:
        r = requests.get(url, timeout=10)
        if r.status_code != 200 or not r.text.startswith(">"):
            NP_CACHE[refseq_id] = (None, None)
            return None, None
        handle = StringIO(r.text)
        record = SeqIO.read(handle, "fasta")
        seq = str(record.seq)
        desc = record.description
        NP_CACHE[refseq_id] = (seq, desc)
        time.sleep(0.3)  # avoid hammering NCBI
        return seq, desc
    except Exception:
        NP_CACHE[refseq_id] = (None, None)
        return None, None

def mutate_sequence(seq, position, aa_from_3, aa_to_3):
    """Mutate a sequence at 1-based position if matches expected AA"""
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

def process_row(row):
    used_seq = ""
    orig_seq = ""
    mutated_seq = ""
    error = ""

    # -------------------
    # Try NP sequences first
    # -------------------
    np_field = row.get("ProteinHGVS")
    np_ids = str(np_field).split(",") if pd.notna(np_field) and str(np_field).strip() != "No NP found" else []
    success = False
    for np_id in np_ids:
        np_id = np_id.strip()
        seq, desc = fetch_np_sequence(np_id)
        if seq:
            mutated_seq, err = mutate_sequence(seq, row['AA Position'], row['Previous AA'], row['New AA'])
            if mutated_seq:
                used_seq = f"NP: {desc}"
                orig_seq = seq
                error = err
                success = True
                break
    # -------------------
    # Fallback to UniProt
    # -------------------
    if not success:
        uid = str(row.get("UniProtKB"))
        if uid in uniprot_seqs:
            seq, desc = uniprot_seqs[uid]
            mutated_seq, err = mutate_sequence(seq, row['AA Position'], row['Previous AA'], row['New AA'])
            if mutated_seq:
                used_seq = f"UniProt: {desc}"
                orig_seq = seq
                error = err
            else:
                used_seq = f"UniProt: {desc}"
                orig_seq = seq
                mutated_seq = ""
                error = err or "AA mismatch at position"
        else:
            used_seq = "No NP, UniProt fallback"
            orig_seq = ""
            mutated_seq = ""
            error = "No sequence found"

    return pd.Series([used_seq, orig_seq, mutated_seq, error])

# -------------------
# Main processing
# -------------------
def main():
    print(f"Running batch {START} to {END}")
    df = pd.read_csv(INPUT_CSV, low_memory=False)
    batch_df = df.iloc[START-1:END].copy()

    # Apply row processing
    results = batch_df.apply(process_row, axis=1)
    batch_df[['Used Sequence','Original Sequence','Mutated Sequence','Error']] = results

    # Save batch output
    batch_df.to_csv(OUTPUT_FILE, index=False)
    print(f"Batch {START}-{END} saved to {OUTPUT_FILE}")

if __name__ == "__main__":
    main()

