import pandas as pd
import requests
import time

def download_uniprot_sequences_single(csv_file, column_name, output_file, delay=1, max_retries=3):
    """Download FASTA sequences from UniProt one ID at a time to avoid server errors."""
    df = pd.read_csv(csv_file)
    accessions = df[column_name].dropna().astype(str).tolist()

    if not accessions:
        print("No UniProt IDs found.")
        return

    print(f"Total UniProt IDs: {len(accessions)}")

    with open(output_file, 'w') as outfile:
        for idx, acc in enumerate(accessions, 1):
            for attempt in range(max_retries):
                try:
                    url = f"https://rest.uniprot.org/uniprotkb/{acc}.fasta"
                    response = requests.get(url)
                    response.raise_for_status()
                    outfile.write(response.text + '\n')
                    print(f"[{idx}/{len(accessions)}] Downloaded {acc}")
                    time.sleep(delay)
                    break  # success, break retry loop
                except requests.exceptions.RequestException as e:
                    print(f"Attempt {attempt+1} failed for {acc}: {e}")
                    time.sleep(2)  # wait before retry
            else:
                print(f"Failed to download {acc} after {max_retries} attempts.")

    print(f"All sequences saved to {output_file}")


# Usage
input_csv = 'uniprots_no_np.csv'
column_name = 'UniProtKB'
output_fasta_file = 'uniprots_sequences.fasta'

download_uniprot_sequences_single(input_csv, column_name, output_fasta_file)

