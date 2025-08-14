import ijson
import pandas as pd
import re

# Paths (renamed to show dedup)
input_json = "/u/home/t/tchhabri/project-kappel/clinvar/full_output_latest_august_2025.json"
input_uniprot_excel = "/u/home/t/tchhabri/project-kappel/with_alphafold_codes.xlsx"
output_excel = "/u/home/t/tchhabri/project-kappel/clinvar/full_rbdpep_uniprot_missense_mutation_dedup.xlsx"
output_uniprot_with_counts = "/u/home/t/tchhabri/project-kappel/clinvar/rbd_pep_with_missense_mutations_count_dedup.xlsx"

# Regex patterns
MC_KEY_RE = re.compile(r"ClassifiedRecord_SimpleAllele_HGVSlist_HGVS_\d+_MolecularConsequence(?:_\d+)?_@Type$")
HGVS_PREFIX_RE = re.compile(r"_MolecularConsequence(?:_\d+)?_@Type$")
PROT_CHANGE_RE = re.compile(r"ClassifiedRecord_SimpleAllele_HGVSlist_HGVS_\d+_ProteinExpression_@change$")
GENE_SYMBOL_KEY_RE = re.compile(r"ClassifiedRecord_SimpleAllele_GeneList_Gene(?:_\d+)?_@Symbol$")

def norm(s: str) -> str:
    s = str(s if s is not None else "").strip().lower()
    return re.sub(r"[\s_]+", " ", s)

def get_variant_genes(flat: dict) -> set:
    genes = set()
    for k, v in flat.items():
        if GENE_SYMBOL_KEY_RE.match(k) and v and str(v).strip().upper() != "NA":
            genes.add(str(v).strip().upper())
    if not genes:
        v = flat.get("ClassifiedRecord_SimpleAllele_GeneList_Gene_@Symbol")
        if v:
            genes.add(str(v).strip().upper())
    return genes

def is_coding_hgvs(prefix: str, flat: dict) -> bool:
    t = norm(flat.get(prefix + "_@Type", ""))
    if "coding" in t:
        return True
    nuc = str(flat.get(prefix + "_NucleotideExpression_@change", ""))
    if nuc.startswith("c."):
        return True
    prot = str(flat.get(prefix + "_ProteinExpression_@change", "")) or str(
        flat.get(prefix + "_ProteinExpression_Expression_#text", "")
    )
    if prot:
        return True
    return False

def has_coding_missense(flat: dict) -> bool:
    for k, v in flat.items():
        if MC_KEY_RE.match(k) and "missense" in norm(v):
            prefix = HGVS_PREFIX_RE.sub("", k)
            if is_coding_hgvs(prefix, flat):
                return True
    return False

def parse_first_coding_missense_protein_change(flat: dict):
    for k, v in flat.items():
        if MC_KEY_RE.match(k) and "missense" in norm(v):
            prefix = HGVS_PREFIX_RE.sub("", k)
            if is_coding_hgvs(prefix, flat):
                prot = str(flat.get(prefix + "_ProteinExpression_@change", ""))
                if prot.startswith("p."):
                    m = re.match(r"p\.([A-Za-z]+)(\d+)([A-Za-z]+)$", prot)
                    if m:
                        return m.group(1), m.group(2), m.group(3)
    return "NA", "NA", "NA"

def get_all_protein_changes(flat: dict) -> str:
    changes = []
    for k, v in flat.items():
        if PROT_CHANGE_RE.match(k):
            if v and v != "NA":
                changes.append(str(v))
    return "; ".join(sorted(set(changes))) if changes else "NA"

def get_indexed_group(flat: dict, base_prefix: str, suffixes):
    results = {}
    for suffix in suffixes:
        vals = []
        pattern = re.compile(re.escape(base_prefix) + r"(\d+)_" + re.escape(suffix) + r"$")
        for key, val in flat.items():
            if pattern.match(key):
                if isinstance(val, list):
                    vals.extend(val)
                else:
                    vals.append(val)
        unique_vals = sorted(set(str(v).strip() for v in vals if v and str(v).strip() != "NA"))
        results[suffix] = "; ".join(unique_vals) if unique_vals else "NA"
    if len(suffixes) == 1:
        return results[suffixes[0]]
    return results

def collect_molecular_consequences(flat: dict) -> str:
    vals = []
    for k, v in flat.items():
        if MC_KEY_RE.match(k) and v and str(v).strip() != "NA":
            vals.append(str(v).strip())
    return "; ".join(sorted(set(vals))) if vals else "NA"

def extract_fields(flat: dict):
    if not has_coding_missense(flat):
        return None
    ref_aa, pos, alt_aa = parse_first_coding_missense_protein_change(flat)
    gene_loc = get_indexed_group(flat, "ClassifiedRecord_SimpleAllele_GeneList_Gene_Location_SequenceLocation_", ["@Chr", "@start", "@stop"])
    var_loc = get_indexed_group(flat, "ClassifiedRecord_SimpleAllele_Location_SequenceLocation_", ["@referenceAlleleVCF", "@alternateAlleleVCF"])
    chromosome = gene_loc.get("@Chr", "NA") if isinstance(gene_loc, dict) else gene_loc
    genomic_start = gene_loc.get("@start", "NA") if isinstance(gene_loc, dict) else "NA"
    genomic_stop = gene_loc.get("@stop", "NA") if isinstance(gene_loc, dict) else "NA"
    ref_allele = var_loc.get("@referenceAlleleVCF", "NA") if isinstance(var_loc, dict) else var_loc
    alt_allele = var_loc.get("@alternateAlleleVCF", "NA") if isinstance(var_loc, dict) else "NA"
    diseases = "; ".join(sorted(set(v for k, v in flat.items() if "ClassifiedConditionList_ClassifiedCondition_#text" in k and v and v != "NA"))) or "NA"
    review_status = "; ".join(sorted(set(v for k, v in flat.items() if "GermlineClassification_ReviewStatus_#text" in k and v and v != "NA"))) or "NA"
    germline_classification = "; ".join(sorted(set(v for k, v in flat.items() if "GermlineClassification_Description_#text" in k and v and v != "NA"))) or "NA"
    condition_list = "; ".join(sorted(set(v for k, v in flat.items() if "GermlineClassification_ConditionList" in k and v and v != "NA"))) or "NA"
    pubmed_ids = "; ".join(sorted(set(v for k, v in flat.items() if "Citation_ID_#text" in k and v and v != "NA"))) or "NA"
    molecular_consequence = collect_molecular_consequences(flat)
    return {
        "VariationID": flat.get("@VariationID", "NA"),
        "Accession": flat.get("@Accession", "NA"),
        "VariationName": flat.get("@VariationName", "NA"),
        "Species": "; ".join(sorted(set(str(v).strip() for k, v in flat.items() if k.startswith("Species_#text") and v))) or "NA",
        "GeneFullName": flat.get("ClassifiedRecord_SimpleAllele_GeneList_Gene_@FullName", "NA"),
        "GeneID": flat.get("ClassifiedRecord_SimpleAllele_GeneList_Gene_@GeneID", "NA"),
        "RelationshipType": flat.get("ClassifiedRecord_SimpleAllele_GeneList_Gene_@RelationshipType", "NA"),
        "Chromosome": chromosome,
        "GenomicStart": genomic_start,
        "GenomicStop": genomic_stop,
        "RefAllele": ref_allele,
        "AltAllele": alt_allele,
        "ProteinChange": get_all_protein_changes(flat),
        "MutatedFrom": ref_aa,
        "ProteinPosition": pos,
        "MutatedTo": alt_aa,
        "MolecularConsequence": molecular_consequence,
        "VariantType": flat.get("@VariationType", "NA"),
        "NumberOfSubmissions": flat.get("@NumberOfSubmissions", "NA"),
        "NumberOfSubmitters": flat.get("@NumberOfSubmitters", "NA"),
        "ReviewStatus": review_status,
        "GermlineClassification": germline_classification,
        "ConditionList": condition_list,
        "Diseases": diseases,
        "PubMedIDs": pubmed_ids,
    }

def stream_variants(json_path: str):
    with open(json_path, 'r', encoding='utf-8') as f:
        for variant in ijson.items(f, 'item'):
            yield variant

def main():
    print("Loading Excel (ProtID, Gene)...")
    df = pd.read_excel(input_uniprot_excel)
    if "ProtID" not in df.columns or "Gene" not in df.columns:
        raise ValueError("Input Excel must have 'ProtID' and 'Gene' columns.")
    gene_to_uniprots = {}
    for _, row in df.iterrows():
        gene = str(row["Gene"]).strip().upper()
        prot = str(row["ProtID"]).strip()
        if gene and gene != "NA" and prot and prot != "NA":
            gene_to_uniprots.setdefault(gene, []).append(prot)
    tracked_genes = set(gene_to_uniprots.keys())
    seen = set()  # (UniProtID, VariationID) to avoid duplicates
    output_rows = []
    missense_counts = {uid: 0 for uid in df["ProtID"].dropna().unique()}
    variant_count = 0
    for flat in stream_variants(input_json):
        variant_count += 1
        if variant_count % 1000 == 0:
            print(f"Processed {variant_count} variants; collected {len(output_rows)} unique matches.")
        genes_in_variant = get_variant_genes(flat)
        matched_genes = genes_in_variant & tracked_genes
        if not matched_genes or not has_coding_missense(flat):
            continue
        rec = extract_fields(flat)
        if not rec:
            continue
        for g in sorted(matched_genes):
            for uid in gene_to_uniprots.get(g, []):
                key = (uid, rec["VariationID"])
                if key in seen:
                    continue
                seen.add(key)
                rec_copy = rec.copy()
                rec_copy["GeneSymbol"] = g
                rec_copy["UniProtID"] = uid
                output_rows.append(rec_copy)
                missense_counts[uid] = missense_counts.get(uid, 0) + 1
    if not output_rows:
        print("No variants found for provided genes/UniProt IDs.")
        return
    df["MissenseVariantCount"] = df["ProtID"].map(missense_counts).fillna(0).astype(int)
    df.to_excel(output_uniprot_with_counts, index=False)
    print(f"Saved updated counts -> {output_uniprot_with_counts}")
    df_out = pd.DataFrame(output_rows)
    df_columns = ["UniProtID", "GeneSymbol", "VariationID", "Accession", "VariationName", "Species", 
                  "GeneFullName", "GeneID", "RelationshipType", "Chromosome", "GenomicStart", "GenomicStop", 
                  "RefAllele", "AltAllele", "ProteinChange", "MutatedFrom", "ProteinPosition", "MutatedTo", 
                  "MolecularConsequence", "VariantType", "NumberOfSubmissions", "NumberOfSubmitters", 
                  "ReviewStatus", "GermlineClassification", "ConditionList", "Diseases", "PubMedIDs"]
    df_out = df_out[[c for c in df_columns if c in df_out.columns]].sort_values(by=["UniProtID", "GeneSymbol", "VariationID"], kind="stable")
    df_out.to_excel(output_excel, index=False)
    print(f"Wrote filtered unique variants -> {output_excel}")
    print("\nMissense variant counts per UniProt ID:")
    for uid, count in sorted(missense_counts.items()):
        print(f"{uid}: {count}")
    print("Done.")

if __name__ == "__main__":
    main()

