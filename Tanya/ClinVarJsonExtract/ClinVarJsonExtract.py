#!/usr/bin/env python3
import ijson
import pandas as pd
import re
import json

# --- Paths ---
input_json = "/u/home/t/tchhabri/project-kappel/clinvar/full_output_latest_august_2025.json"
input_uniprot_excel = "/u/home/t/tchhabri/project-kappel/with_alphafold_codes.xlsx"
output_excel = "/u/home/t/tchhabri/project-kappel/clinvar/ClinvarWithAllClassifications4.xlsx"
output_uniprot_with_counts = "/u/home/t/tchhabri/project-kappel/clinvar/ClinvarWithAllClassifications4_missenseCount.xlsx"

# --- Regex patterns ---
MC_KEY_RE = re.compile(r"ClassifiedRecord_SimpleAllele_HGVSlist_HGVS_\d+_MolecularConsequence(?:_\d+)?_@Type$")
HGVS_PREFIX_RE = re.compile(r"_MolecularConsequence(?:_\d+)?_@Type$")
PROT_CHANGE_RE = re.compile(r"ClassifiedRecord_SimpleAllele_HGVSlist_HGVS_\d+_ProteinExpression_@change$")
PROT_EXPR_TEXT_RE = re.compile(r"ClassifiedRecord_SimpleAllele_HGVSlist_HGVS_\d+_ProteinExpression_Expression_#text$")
GENE_SYMBOL_KEY_RE = re.compile(r"ClassifiedRecord_SimpleAllele_GeneList_Gene(?:_\d+)?_@Symbol$")

# --- FIXED regexes: allow optional index ---
RCV_RE = re.compile(r"ClassifiedRecord_RCVList_RCVAccession(?:_(\d+))?_@Accession$")
CLASSIFICATION_RE = re.compile(r"ClassifiedRecord_RCVList_RCVAccession(?:_(\d+))?_RCVClassifications_(\w+?)Classification")

# --- Helpers ---
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

def get_all_protein_hgvs_expressions(flat: dict):
    fulls = []
    codes = []
    for k, v in flat.items():
        if PROT_EXPR_TEXT_RE.match(k):
            if v and v != "NA":
                fulls.append(str(v))
                code = v.split(":", 1)[0] if ":" in v else v
                codes.append(code)
    return ("; ".join(sorted(set(fulls))) if fulls else "NA",
            "; ".join(sorted(set(codes))) if codes else "NA")

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

def extract_gene_info(flat: dict, target_gene_symbol: str, assembly="GRCh38"):
    info = {"GeneFullName": "NA", "GeneID": "NA", "RelationshipType": "NA",
            "Chromosome": "NA", "GenomicStart": "NA", "GenomicStop": "NA"}
    target_norm = target_gene_symbol.strip().upper()
    for k, v in flat.items():
        if k.endswith("@Symbol") and str(v).strip().upper() == target_norm:
            loc_keys = [sk for sk in flat if sk.startswith(k.rsplit("_@Symbol",1)[0]+"_Location_SequenceLocation_")]
            idx = None
            for sk in loc_keys:
                if sk.endswith(f"_@Assembly") and flat[sk] == assembly:
                    m = re.search(r"_Location_SequenceLocation_(\d+)_@Assembly$", sk)
                    if m:
                        idx = m.group(1)
                        break
            if idx is None:
                idx = "0"
            prefix = k.rsplit("_@Symbol", 1)[0]
            info["GeneFullName"] = flat.get(prefix + "_@FullName", "NA")
            info["GeneID"] = flat.get(prefix + "_@GeneID", "NA")
            info["RelationshipType"] = flat.get(prefix + "_@RelationshipType", "NA")
            info["Chromosome"] = flat.get(prefix + f"_Location_SequenceLocation_{idx}_@Chr", "NA")
            info["GenomicStart"] = flat.get(prefix + f"_Location_SequenceLocation_{idx}_@start", "NA")
            info["GenomicStop"] = flat.get(prefix + f"_Location_SequenceLocation_{idx}_@stop", "NA")
            return info
    return info

# --- FIXED classification extraction (robust to indexed/unindexed + multi-conditions) ---
def extract_all_classifications(flat: dict):
    classifications = {
        "Germline": [],
        "SomaticClinicalImpact": [],
        "Oncogenicity": []
    }

    # Detect whether we have indexed RCVs; only use unindexed if no indexed exist
    indexed = set()
    has_unindexed = False
    for k in flat.keys():
        m = RCV_RE.match(k)
        if not m:
            continue
        if m.group(1) is None:
            has_unindexed = True
        else:
            indexed.add(m.group(1))

    if indexed:
        # sort numeric indices
        rcv_indices = sorted(indexed, key=lambda x: int(x))
        use_unindexed = False
    else:
        rcv_indices = [None] if has_unindexed else []
        use_unindexed = has_unindexed

    def base_prefix_for(i):
        # i=None -> unindexed base, else indexed
        return "ClassifiedRecord_RCVList_RCVAccession_" if i is None else f"ClassifiedRecord_RCVList_RCVAccession_{i}_"

    def collect_conditions_for(i):
        base = base_prefix_for(i)
        conds = []
        # Handles ..._ClassifiedCondition_#text and ..._ClassifiedCondition_0_#text, _1_, ...
        cond_pat = re.compile(re.escape(base) + r"ClassifiedConditionList_ClassifiedCondition(?:_\d+)?_#text$")
        for k, v in flat.items():
            if cond_pat.match(k) and v and str(v).strip() != "NA":
                conds.append(str(v).strip())
        return "; ".join(sorted(set(conds))) if conds else "NA"

    # Pre-index the classification types by index to avoid cross-talk
    # key: idx (str or None) -> set of class types for that idx
    idx_to_types = {}
    for k in flat.keys():
        m = CLASSIFICATION_RE.match(k)
        if not m:
            continue
        idx = m.group(1)  # may be None
        # honor unindexed only when we are using unindexed
        if idx is None and not use_unindexed:
            continue
        if idx is not None and idx not in rcv_indices:
            continue
        idx_to_types.setdefault(idx, set()).add(m.group(2))  # Germline, SomaticClinicalImpact, Oncogenicity

    for i in rcv_indices:
        condition = collect_conditions_for(i)
        base = base_prefix_for(i)
        class_types = idx_to_types.get(i, set())

        for class_type in sorted(class_types):
            class_prefix = f"{base}RCVClassifications_{class_type}Classification_"

            review_status = flat.get(class_prefix + "ReviewStatus_#text", "NA")
            description = flat.get(class_prefix + "Description_#text", "NA")
            submission_count = flat.get(class_prefix + "Description_@SubmissionCount", "NA")
            date_last_evaluated = flat.get(class_prefix + "Description_@DateLastEvaluated", "NA")

            record = {
                "condition": condition,
                "review_status": review_status,
                "description": description,
                "submission_count": submission_count
            }
            if date_last_evaluated != "NA":
                record["date_last_evaluated"] = date_last_evaluated

            if class_type == "SomaticClinicalImpact":
                record["clinical_impact_type"] = flat.get(class_prefix + "Description_@ClinicalImpactAssertionType", "NA")
                record["clinical_impact_significance"] = flat.get(class_prefix + "Description_@ClinicalImpactClinicalSignificance", "NA")

            if class_type in classifications:
                classifications[class_type].append(record)

    # Dump lists as JSON strings for Excel cells
    for k in classifications:
        classifications[k] = json.dumps(classifications[k])

    return classifications

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
    all_classifications = extract_all_classifications(flat)
    molecular_consequence = collect_molecular_consequences(flat)
    prot_hgvs_full, prot_hgvs_code = get_all_protein_hgvs_expressions(flat)
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
        "ProteinHGVSFull": prot_hgvs_full,
        "ProteinHGVSCode": prot_hgvs_code,
        "MutatedFrom": ref_aa,
        "ProteinPosition": pos,
        "MutatedTo": alt_aa,
        "MolecularConsequence": molecular_consequence,
        "VariantType": flat.get("@VariationType", "NA"),
        "Germline_Classifications": all_classifications["Germline"],
        "Somatic_Classifications": all_classifications["SomaticClinicalImpact"],
        "Oncogenic_Classifications": all_classifications["Oncogenicity"],
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
    seen = set()
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
                gene_info = extract_gene_info(flat, g, assembly="GRCh38")
                for col, val in gene_info.items():
                    if rec_copy.get(col, "NA") == "NA" and val != "NA":
                        rec_copy[col] = val
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
                  "RefAllele", "AltAllele", "ProteinChange", "ProteinHGVSFull", "ProteinHGVSCode",
                  "MutatedFrom", "ProteinPosition", "MutatedTo", 
                  "MolecularConsequence", "VariantType", "Germline_Classifications", "Somatic_Classifications", "Oncogenic_Classifications"]
    df_out = df_out[[c for c in df_columns if c in df_out.columns]].sort_values(
        by=["UniProtID", "GeneSymbol", "VariationID"], kind="stable"
    )
    df_out.to_excel(output_excel, index=False)
    print(f"Wrote filtered unique variants -> {output_excel}")

    print("\nMissense variant counts per UniProt ID:")
    for uid, count in sorted(missense_counts.items()):
        print(f"{uid}: {count}")
    print("Done.")

if __name__ == "__main__":
    main()

