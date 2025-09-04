After downloading a Clinvar VCV .xml file (and unzipping), which I did on hoffman using streamingParse.py, I wanted to extract certain information from the resultant .json. 
This code will output the "UniProtID", "GeneSymbol", "VariationID", "Accession", "VariationName", "Species", 
                  "GeneFullName", "GeneID", "RelationshipType", "Chromosome", "GenomicStart", "GenomicStop", 
                  "RefAllele", "AltAllele", "ProteinChange", "MutatedFrom", "ProteinPosition", "MutatedTo", 
                  "MolecularConsequence", "VariantType", "GermlineClassification", "Somatic_Classifications", "Oncogenic_Classifications:
for each variant on a spreadsheet.

Several variant (each denoted by VariationArchive have several HGVS, this code only looks at the coding mutations and missense. You can change this filtering by 
editing the molecular_consequence_@type =="coding" to nonsense or deletion, but check how the .json has it!
