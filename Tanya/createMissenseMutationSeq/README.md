Given an input of HGVS prot code and location of missense mutation it will create new columns for HGVS code used + description of isoform, original seq, and the new mutated seq.

Only handles single point mutations, if there are two missense mutations (ex) Ser243_Hist244_PheAla it will not be able to handle it (so will have to do manually ~.05%
