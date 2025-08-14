The PymolRBD_PDB_Statistics.py uses the PDB_to_run.xlsx as which PDBs to load in. 
Then it loads a PDB, maps each chain to each Uniprot ID from another spreadsheet which was generated using another script (Uniprot_mappings.py).
It then colors every protein grey and then colors RNA and then colors in the resolved residues green, blue is close to RNA but not a fragment from RBDpep (so not LysC), cyan is N-link but not close,
Yellow is N-link and close to RNA, Magenta is X-link but not close, and Red is X-link and close to RNA. Then the remaining residues are green, meaning they are not close nor from RBDpep data.


Note: Close is defined as 10A from the current script, can be changed.

For the color encoding:

1 - green

2 - blue

3 - cyan

4 - yellow

5 - magenta

6 - red

0 - unresolved.

There are two sets of statistics calculated, one where N-links are onsidered part of the RBD and one with N-links not.

If we count N-links:

Blue: False Neg

Cyan: False Pos

Yellow: True Pos

Magenta: False Pos

Red: True Pos

Green: True Neg


If we do not count N-links as part of the RBD:

Blue: False Neg

Cyan: True Neg

Yellow: False Neg

Magenta: False Pos

Red: True Pos

Green: True neg

TPR, FPR, TNR, FNR were calculated for both conditions.

I would like to note that based on the experiment the X-links (LysC - MS) could either be cross links or could not have been analyzed by the MS for other reasons (like too high ARG or LYS) so be mindful.

The spreadsheet has one row for each cchain of each PDB which contains the following information (parenthesis are the column names): Chain, UniProt, # of Resolved Residues 
(All Residues), # of Near Residues(Near Residues), # of Red Residues, # of Magenta Residues, # of Yellow Residues, # of Cyan Residues, # of Blue Residues, # of Green 
Residues, # of True Positive Residues with N-links counted as part of RBD (True Positive with N-links), # of True Negative Residues with N-links counted as part of RBD 
(True Negative with N-links), # of False Positive Residues with N-links counted as part of RBD (False Positive with N-links), # of False Negative Residues with N-links 
counted as part of RBD (False Negative with N-links), True Postive Rate with N-links counted (True Positive Rate with N-links), True Negative Rate with N-links counted 
(True Negative Rate with N-links), False Postive Rate with N-links counted (False Positive Rate with N-links), False Negative Rate with N-links counted (False Negative Rate 
with N-links), True Postive Rate with N-links not counted as part of RBD(True Positive Rate), True Negative Rate with N-links not counted as part of RBD(True Negative Rate), 
False Postive Rate with N-links not counted as part of RBD(False Positive Rate), False Negative Rate with N-links not counted as part of RBD(False Negative Rate), Red 
Residues Positioning with what Fragment from using start and stop of the LysC fragment (Red Residues Positioning), Residue position (Residue Numbers), Color encoding 
that matches with residue position, so first in color enoding has the position denoted by Residue Numbers (Residue_Colors), # of Unresolved (Grey) Residues
                

There are some empty columns since the uniprot mapping on the PDB is incomplete, meaning some chains are not said to be a certain chain (ex. chain G on 7OKY), so for those
we did not calculate statistics as we are not sure if they would be a uniprot on the RBDpep table or not. It accounts for 243 lines of the total 20,822 lines, so we decided to ignore them as it is 1.16%.

When running on hoffman using the PDB_to_run.xlsx, do in batches of a 100 each (apart from 200-300, split up from 200-240, 240-260, 260-300) in order to finish in 24 hours 
it uses ~ 7GB, but I requested 32 GB since didn't know. Also the code will check if the LysC frag matches up with the AA sequence on the chain, and in some cases it will not
line up, so it outputs the mismatch on the output file, so manually transfer or edit, some fragments do not match really frequently, but not in every case of that fragment (so on some PDBs it doesn't match up, others it does).

I have the current combined output saved as full_pdb_analysis_pymol.xlsx on the project directory in my folder on hoffman2.



