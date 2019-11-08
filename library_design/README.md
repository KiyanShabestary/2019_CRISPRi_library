## library_design demo

All files needed for the script are located in the input folder. 

library_design.py and library_design_ncRNAs.py were built to make a sgRNA library targeting ORFs and ncRNAs respectively. The script library_design.py was designed to enable parallelization so that Synechocystis ORFs were first divided into subsets. The script is then run for each subset in parallel before all results are merged together to give the final library. This typically runs overnight. For simplicity, a test_subset is created here for demonstration. The execution of the script takes approximately XX min.

To run a script, open a terminal and go to 2019_CRISPRi_library/library_design, then type:

`python library_design.py test`

Result of this script is located in results folder.

To make the sgRNA library targeting ncRNAs, type:

`python library_design_ncRNAs.py`

Not that this script runs on the whole ncRNA dataset here, so it will take much longer to complete.
