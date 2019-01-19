# Hydro_blast
"""
Python script for blasting against a protein database and parsing output. 

The script relies heavily on the Biopython package which may be installed directly or through Anaconda (conda install -c bioconda biopython). It is also neccessary for the user to already have the Blast package installed.

This script has been written with particular parameters and purpose in mind which may not apply in most cases. Specifically, one of the output files this script produces is a fasta file which includes the full sequences for queries aligned to subject sequences along with the genes immediately up- and downstream of the hit for all faa files in a given input folder. If adjacent genes are not required, this script is inefficient and should be modified or the use of another script should be considered. The author strongly recommends visiting the Biopython Cookbook for a quick intro to using Blast and parsing output with Python.

"""
