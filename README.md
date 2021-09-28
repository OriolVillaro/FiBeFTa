# FiBeFTa
An online tool for virtual screening. Avilable as a web server and as a terminal tool

DESCRIPTION:

FiBeFTa compares the performance of 10 different molecular fingerprints, according to their ability to differentiate between an active set of molecules and a decoy one. This program outputs a metrics.csv, comparing the performance of all fingerprints according to Enrichment Factor 1 and 10, the Area Under the Curve and BEDROC. It also generates one file per fingerprint used, with the list of molecules ordered by Tanimoto score with the active molecule they most resemble.

PREREQUISITES:

The only requisite for the use of this software is the environment manager "Conda", in order to create the environment with the file "req.txt". This files includes all the packages and the versions used by fibefta, and cab be used to replicate it by doing:

		conda create --name <env> --file <this file>

USE:

The use of the terminal tool follows the following structure:
   
		python fibefta.py -f active_file.sdf decoy_file.smi
  
This is the minimal instruction, with the required argument --files, to specify the two sets of molecules. The format accepted are SDF and SMILES.
  
The optional arguments are --help and --destination. --Destination, or -d, is used to specify a path to save the results.
