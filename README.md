# FiBeFTa

An online tool for finding the best fingerprint for a specific target and discern between active and decoy molecules for that target. Available as a terminal tool here and also as a web server at http://fibefta.urv.cat.


## DESCRIPTION:

FiBeFTa compares the performance of 10 different molecular fingerprints, according to their ability to differentiate between an active set of molecules and a decoy one. The main output of the program is the "metrics.csv" file, that compares the performance of all fingerprints according to Enrichment Factor 1 and 10, the AUC and BEDROC. It also generates one file per fingerprint used, with the list of molecules ordered by Tanimoto score with the active molecule they most resemble.

The fingerprints used by FiBeFTa are:

-OpenBabel: FP2, FP3, ChemFP-Substruct

-RDKit: AtomPair, Avalon, Fingerprint, MACCS166, Morgan, Pattern, Torsion


## PREREQUISITES:

The easiest way to install the libraries needed to run a local copy of FiBeFTa is the use of the environment manager "Conda" (https://docs.conda.io/en/latest/). To create a new environment (<env>) that includes all the packages with their respective versions used by FiBeFTa use:

		conda create --name <env> --file req.txt

For the correct calculation of molecular fingerprints the recommended Python version is 2.7. 


## USE:

If an environment has been created with Conda, remember to activate it before the execution by doing:

		conda activate <env>
	
To run the terminal tool use:
   
		python fibefta.py -f active_file.sdf decoy_file.smi
	
The required argument --files (-f) is used to specify the files with the two sets of molecules (actives and decoys). The formats accepted for these files are SDF and SMILES.
  
The optional arguments are --help and --destination. --Destination, or -d, is used to specify a path to save the results. The default option saves the "metrics.csv" in the current location of the terminal, and creates a folder named "FPs", with the list of molecules according to each fingerprint.
