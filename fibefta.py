# -*- coding: utf-8 -*-

import sys
import os
from  rdkit.ML.Scoring import Scoring
import chemfp
from chemfp import bitops
import pandas as pd
import sklearn
from sklearn import metrics
from openbabel import pybel
import pathos.multiprocessing
from pathos.multiprocessing import ProcessingPool as Pool
import argparse

pybel.ob.obErrorLog.SetOutputLevel(-1)

fptypes = (
	'OpenBabel-FP2', 'OpenBabel-FP3', 'ChemFP-Substruct-OpenBabel', 'RDKit-AtomPair', 
	'RDKit-Avalon', 'RDKit-Fingerprint', 'RDKit-MACCS166', 'RDKit-Morgan', 
    'RDKit-Pattern', 'RDKit-Torsion' 
	
    )


def crearLlistaTuple(list_of_molecules, esAct):
	
	
	return [(mol, esAct) for mol in list_of_molecules]
	

def ordenarTanimotos(llista_maxims):
	
	llista_maxims.sort(key=lambda x: x[2],reverse = True)
	
	return llista_maxims

def calcularEF(factor, llistaTuplesOrdenada, num_actius):
	
	actius_trobats = 0
	ef_llargada = len(llistaTuplesOrdenada)*factor/100
	
	for i in range(ef_llargada):
		if llistaTuplesOrdenada[i][3]:
			 actius_trobats+=1
	
	a = actius_trobats/float(ef_llargada)
	b = num_actius/float(len(llistaTuplesOrdenada))

	ef = a / b
	
	return ef
	
def calcularBEDROC(llistaTuplesOrdenada):
	
	llista_scores = [(1-el[2], el[3]) for el in llistaTuplesOrdenada]
	bedroc = Scoring.CalcBEDROC(llista_scores, 1, 20) 

	return bedroc

def eliminar_repetits(mfile):
	
	if(mfile.endswith(".sdf")):
		unique_mols = {mol.write("inchi") : mol for mol in pybel.readfile("sdf", mfile)}
		
	elif(mfile.endswith(".smi")):
			unique_mols = {mol.write("inchi") : mol for mol in pybel.readfile("smi", mfile)}
	
	else: raise Exception("Sorry, the formats supported are SDF and SMILES")
			
	outputsdf = pybel.Outputfile("sdf", str(mfile[:-4])+"_uniques.sdf", overwrite=True) 

	for mol in unique_mols.itervalues(): 
		outputsdf.write(mol) 
	
	outputsdf.close() 
		
def funcio_general(fingerprint):
	
	fptype=chemfp.get_fingerprint_type(fingerprint)
	T=fptype.toolkit
	
	mfile=(args.files[0][:-4])+"_uniques.sdf"
	
	with T.read_molecules(mfile) as reader:
		actives=[T.copy_molecule(mol) for mol in reader]

	mfile=(args.files[1][:-4])+"_uniques.sdf"
	
	with T.read_molecules(mfile) as reader:
		decoys=[T.copy_molecule(mol) for mol in reader]

		
	llistaActius = crearLlistaTuple(actives,1)
	llistaTotal = crearLlistaTuple(decoys,0)

	llistaTotal = llistaTotal + llistaActius
	
	df_max = pd.DataFrame()
	df_max["Molècula"] = llistaTotal	
	
	maxims = [[0 for x in range(6)] for y in range(len(llistaTotal))]	
	
	llistaFPA = [None] * len(llistaActius)
	inchi_actius = [None] * len(llistaActius)
	llistaFPT = [None] * len(llistaTotal)
	inchi_total = [None] * len(llistaTotal)

	for a in range(len(llistaActius)):
		llistaFPA[a] = fptype.compute_fingerprint(llistaActius[a][0])
		inchi_actius[a] = T.create_string(llistaActius[a][0],"inchistring")

	for b in range(len(llistaTotal)):
		
		llistaFPT[b] = fptype.compute_fingerprint(llistaTotal[b][0])
		inchi_total[b] = T.create_string(llistaTotal[b][0],"inchistring")

		
	for i in range(len(llistaTotal)):
		
		maxims[i][0] = T.get_id(llistaTotal[i][0])
		maxims[i][1] = T.create_string(llistaTotal[i][0],"smistring")
		maxims[i][2] = 0
		
		for j in range(len(llistaActius)):
			tan=chemfp.bitops.byte_tanimoto(llistaFPT[i],llistaFPA[j])
			if ((tan > maxims[i][2]) and ((inchi_total[i]) != inchi_actius[j])): #Descartar si s'està comparant una molècula amb ella mateixa
				maxims[i][2] = tan
				maxims[i][3] = llistaTotal[i][1]
				maxims[i][4] = T.get_id(llistaActius[j][0])
				maxims[i][5] = T.create_string(llistaActius[j][0],"smistring")
				
	maxims = ordenarTanimotos(maxims)
	
	df_max = pd.DataFrame(maxims, columns =['Molecule ID','Molecule SMILES','Tanimoto', 'Is Active', 'Closest Active ID', 'Closest Active SMILES'])
		
	df_max.to_csv(args.path+'/FPs/'+str(fingerprint)+'.csv', index=False)
	print(str(fingerprint)+" COMPLETED")
	
	metriques[0] = fingerprint
	metriques[1] = round((calcularEF(1,maxims,len(llistaActius))),2)
	metriques[2] = round((calcularEF(10,maxims,len(llistaActius))),2)
	metriques[3] = round(sklearn.metrics.roc_auc_score(df_max['Is Active'],df_max['Tanimoto']),4)
	metriques[4] = round(calcularBEDROC(maxims),4)
	
	return(metriques)
	



if __name__ == '__main__':
	
	
		
	parser=argparse.ArgumentParser()
	requiredArg = parser.add_argument_group('required arguments')
	requiredArg.add_argument('-f', '--files',
						required=True,
						dest='files',
						help='Files for molecule sets in format SDF or SMILES, first specifiying the active file set and then the decoys',
						nargs='+',
						type=str
						)			
	
	parser.add_argument('-d', '--destination',
						dest='path',
						help= 'Specifiy path to save results',
						type=str,
						default='./FiBeFTa/FPs'
						)
											
	args=parser.parse_args()

	directory=os.path.expanduser(args.path)
	if not os.path.exists(directory):
		os.makedirs(directory)
	os.makedirs(directory+'/FPs')

	eliminar_repetits(args.files[0])
	eliminar_repetits(args.files[1])
	
	metriques = [0 for x in range(5)]
	
	#resultats=funcio_general(fptypes[0])
#	pool = multiprocessing.Pool()				
#	resultats=pool.map(funcio_general,fptypes)

	resultats = Pool(len(fptypes)).map(funcio_general,fptypes)
	
	df_met=pd.DataFrame(resultats, columns =['Fingerprint','EF1%','EF10%','AUC','BEDROC'])
	print("\n")
	print(df_met)
	print("\n")

	df_met.to_csv(args.path+'/metrics.csv', index=False)
	
	mfile=args.files[0]
	os.remove((mfile[:-4])+"_uniques.sdf")
	mfile=args.files[1]
	os.remove((mfile[:-4])+"_uniques.sdf")

