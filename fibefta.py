# -*- coding: utf-8 -*-

import sys
from  rdkit.ML.Scoring import Scoring
import chemfp
from chemfp import bitops
import pandas as pd
import sklearn
from sklearn import metrics
from openbabel import pybel
#import multiprocessing
#from multiprocessing import Pool
import pathos.multiprocessing
from pathos.multiprocessing import ProcessingPool as Pool
#import pathos.multiprocessing
#from pathos.multiprocessing import ProcessingPool as Pool


pybel.ob.obErrorLog.SetOutputLevel(-1)
#chemfp.ob.obErrorLog.SetOutputLevel(-1)

fptypes = (
	'OpenBabel-FP2', 'OpenBabel-FP3', 'ChemFP-Substruct-OpenBabel', 'RDKit-AtomPair', 
	'RDKit-Avalon', 'RDKit-Fingerprint', 'RDKit-MACCS166', 'RDKit-Morgan', 
    'RDKit-Pattern', 'RDKit-Torsion' 
	
    )
    
#'RDMACCS-RDKit','RDKit-Morgan', 
       
"""fptypes = (
    'RDKit-Pattern', 'OpenEye-Path', 'OpenBabel-MACCS', 'RDKit-Avalon',
'RDKit-AtomPair', 'RDKit-Fingerprint', 'OpenEye-SMARTSScreen',
'OpenBabel-ECFP2', 'RDKit-SECFP', 'RDKit-Torsion',
'OpenBabel-ECFP8', 'ChemFP-Substruct-RDKit', 'RDMACCS-OpenEye',
'OpenBabel-ECFP6', 'RDMACCS-OpenBabel', 'OpenEye-MDLScreen',
'OpenEye-MACCS166', 'RDMACCS-RDKit', 'OpenBabel-FP4',
'OpenEye-Tree', 'RDKit-Morgan', 'ChemFP-Substruct-OpenEye',
'OpenBabel-FP3', 'OpenBabel-FP2', 'OpenBabel-ECFP0',
'ChemFP-Substruct-OpenBabel', 'OpenEye-Circular',
'OpenBabel-ECFP10', 'OpenBabel-ECFP4', 'OpenEye-MoleculeScreen',
'RDKit-MACCS166'
    )"""


dianes2 = (
'ace', 'ache', 'ada', 'alr2', 'ampc', 'ar', 'cdk2', 'comt', 'cox1',
'cox2', 'dhfr', 'egfr', 'er_agonist', 'er_antagonist', 'fgfr1', 'fxa',
'gart', 'gpb', 'gr', 'hivpr', 'hivrt', 'hmga', 'hsp90', 'inha', 'mr', 'na',
'p38', 'parp', 'pde5', 'pdgfrb', 'pnp', 'ppar_gamma', 'pr', 'rxr_alpha',
'sahh', 'src', 'thrombin', 'tk', 'trypsin', 'vegfr2' 
	)

dianes = (
'mpro'
	)

"""actius = [mol for mol in pybel.readfile("sdf", sys.argv[1])]
problemes = [mol for mol in pybel.readfile("sdf", sys.argv[2])]



print("\nFITXER ACT:\t "+sys.argv[1]+"\t"+str(len(actius))+" molècules")
print("FITXER DEC:\t "+sys.argv[2]+"\t"+str(len(problemes))+" molècules")"""

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
		#mols=[mol for mol in pybel.readfile("sdf", mfile)]
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
	
	#mfile='/home/ori/Ophidian/dud_ligands2006_sdf/'+dianes[i_dianes]+'_ligands.sdf'
	mfile=sys.argv[1]
	
	with T.read_molecules(mfile) as reader:
		actives=[T.copy_molecule(mol) for mol in reader]

	
	#mfile='/home/ori/Ophidian/dud_decoys2006_sdf/'+dianes[i_dianes]+'_decoys.sdf'
	mfile=sys.argv[2]
	
	with T.read_molecules(mfile) as reader:
		decoys=[T.copy_molecule(mol) for mol in reader]

		
	llistaActius = crearLlistaTuple(actives,1)
	llistaTotal = crearLlistaTuple(decoys,0)

	llistaTotal = llistaTotal + llistaActius
	
	df_max = pd.DataFrame()
	df_max["Molècula"] = llistaTotal	
	
	maxims = [[0 for x in range(6)] for y in range(len(llistaTotal))]	
	
	#llistaFPA = [fptype.compute_fingerprint(llistaActius[a][0]) for a in range(len(llistaActius))]
	#llistaFPT = [fptype.compute_fingerprint(llistaTotal[b][0]) for b in range(len(llistaTotal))]
	
	#inchi_actius = [llistaActius[a][0].write("inchi"): mol for mol in 
	llistaFPA = [None] * len(llistaActius)
	inchi_actius = [None] * len(llistaActius)
	llistaFPT = [None] * len(llistaTotal)
	inchi_total = [None] * len(llistaTotal)

	for a in range(len(llistaActius)):
		print(llistaActius[a][0])
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
	
	df_max = pd.DataFrame(maxims, columns =['Molecule ID','Molecule SMILES','Tanimoto', 'Is Active', 'Closer Active ID', 'Closer Active SMILES'])
		
	df_max.to_csv(r'/home/ori/Ophidian/Results/Mpro_'+str(fingerprint)+'.csv')
	print(str(fingerprint)+" COMPLETED")
	
	metriques[0] = fingerprint
	metriques[1] = round((calcularEF(1,maxims,len(llistaActius))),2)
	metriques[2] = round((calcularEF(10,maxims,len(llistaActius))),2)
	metriques[3] = round(sklearn.metrics.roc_auc_score(df_max['Is Active'],df_max['Tanimoto']),4)
	metriques[4] = round(calcularBEDROC(maxims),4)
	
	return(metriques)
	



if __name__ == '__main__':
	
	"""df_ef1 = pd.DataFrame(columns = fptypes)
	df_ef10 = pd.DataFrame(columns = fptypes)
	df_auc = pd.DataFrame(columns = fptypes)
	df_bed = pd.DataFrame(columns = fptypes)
	
	ef1 = [0 for x in range(len(fptypes))]
	ef10 = [0 for x in range(len(fptypes))]
	auc = [0 for x in range(len(fptypes))]
	bedroc = [0 for x in range(len(fptypes))]"""
	

	#for i_dianes in range(len(dianes)):
		
#	eliminar_repetits(sys.argv[1])
	#eliminar_repetits(sys.argv[2])
	
	metriques = [0 for x in range(5)]
	
	#resultats=funcio_general(fptypes[0])
#	pool = multiprocessing.Pool()				
#	resultats=pool.map(funcio_general,fptypes)

	resultats = Pool(len(fptypes)).map(funcio_general,fptypes)
	#print resultats
	#resultats.get()
	
	df_met=pd.DataFrame(resultats, columns =['Fingerprint','EF1%','EF10%','AUC','BEDROC'])
	print(df_met)

	df_met.to_csv(r'/home/ori/Ophidian/Resultats/metrics.csv')

	"""for i in range(len(resultats)):
		ef1[i]=round(resultats[i][1],2)
		ef10[i]=round(resultats[i][2],2)
		auc[i]=round(resultats[i][3],4)
		bedroc[i]=round(resultats[i][4],4)
																
		
		df_ef1.loc[0] = ef1
		df_ef10.loc[0] = ef10
		df_auc.loc[0] = auc
		df_bed.loc[0] = bedroc
		
		#print (dianes[i_dianes]+" done")
	
	df_ef1.to_csv(r'/home/ori/Ophidian/Results2/EF1.csv')
	df_ef10.to_csv(r'/home/ori/Ophidian/Results2/EF10.csv')
	df_auc.to_csv(r'/home/ori/Ophidian/Results2/AUC.csv')
	df_bed.to_csv(r'/home/ori/Ophidian/Results2/BEDROC.csv')"""
	
