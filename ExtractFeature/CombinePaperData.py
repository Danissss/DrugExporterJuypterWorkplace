import os, sys
import csv

# multiprocessing module in python
import signal
import time
import multiprocessing

# RDKit module
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import rdMolDescriptors
from rdkit.ML.Descriptors import Descriptors
from rdkit.Chem import MACCSkeys
from rdkit.Chem import Descriptors3D
from rdkit.Chem import Lipinski
from rdkit.Chem.rdPartialCharges import ComputeGasteigerCharges
from rdkit.ML.Descriptors.MoleculeDescriptors import MolecularDescriptorCalculator
from ExtractMolecularFeature import *


mol3D = []


def GetCallBackResult(result):
	molObject = result[0]
	propDic   = result[1]
	for key, value in propDic.items():
		molObject.SetProp(key,str(value))
	mol3D.append(molObject)



def Get3DMolFromSMILES(CSVfile):
	"""
	read smiles from csv and save 3d structure to sdf file 
	using 4 processes

	Args:
	Returns:
	Raise:
	"""
	print("==========================================================")
	print(CSVfile)
	MOLFile = CSVfile.replace(".csv","_3DFile.sdf")
	w = Chem.SDWriter(MOLFile)
	p = multiprocessing.Pool()

	mols = []
	with open(CSVfile) as fd:
		rd = csv.reader(fd,delimiter=",")
		for row in rd:

			smiles = row[0]
			if "." in smiles:
				smilesList = smiles.split(".")
				first_smiles = smilesList[0]
				second_smiles = smilesList[1]
				if len(first_smiles) > len(second_smiles):
					smiles = first_smiles
				else:
					smiles = second_smiles

			roles = row[1]
			try:
				mol = GetMolFromSmiles(smiles,dimension=2)
				mol.SetProp("Roles",roles)
				mol.SetProp("SMILES",smiles)
				mols.append(mol)
			except:
				print(smiles)

	

	for mol in mols:
		p.apply_async(Get3DMolFromMol,[mol,mol.GetPropsAsDict()],callback=GetCallBackResult)

	p.close()
	p.join()


	
	for mol in mol3D:
		w.write(mol)
	del mol3D[:]

	w.close()




def GetMolFromMolWithPro(molObject,molPro,dimension=2):
	"""
	Read molecule from mol file and generate 3d coordinate

	Args:
		param1 (string): molFile
	Returns:
		mol object 
	Raise:
		Exceptions
	"""
	try:
		if dimension == 2:
			AllChem.Compute2DCoords(molObject)
			return m2
		else:
			# propDic = molObject.GetPropsAsDict()
			# print(propDic)
			compoundname = molObject.GetName()
			
			m3 = Chem.AddHs(molObject)
			AllChem.EmbedMolecule(m3)
			m3 = Chem.RemoveHs(m3)
			m3.SetName(compoundname)


			# for key, value in molPro.items():
			# 	m3.SetProp(key,value)

			return m3
	except Exception as error:
		print(error)


	return None


def GenerateSMILESCSVFromSubstrateSDF(sdfFile):

	suppl = Chem.SDMolSupplier(sdfFile)
	mols = [x for x in suppl]

	FeatureValueCsv = sdfFile.replace(".sdf","_pgpsubstrate.csv")
	FeatureValueFile = open(FeatureValueCsv, 'w', newline='')
	writer = csv.writer(FeatureValueFile,quoting=csv.QUOTE_ALL)

	for mol in mols:
		substrate  = mol.GetProp("substrate")
		smiles = Chem.MolToSmiles(mol)
		if substrate == "1":
			tmplist = [smiles,"substrate"]
			writer.writerow(tmplist)
		else:
			tmplist = [smiles,"non-substrate"]
			writer.writerow(tmplist)


def GenerateSMILESCSVFromSubstrateSDFBCRP(sdfFile1,sdfFile2):

	suppl = Chem.SDMolSupplier(sdfFile1)
	suppl2 = Chem.SDMolSupplier(sdfFile2)
	mols = [x for x in suppl]
	mols2 = [x for x in suppl2]


	FeatureValueCsv = sdfFile1.replace(".sdf","_BCRPsubstrate.csv")
	FeatureValueFile = open(FeatureValueCsv, 'w', newline='')
	writer = csv.writer(FeatureValueFile,quoting=csv.QUOTE_ALL)

	for mol in mols:
		try:
			smiles = Chem.MolToSmiles(mol)
			tmplist = [smiles,"substrate"]
			writer.writerow(tmplist)
		except:
			continue

	for mol in mols2:
		try:
			smiles = Chem.MolToSmiles(mol)
			tmplist = [smiles,"non-substrate"]
			writer.writerow(tmplist)
		except:
			continue









def main():
	# GenerateSMILESCSVFromSubstrateSDF("/Users/xuan/Desktop/RDKITWorkPlace/NewData/MDR1/pgpsubstrate.sdf")
	GenerateSMILESCSVFromSubstrateSDFBCRP("/Users/xuan/Desktop/RDKITWorkPlace/NewData/BCRP/EszterHazaiBCRPsubstrate.sdf",
										"/Users/xuan/Desktop/RDKITWorkPlace/NewData/BCRP/EszterHazaiBCRPnonsubstrate.sdf")




if __name__ == '__main__':
	main()
