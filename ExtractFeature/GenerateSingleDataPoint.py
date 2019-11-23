import os, sys
import argparse

from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import rdMolDescriptors
from rdkit.ML.Descriptors import Descriptors
from rdkit.Chem import MACCSkeys
from rdkit.Chem import Descriptors3D
from rdkit.Chem import Lipinski
from rdkit.Chem.rdPartialCharges import ComputeGasteigerCharges
from rdkit.ML.Descriptors.MoleculeDescriptors import MolecularDescriptorCalculator

def GetMolFromMol(molFile,dimension=2):
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
			m2 = Chem.MolFromMolFile(molFile)
			AllChem.Compute2DCoords(m2)
			return m2
		else:
			m2 = Chem.MolFromMolFile(molFile)
			AllChem.Compute2DCoords(m2)
			m3 = Chem.AddHs(m2)
			AllChem.EmbedMolecule(m3)
			m3 = Chem.RemoveHs(m3)
			return m3
	except Exception as error:
		print(error)


	return None


def GetMolFromSmiles(smiles,dimension=2):
	"""
	Read molecule from smiles and generate 3d coordinate
	Args:
		param1 (string): smiles
		param2 (int): define the dimension
	Returns:
		mol object 
	Raise:
		Exceptions
	"""
	try:
		if dimension == 2:

			m2 = Chem.MolFromSmiles(smiles)
			AllChem.Compute2DCoords(m2)
			return m2
		else:
			m2 = Chem.MolFromSmiles(smiles)
			AllChem.Compute2DCoords(m2)
			m3 = Chem.AddHs(m2)
			AllChem.EmbedMolecule(m3)
			m3 = Chem.RemoveHs(m3)
			return m3
	except Exception as error:
		print(error)

	return None


def GetMolecularDescriptor(molObject,descriptorName):
	"""
	Read molecule from smiles and generate 3d coordinate
	Args:
		param1 (mol object): rdkit mol object 
		param2 (list): list of descriptor name
	Returns:
		list of descriptor value 
	Raise:
		Exceptions
	"""
	calc = MolecularDescriptorCalculator(descriptorName)
	descrs = calc.CalcDescriptors(molObject)
	return list(descrs)



def GetMACCFP(molObject):
	"""
	Get MACC 166 fingerprint 
	Args:
		param1 (mol object): molObject
	Returns:
		integer array
	Raise:
		No Exceptions
	"""

	maccs = MACCSkeys.GenMACCSKeys(molObject)
	maccs_fingerprint = []
	for i in range(len(maccs)):
		maccs_fingerprint.append(maccs[i])

	return maccs_fingerprint




def WriteDesValueToFile(csvFile):
	"""
	Read molecule from smiles and generate 3d coordinate
	Args:
		param1 (string): smiles
		param2 (int): define the dimension
	Returns:
		mol object 
	Raise:
		Exceptions
	"""
	

	DescValue = GetMolecularDescriptor(molObject,descriptorName)
	print(DescValue)


def GetDescriptorName():
	"""
	Get all descriptor name in list
	Args:

	Returns:
		List
	Raise:
		Exceptions
	"""

	descriptorName = ['BalabanJ','BertzCT','Ipc','HallKierAlpha','Kappa1',
						'Kappa2','Kappa3','Chi0', 'Chi1','Chi0n','Chi1n','Chi2n','Chi3n',
						'Chi4n','Chi0v','Chi1v','Chi2v','Chi3v','Chi4v','MolLogP','MolMR',
						'MolWt','ExactMolWt','HeavyAtomCount','HeavyAtomMolWt','NHOHCount',
						'NOCount','NumHAcceptors','NumHDonors','NumHeteroatoms','NumRotatableBonds',
						'NumValenceElectrons','NumAromaticRings','NumSaturatedRings',
						'NumAliphaticRings','NumAromaticHeterocycles','NumSaturatedHeterocycles',
						'NumAliphaticHeterocycles','NumAromaticCarbocycles','NumSaturatedCarbocycles',
						'NumAliphaticCarbocycles','RingCount','FractionCSP3','TPSA','LabuteASA',]

	for i in range(1,15):
		descriptorName.append('PEOE_VSA{0}'.format(i))
	for i in range(1,11):
		descriptorName.append('SMR_VSA{0}'.format(i))
	for i in range(1,13):
		descriptorName.append('SlogP_VSA{0}'.format(i))
	for i in range(1,12):
		descriptorName.append('EState_VSA{0}'.format(i))
	for i in range(1,11):
		descriptorName.append('VSA_EState{0}'.format(i))
	return descriptorName

def Get3dDescriptor(molObject):
	"""
	Get all 3D descriptor
	Args:

	Returns:
		List
	Raise:
		Exceptions
	"""
	value_list = []
	value_list.append(Descriptors3D.Asphericity(molObject))
	value_list.append(Descriptors3D.Eccentricity(molObject))
	value_list.append(Descriptors3D.InertialShapeFactor(molObject))
	value_list.append(Descriptors3D.NPR1(molObject))
	value_list.append(Descriptors3D.NPR2(molObject))
	value_list.append(Descriptors3D.PMI1(molObject))
	value_list.append(Descriptors3D.PMI2(molObject))
	value_list.append(Descriptors3D.PMI3(molObject))
	value_list.append(Descriptors3D.RadiusOfGyration(molObject))	# Radius of gyration
	value_list.append(Descriptors3D.SpherocityIndex(molObject))		# Spherocity Index
	value_list.append(rdMolDescriptors.CalcPBF(molObject)) 			# Returns the PBF (plane of best fit) descriptor
	value_list += rdMolDescriptors.CalcAUTOCORR3D(molObject)
	value_list += rdMolDescriptors.CalcRDF(molObject)
	value_list += rdMolDescriptors.CalcMORSE(molObject)
	value_list += rdMolDescriptors.CalcWHIM(molObject)
	value_list += rdMolDescriptors.CalcGETAWAY(molObject)

	return value_list

def CalculateStandAloneDescriptor(molObject):
	"""
	Get all standaloneDescriptor
	Args:

	Returns:
		List
	Raise:
		Exceptions
	"""
	value_list = []

	if AllChem.ComputeGasteigerCharges(molObject) == None:
		value_list.append(0.0)
	else:
		value_list.append(1.0)

	value_list.append(rdMolDescriptors.CalcNumAmideBonds(molObject))
	value_list.append(rdMolDescriptors.CalcNumSpiroAtoms(molObject))
	value_list.append(rdMolDescriptors.CalcNumBridgeheadAtoms(molObject))
	value_list += rdMolDescriptors.MQNs_(molObject)

	return value_list



def GenerateAllDescriptor(molObject,descriptorName):
	"""
	Get all descriptor value

	Args:
		param1 (mol object)
		descriptorName(list) generated from GetDescriptorName
	Returns:
		List list of descriptor value
	Raise:
		Exceptions
	"""
	final_list = []
	final_list += GetMACCFP(molObject)
	final_list += GetMolecularDescriptor(molObject,descriptorName)
	final_list += CalculateStandAloneDescriptor(molObject)
	final_list += Get3dDescriptor(molObject)
	return final_list

def Generate2dDescriptor(molObject,descriptorName):
	"""
	Get all descriptor value

	Args:
		param1 (mol object)
		descriptorName(list) generated from GetDescriptorName
	Returns:
		List list of descriptor value
	Raise:
		Exceptions
	"""
	final_list = []
	final_list += GetMACCFP(molObject)
	final_list += GetMolecularDescriptor(molObject,descriptorName)
	final_list += CalculateStandAloneDescriptor(molObject)
	return final_list

def GenerateAllDescriptorHere(molObject,descriptorName):
	"""
	Get all descriptor value

	Args:
		param1 (mol object)
		descriptorName(list) generated from GetDescriptorName
	Returns:
		List list of descriptor value
	Raise:
		Exceptions
	"""
	final_list = []
	final_list += GetMACCFP(molObject)
	final_list += GetMolecularDescriptor(molObject,descriptorName)
	final_list += CalculateStandAloneDescriptor(molObject)
	final_list += Get3dDescriptor(molObject)
	return final_list


def GenerateTrainingSet3DFromMol(sdfFile):
	"""
	Get mol file (contain 3d conformation); calculate the 3d descriptor + all 2d Descriptor + maccs fingerprint
	Have to utilize multiprocessing otherwise some molecule takes forever to run

	"""
	suppl = Chem.SDMolSupplier(sdfFile)
	mols = [x for x in suppl]
	descriptorName = GetDescriptorName()
	return GenerateAllDescriptorHere(mols[0],descriptorName)




def GenerateTrainingSet3DFromSmile(SMILES):
	molObject = Chem.MolFromSmiles(SMILES)
	m3 = Chem.AddHs(molObject)
	AllChem.EmbedMolecule(m3)
	m3 = Chem.RemoveHs(m3)
	descriptorName = GetDescriptorName()
	return GenerateAllDescriptorHere(m3,descriptorName)



def main():
	parser = argparse.ArgumentParser()
	parser.add_argument("-f", "--file", action="store", dest="sdfFile", help="input file with format as sdf")
	parser.add_argument("-s","--smiles",action="store", dest="SMILES", help="standard smiles string")

	args = parser.parse_args()
	# print(args)
	if args.sdfFile != None and args.SMILES != None:
		print("Only take one options {SDF file or SMILES}")

	elif args.sdfFile != None and args.SMILES == None:

		# print(args.sdfFile)
		result = GenerateTrainingSet3DFromMol(args.sdfFile)
		print(result)
	elif args.sdfFile == None and args.SMILES != None:
		# print(args.SMILES)
		result = GenerateTrainingSet3DFromSmile(args.SMILES)
		print(result)
	else:
		parser.print_help()

















if __name__ == '__main__':
	main()