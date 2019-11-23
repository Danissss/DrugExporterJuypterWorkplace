import os, sys
import csv
from ExtractMolecularFeature import *
from rdkit import Chem
from rdkit.Chem import AllChem
import multiprocessing


descriptor_list = []


def GetAttributeNameFake(size):
	attribute_name = []
	for i in range(size):
		attribute_name.append("Attribute_{0}".format(i))
	attribute_name.append("Class")
	return attribute_name



def GetResultCallBack(result):
	"""
	call back function for apply_asyc

	"""

	descriptor_list.append(result)




def GenerateTrainingSet(csvFile):
	"""



	"""

	FeatureValueCsv = csvFile.replace(".","_value.")
	FeatureValueFile = open(FeatureValueCsv, 'w', newline='')
	writer = csv.writer(FeatureValueFile,quoting=csv.QUOTE_ALL)
	descriptorName = GetDescriptorName()


	attribute_name = GetAttributeNameFake(315)
	writer.writerow(attribute_name)


	# descriptorValue = []
	with open(csvFile) as fd:
		rd = csv.reader(fd, delimiter=",")

		for row in rd:
			smiles = row[0]
			# check the dot
			if "." in smiles:
				smilesList = smiles.split(".")
				first_smiles = smilesList[0]
				second_smiles = smilesList[1]
				if len(first_smiles) > len(second_smiles):
					smiles = first_smiles
				else:
					smiles = second_smiles


			roles = row[1]
			# mol = GetMolFromSmiles(smiles,dimension=3)
			# tmp_list = GenerateAllDescriptor(mol,descriptorName)
			# tmp_list.append(row[1])
			# print(tmp_list)
			# writer.writerow(tmp_list)

			try:
				mol = GetMolFromSmiles(smiles,dimension=2)
				tmp_list = Generate2dDescriptor(mol,descriptorName)
				tmp_list.append(row[1])
				writer.writerow(tmp_list)
			except:
				print(smiles)
				continue


	FeatureValueFile.close()

	return None



def GenerateAllDescriptorHere(molObject,descriptorName,role):
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
	final_list.append(role)
	return final_list


def GenerateTrainingSet3DFromMol(sdfFile):
	"""
	Get mol file (contain 3d conformation); calculate the 3d descriptor + all 2d Descriptor + maccs fingerprint
	Have to utilize multiprocessing otherwise some molecule takes forever to run

	"""
	suppl = Chem.SDMolSupplier(sdfFile)
	mols = [x for x in suppl]

	FeatureValueCsv = sdfFile.replace(".sdf","_3D_descriptor_value.csv")
	FeatureValueFile = open(FeatureValueCsv, 'w', newline='')
	writer = csv.writer(FeatureValueFile,quoting=csv.QUOTE_ALL)
	descriptorName = GetDescriptorName()

	p = multiprocessing.Pool()

	attribute_name = GetAttributeNameFake(1227)
	writer.writerow(attribute_name)
	
	for mol in mols:
		p.apply_async(GenerateAllDescriptorHere,[mol,descriptorName,mol.GetProp("Roles")], callback=GetResultCallBack)
	p.close()
	p.join()

	for i in descriptor_list:
		writer.writerow(i)


	FeatureValueFile.close()
	del descriptor_list[:]

# def GenerateTrainingSet3DFromMol(sdfFile):
# 	"""
# 	Get mol file (contain 3d conformation); calculate the 3d descriptor + all 2d Descriptor + maccs fingerprint
# 	Have to utilize multiprocessing otherwise some molecule takes forever to run

# 	"""
# 	suppl = Chem.SDMolSupplier(sdfFile)
# 	mols = [x for x in suppl]

# 	FeatureValueCsv = sdfFile.replace(".sdf","_3D_descriptor_value.csv")
# 	FeatureValueFile = open(FeatureValueCsv, 'w', newline='')
# 	writer = csv.writer(FeatureValueFile,quoting=csv.QUOTE_ALL)
# 	descriptorName = GetDescriptorName()

# 	p = multiprocessing.Pool()

# 	attribute_name = GetAttributeNameFake(1227)
# 	writer.writerow(attribute_name)

# 	for mol in mols:
# 		try:
# 			roles = mol.GetProp("Roles")
# 			result = p.apply_async(GenerateAllDescriptor,[mol,descriptorName])

# 			tmp_list = result.get(timeout=100)

# 			tmp_list.append(roles)

# 			writer.writerow(tmp_list)
# 		except multiprocessing.TimeoutError:
# 			print(mol.GetProp("SMILES"))
# 			continue

# 	FeatureValueFile.close()







def GenerateTrainingSet3D(csvFile):
	"""



	"""

	FeatureValueCsv = csvFile.replace(".","_3D_value.")
	FeatureValueFile = open(FeatureValueCsv, 'w', newline='')
	writer = csv.writer(FeatureValueFile,quoting=csv.QUOTE_ALL)
	descriptorName = GetDescriptorName()


	attribute_name = GetAttributeNameFake(315)
	writer.writerow(attribute_name)


	# descriptorValue = []
	with open(csvFile) as fd:
		rd = csv.reader(fd, delimiter=",")

		for row in rd:
			smiles = row[0]
			roles = row[1]
			# mol = GetMolFromSmiles(smiles,dimension=3)
			# tmp_list = GenerateAllDescriptor(mol,descriptorName)
			# tmp_list.append(row[1])
			# print(tmp_list)
			# writer.writerow(tmp_list)

			try:
				mol = GetMolFromSmiles(smiles,dimension=2)
				tmp_list = Generate2dDescriptor(mol,descriptorName)
				tmp_list.append(row[1])
				writer.writerow(tmp_list)
			except:
				print(smiles)
				continue


	FeatureValueFile.close()

	return None





def GenerateTestSetForLipid(csvFile):

	FeatureValueCsv = csvFile.replace(".","_value.")
	FeatureValueFile = open(FeatureValueCsv, 'w', newline='')
	writer = csv.writer(FeatureValueFile,quoting=csv.QUOTE_ALL)
	descriptorName = GetDescriptorName()

	attribute_name = GetAttributeNameFake(315)
	writer.writerow(attribute_name)

	with open(csvFile) as fd:
		rd = csv.reader(fd, delimiter=",")

		for row in rd:
			smiles = row[0].replace("\'","")
			roles = "?"
			# mol = GetMolFromSmiles(smiles,dimension=3)
			# tmp_list = GenerateAllDescriptor(mol,descriptorName)
			# tmp_list.append(row[1])
			# print(tmp_list)
			# writer.writerow(tmp_list)

			try:
				mol = GetMolFromSmiles(smiles,dimension=2)
				tmp_list = Generate2dDescriptor(mol,descriptorName)
				tmp_list.append(roles)
				writer.writerow(tmp_list)
			except:
				print(smiles)
				continue


	FeatureValueFile.close()




def main():
	
	GenerateTrainingSet3DFromMol(r'C:\Users\Danis\Desktop\GetDescriptorsBindingData\BCRP_binding_dataset_3DFile.sdf')
	GenerateTrainingSet3DFromMol(r'C:\Users\Danis\Desktop\GetDescriptorsBindingData\MDR1_binding_dataset_3DFile.sdf')
	GenerateTrainingSet3DFromMol(r'C:\Users\Danis\Desktop\GetDescriptorsBindingData\MRP1_binding_dataset_3DFile.sdf')
	GenerateTrainingSet3DFromMol(r'C:\Users\Danis\Desktop\GetDescriptorsBindingData\MRP2_binding_dataset_3DFile.sdf')

















if __name__ == '__main__':
	main()