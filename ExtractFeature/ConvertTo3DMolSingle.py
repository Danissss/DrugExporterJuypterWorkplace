import os
import argparse
import random
from rdkit import Chem
from rdkit.Chem import AllChem






def Get3DMolFromMol(molObject):
	"""
	add 3d confirmation and save to file

	Args:
		param1 (string): molFile
	Returns:
		mol object 
	Raise:
		Exceptions
	"""
	pwd = os.getcwd()

	# generate file name 
	ALPHABET = "0123456789abcdefghijklmnopqrstuvwxyzABCDEFGHIJKLMNOPQRSTUVWXYZ"
	chars=[]
	for i in range(8):
		chars.append(random.choice(ALPHABET))
	file_name = "".join(chars)


	try:
		m3 = Chem.AddHs(molObject)
		AllChem.EmbedMolecule(m3)
		print(Chem.MolToMolBlock(m3),file=open('{0}/{1}.sdf'.format(pwd,file_name),'w+'))

	except:
		m3 = "Error Occur"
		return m3

	m3 = '{0}.sdf'.format(file_name)

	return m3

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("-f", "--file", action="store", dest="sdfFile", help="input file with format as sdf")
    parser.add_argument("-s","--smiles",action="store", dest="SMILES", help="standard smiles string")

    args = parser.parse_args()
    # print(args)
    if args.sdfFile != None and args.SMILES != None:
    	
    	print("Only take one options {SDF file or SMILES}")

    elif args.sdfFile != None and args.SMILES == None:

    	m2 = Chem.MolFromMolFile(args.sdfFile)
    	result = Get3DMolFromMol(m2)
    	print(result)

    elif args.sdfFile == None and args.SMILES != None:
    	
    	m2 = Chem.MolFromSmiles(args.SMILES)
    	result = Get3DMolFromMol(m2)
    	print(result)

    else:
    	parser.print_help()



if __name__ == '__main__':
	main()
