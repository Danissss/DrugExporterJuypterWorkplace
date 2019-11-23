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


def Get3DMolFromMol(molObject):
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

		m3 = Chem.AddHs(molObject)
		AllChem.EmbedMolecule(m3)
		m3 = Chem.RemoveHs(m3)
		return m3
	except Exception as error:
		print(error)


	return None

def GetMolFromMol(molObject,dimension=2):
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
			m3 = Chem.AddHs(molObject)
			AllChem.EmbedMolecule(m3)
			m3 = Chem.RemoveHs(m3)
			return m3
	except Exception as error:
		print(error)


	return None


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

def Construct3DMolToFileMultiprocess(fileName,writeFile):
    """
	Read molecule from mol file and generate 3d coordinate and write to file

	Args:
		param1 (string): file for reading
        param2 (string): file for writing
	Returns:
		void 
	Raise:
		Exceptions
	"""
    # Writing sets of molecules

    suppl = Chem.SDMolSupplier(fileName)
    w = Chem.SDWriter(writeFile)
    mols = [x for x in suppl]
    with Pool(processes=2) as pool:
    	for mol in mols:
    		propDic = mol.GetPropsAsDict()
    		result = pool.apply_async(GetMolFromMolWithPro,(mol,propDic,3))
    		
    		try:
    			moleObject = result.get(timeout=100)
    			w.write(moleObject)
    		except:
    			continue

    		

    w.close()


def Construct3DMolToFileMultiprocess2(fileName,writeFile):
    """
	Read molecule from mol file and generate 3d coordinate and write to file

	Args:
		param1 (string): file for reading
        param2 (string): file for writing
	Returns:
		void 
	Raise:
		Exceptions
	"""
    # Writing sets of molecules

    suppl = Chem.SDMolSupplier(fileName)
    w = Chem.SDWriter(writeFile)

    mols = [x for x in suppl]
    p = Pool(processes=2)
    result = p.map(Get3DMolFromMol,mols)
    p.close()
    for i in result:
    	w.write(i)


    w.close()

def Construct3DMolToFile(fileName,writeFile):
    """
	Read molecule from mol file and generate 3d coordinate and write to file
	use timer to set the limit for preventing calculating large molecule
	Args:
		param1 (string): file for reading
        param2 (string): file for writing
	Returns:
		void 
	Raise:
		Exceptions
	"""
    # Writing sets of molecules
    

    w = Chem.SDWriter(writeFile)
    suppl = Chem.SDMolSupplier(fileName)
    mols = [x for x in suppl]
    for mol in mols:
    	# print(mol.GetProp("Solvent"))
    	# print(mol.GetPropNames)
    	signal.signal(signal.SIGALRM, handler)
    	signal.alarm(100)
    	try:
    		mol3d = GetMolFromMol(mol,dimension=3)
    		w.write(mol3d)
    	except Exception:
    		mol3d = mol
    		w.write(mol3d)
    		# print(mol.GetPropsAsDict())


    w.close()




def handler(signum, frame):
	# print("Molecule 3D calculated too long!")
	raise Exception("Molecule 3D calculated too long!")


def main():
    # os.environ["PATH"] += r";C:\Users\me\\Anaconda3\envs\Tensorflow"
    # os.environ["PATH"] += r";C:\Users\me\Anaconda3\envs\Tensorflow\Library\mingw-w64\bin"
    # os.environ["PATH"] += r";C:\Users\me\Anaconda3\envs\Tensorflow\Library\usr\bin"
    # os.environ["PATH"] += r";C:\Users\me\Anaconda3\envs\Tensorflow\Library\bin"
    # os.environ["PATH"] += r";C:\Users\me\Anaconda3\envs\Tensorflow\Scripts"
    absolute_path = "/Users/xuan/Desktop/RDKITWorkPlace/NMRDataSet/test.sdf"
    write_path = "/Users/xuan/Desktop/RDKITWorkPlace/NMRDataSet/test3d.sdf"
    Construct3DMolToFile(absolute_path,write_path)
    # Construct3DMolToFileMultiprocess(absolute_path,write_path)
    # Construct3DMolToFileMultiprocess2(absolute_path,write_path)

if __name__ == '__main__':
	main()
