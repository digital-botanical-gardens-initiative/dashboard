from rdkit import Chem
import sys

# Get SMILES string from command line argument
smiles_string = sys.argv[1]

# Convert the SMILES string to a molecule object
molecule = Chem.MolFromSmiles(smiles_string)

# Convert the molecule object to an InChI string
inchi_string = Chem.MolToInchiKey(molecule)

print(inchi_string , end='')
