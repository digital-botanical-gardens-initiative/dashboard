import sys
from rdkit import Chem
from rdkit.Chem import rdFMCS

try:
    print("Getting query SMILES from command line arguments", file=sys.stderr)
    query_smiles = sys.argv[1]

    print("Getting database SMILES strings from stdin", file=sys.stderr)
    database_smiles_list = [line.strip() for line in sys.stdin]

    print("Converting the query SMILES string to a molecule object", file=sys.stderr)
    query_molecule = Chem.MolFromSmiles(query_smiles)

    print("Performing the substructure search and storing matching SMILES strings", file=sys.stderr)
    matching_smiles = [smiles for smiles in database_smiles_list 
                       if Chem.MolFromSmiles(smiles).HasSubstructMatch(query_molecule)]

    print(f"Found {len(matching_smiles)} matches", file=sys.stderr)

    print("Converting matching SMILES strings to molecules for MCS finding", file=sys.stderr)
    matching_molecules = [Chem.MolFromSmiles(smiles) for smiles in matching_smiles]

    print("Finding the Maximum Common Substructure (MCS)", file=sys.stderr)
    mcs_result = rdFMCS.FindMCS(matching_molecules)

    print("Printing the MCS to stderr for debugging", file=sys.stderr)
    print(mcs_result.smartsString, file=sys.stderr)

    print("Printing matching SMILES strings to stdout", file=sys.stderr)
    print('\n'.join(matching_smiles))

except Exception as e:
    print(f"Error: {e}", file=sys.stderr)
