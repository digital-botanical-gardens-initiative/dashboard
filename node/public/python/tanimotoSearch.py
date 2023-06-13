import itertools
import pandas as pd
import numpy as np
from rdkit.Chem import AllChem, DataStructs, rdFMCS
from rdkit import Chem
import sys


try:
    print("Getting query SMILES from command line arguments", file=sys.stderr)
    query_smiles = sys.argv[1]
    threshold_similarity = float(sys.argv[2])

    print("Getting database SMILES strings from stdin", file=sys.stderr)
    database_smiles_list = [line.strip() for line in sys.stdin]

    print("Converting the query SMILES string to a molecule object", file=sys.stderr)
    query_molecule = Chem.MolFromSmiles(query_smiles)
    query_molecule_fps = AllChem.GetMorganFingerprintAsBitVect(query_molecule, radius=2, nBits=1024)

    print("Performing the substructure search and storing matching SMILES strings", file=sys.stderr)
    db_molecules = [Chem.MolFromSmiles(smiles) for smiles in database_smiles_list]
    db_molecules_fps = [AllChem.GetMorganFingerprintAsBitVect(mol, radius=2, nBits=1024) for mol in db_molecules]

    matching_smiles = []

    print("Find similar smiles", file=sys.stderr)
    # Iterate over each pair of molecules
    for i, smiles in enumerate(database_smiles_list):
        # Calculate the Tanimoto Similarity
        similarity = DataStructs.TanimotoSimilarity(db_molecules_fps[i], query_molecule_fps)
        # If the similarity is above threshold_similarity, add the smiles to matching_smiles
        if similarity >= threshold_similarity:
            matching_smiles.append(smiles)

    print(f"Found {len(matching_smiles)} matches", file=sys.stderr)

    print("Printing matching SMILES strings to stdout", file=sys.stderr)
    print('\n'.join(matching_smiles))

except Exception as e:
    print(f"Error: {e}", file=sys.stderr)
