from rdkit import Chem
from rdkit.Chem import AllChem
import pandas as pd
import os

failures = []

def create_ligand_from_smiles(smiles: str, ligand_name: str):
    '''
    Create a ligand SDF file of multiple conformers from a SMILES string.
    
    Parameters:
        - smiles (str): The SMILES representation of the ligand.
        - ligand_name (str): The name of the ligand file to be created.
    
    Outputs:
        - An SDF file named '<ligand_name>.sdf' containing the 3D structure of the ligand.
    '''
    sdf_filename = f"/home/mhilali/projects/def-bingalls/mhilali/natural_agis/ligands/{ligand_name}.sdf"
    if os.path.exists(sdf_filename):
        print(f'{ligand_name}.sdf exists already.')
        return 1
    
    # Convert SMILES to RDKit molecule
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        raise ValueError("Invalid SMILES string provided.")
    mol = Chem.AddHs(mol)
    AllChem.EmbedMolecule(mol, randomSeed=42)
    
    # Optimize molecular force field
    try:
        AllChem.MMFFOptimizeMolecule(mol)
    except Exception as e:
        print(f'Could not optimize {ligand_name}.')
        failures.append(ligand_name)
        return 1
      
    # Write out SDF
    try:
        writer = Chem.SDWriter(sdf_filename)
        writer.write(mol)
        writer.close()
        print(f"Ligand SDF file '{sdf_filename}' created successfully.")
    except Exception as e:
        print(f'Could not not create {ligand_name}. Error: {e}.')
    

high_smililarity_candidates = pd.read_csv('high_similarity_candidates.csv')
ligands = high_smililarity_candidates[['identifier', 'canonical_smiles']].values.tolist()

for identifier, smiles in ligands:
    create_ligand_from_smiles(smiles, identifier)

with open("fails.txt", "w") as f:
    for item in failures:
        f.write(f'{item}\n')