from rdkit import Chem
from rdkit.Chem import rdMolDescriptors, AllChem
import numpy as np

def usr_from_smiles(smiles):
    '''Generate USR fingerprint vector from a SMILES string.
    Args:
        smiles (str): SMILES representation of the molecule.
    Returns:
        list or None: USR fingerprint vector if molecule is valid and computed within the timeout, else None.
    '''
    # first, generate + optimize conformers from SMILES

    mol = Chem.MolFromSmiles(smiles)
    if mol is None or mol.GetNumAtoms() <= 2:
        return None
    mol = Chem.AddHs(mol)
    num_conformers = 5
    cids = AllChem.EmbedMultipleConfs(mol, numConfs=num_conformers, params=AllChem.ETKDGv3())
    for cid in cids:
        AllChem.MMFFOptimizeMolecule(mol, confId=cid)
    # generate USR vectors for conformers
    print("Getting USR vectors...")
    usr_vecs = [rdMolDescriptors.GetUSR(mol, confId=conf) for conf in cids]
    usr_array = np.array(usr_vecs)
    
    # take unified mean of vector array
    unified_usr = usr_array.mean(axis=0)
    print("Got unified USR vector.")
    try:
        return unified_usr.tolist()
    except Exception as e:
        print(f"Error in usr_from_smiles: {e}")
        return None
