from rdkit import Chem
from rdkit.Chem import rdFingerprintGenerator, MACCSkeys

# RDKIT Morgan fingerprint generator, ECFP-like
ecfp_gen = rdFingerprintGenerator.GetMorganGenerator(
    radius=4, 
    fpSize=4096, 
    includeChirality=True
)

def ecfp_from_smiles(smiles):
    '''Generate ECFP (Morgan) fingerprint bitstring from a SMILES string.
    Args: 
        smiles (str): SMILES representation of the molecule.
    Returns:
        str or None: Fingerprint bitstring if molecule is valid, else None.
    '''
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return None
    fp = ecfp_gen.GetFingerprint(mol)
    return fp.ToBitString()

# RDKIT Morgan fingerprint generator, with pharmoacophore features
feat_gen = rdFingerprintGenerator.GetMorganGenerator(
    radius=4,
    fpSize=4096, 
    atomInvariantsGenerator=rdFingerprintGenerator.GetMorganFeatureAtomInvGen()
)

def pharmacophore_fp_from_smiles(smiles):
    '''Generate pharmacophore fingerprint bitstring from a SMILES string.
    Args: 
        smiles (str): SMILES representation of the molecule.
    Returns:
        str or None: Fingerprint bitstring if molecule is valid, else None.
    '''
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return None
    fp = feat_gen.GetFingerprint(mol)
    return fp.ToBitString()

# RDKIT MACCS keys fingerprint generator
def maccs_fp_from_smiles(smiles):
    '''Generate MACCS keys fingerprint bitstring from a SMILES string.
    Args: 
        smiles (str): SMILES representation of the molecule.
    Returns:
        str or None: Fingerprint bitstring if molecule is valid, else None.
    '''
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return None
    fp = MACCSkeys.GenMACCSKeys(mol)
    return fp.ToBitString()

ACARBOSE_SMILES = "C[C@@H]1[C@H]([C@@H]([C@H]([C@H](O1)O[C@@H]2[C@H](O[C@@H]([C@@H]([C@H]2O)O)O[C@H]([C@@H](CO)O)[C@@H]([C@H](C=O)O)O)CO)O)O)N[C@H]3C=C([C@H]([C@@H]([C@H]3O)O)O)CO"
acarbose_ecfp = ecfp_from_smiles(ACARBOSE_SMILES)
acarbose_pharmacophore_fp = pharmacophore_fp_from_smiles(ACARBOSE_SMILES)
acarbose_maccs_fp = maccs_fp_from_smiles(ACARBOSE_SMILES)