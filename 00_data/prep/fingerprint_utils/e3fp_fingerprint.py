import e3fp
from e3fp.pipeline import confs_from_smiles, fprints_from_mol
from rdkit import Chem

# E3FP fingerprint generator
def e3fp_from_smiles(smiles):
    
    e3fp_params = {'bits': 4096, 'radius_multiplier': 1.5, 'rdkit_invariants': True}
    confgen_params = {'max_energy_diff': 20.0, 'first': 3}
    mol = confs_from_smiles(smiles, "", confgen_params=confgen_params)
    print(mol)
    fprints = fprints_from_mol(mol, fprint_params=e3fp_params)
    
    # Use the first conformer's fingerprint (is this the best approach?)
    return fprints

ACARBOSE_SMILES = "C[C@@H]1[C@H]([C@@H]([C@H]([C@H](O1)O[C@@H]2[C@H](O[C@@H]([C@@H]([C@H]2O)O)O[C@H]([C@@H](CO)O)[C@@H]([C@H](C=O)O)O)CO)O)O)N[C@H]3C=C([C@H]([C@@H]([C@H]3O)O)O)CO"
acarbose_e3fp = e3fp_from_smiles(ACARBOSE_SMILES)
print(acarbose_e3fp)