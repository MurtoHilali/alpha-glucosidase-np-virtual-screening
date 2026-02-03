#!/usr/bin/env python3

import argparse
import pandas as pd
from pathlib import Path
from rdkit import Chem
from rdkit.Chem import AllChem
from concurrent.futures import ProcessPoolExecutor, as_completed


def create_decoy_from_smiles(args):
    smiles, output_path, seed = args

    if output_path.exists():
        return f"[skip] {output_path.name}"

    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return f"[fail] invalid SMILES"

    mol = Chem.AddHs(mol)

    if AllChem.EmbedMolecule(mol, randomSeed=seed) != 0:
        return f"[fail] embed failed: {output_path.name}"

    try:
        AllChem.MMFFOptimizeMolecule(mol)
    except Exception:
        pass

    writer = Chem.SDWriter(str(output_path))
    writer.write(mol)
    writer.close()

    return f"[ok] {output_path.name}"


def load_query_ligands(txt_file: Path) -> pd.DataFrame:
    with open(txt_file) as f:
        ligands = [line.strip() for line in f if line.strip()]

    return pd.DataFrame({
        "Ligand": ligands,
        "Query": [f"Query_{i+1}" for i in range(len(ligands))]
    })


def main(args):
    query_df = load_query_ligands(Path(args.query_txt))
    decoys_df = pd.read_csv(args.decoys_csv)

    if not {"Query", "SMILE"}.issubset(decoys_df.columns):
        raise ValueError("decoys CSV must contain columns: Query, SMILE")

    merged = pd.merge(query_df, decoys_df, on="Query")

    out_dir = Path(args.output_dir)
    out_dir.mkdir(parents=True, exist_ok=True)

    jobs = []
    for ligand, group in merged.groupby("Ligand"):
        for i, row in enumerate(group.itertuples(index=False), start=1):
            out_name = f"{ligand}_Decoy_{i}.sdf"
            out_path = out_dir / out_name
            jobs.append((row.SMILE, out_path, args.seed))

    print(f"spawning {len(jobs)} jobs across {args.n_jobs} cores")

    with ProcessPoolExecutor(max_workers=args.n_jobs) as executor:
        futures = [executor.submit(create_decoy_from_smiles, job) for job in jobs]
        for fut in as_completed(futures):
            print(fut.result())

    print("done. all decoys generated.")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="Parallel decoy SDF generation with RDKit"
    )
    parser.add_argument("--query-txt", required=True)
    parser.add_argument("--decoys-csv", required=True)
    parser.add_argument("--output-dir", default="decoys")
    parser.add_argument("--seed", type=int, default=42)
    parser.add_argument(
        "--n-jobs",
        type=int,
        default=1,
        help="Number of CPU cores to use",
    )

    args = parser.parse_args()
    main(args)
