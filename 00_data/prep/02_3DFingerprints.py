#!/usr/bin/env python3
import os
import time
import glob
import numpy as np
import pandas as pd

from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem.rdMolDescriptors import GetUSRCAT
from rdkit import RDLogger
RDLogger.DisableLog('rdApp.*')  # disable rdkit warnings

from joblib import Parallel, delayed
from tqdm.auto import tqdm

# --- parquet writer backend ---
PARQUET_ENGINE = "pyarrow"


# ----------------------------
# paths / config
# ----------------------------
IN_CSV   = "/path/to/coconut_cleaned.csv"
OUT_DIR  = "/path/to/coconut_with_usrcat_parts"
CKPT     = "/path/to/coconut_with_usrcat.checkpoint.txt"

USECOLS  = ["identifier", "canonical_smiles", "normalized_smiles"]

N_USR    = 60
USR_COLS = [f"usr_{i}" for i in range(N_USR)]

# streaming + parallel knobs
READ_CHUNKSIZE = 10000
N_JOBS         = int(os.environ.get("SLURM_CPUS_PER_TASK", "8"))
JOBLIB_BACKEND = "multiprocessing"


# ----------------------------
# checkpoint helpers
# ----------------------------
def load_checkpoint() -> int:
    if not os.path.exists(CKPT):
        return 0
    with open(CKPT, "r") as f:
        s = f.read().strip()
    return int(s) if s else 0

def save_checkpoint(row_offset: int):
    tmp = CKPT + ".tmp"
    with open(tmp, "w") as f:
        f.write(str(row_offset))
    os.replace(tmp, CKPT)

def part_path(part_idx: int) -> str:
    return os.path.join(OUT_DIR, f"part_{part_idx:06d}.parquet")


# ----------------------------
# core compute
# ----------------------------
def compute_usrcat_row(smiles: str):
    fail = [np.nan] * N_USR

    if smiles is None:
        return fail
    smiles = str(smiles).strip()
    if not smiles:
        return fail

    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return fail

    try:
        mol = Chem.AddHs(mol)
        status = AllChem.EmbedMolecule(
            mol,
            useExpTorsionAnglePrefs=True,
            useBasicKnowledge=True,
            maxAttempts=20
        )
    except Exception as e:
        print(time.strftime("%Y-%m-%d %H:%M:%S"),
              f"Embed runtime error for {smiles[:10]}..., {e}")
        return fail

    if status != 0 or mol.GetNumConformers() == 0:
        return fail

    try:
        usr = GetUSRCAT(mol)
        return [float(x) for x in usr]
    except Exception:
        return fail


# ----------------------------
# main streaming pipeline
# ----------------------------
def process():
    print("Making output directory...")
    os.makedirs(OUT_DIR, exist_ok=True)
    print("Current time:", time.strftime("%Y-%m-%d %H:%M:%S"))

    start_row = load_checkpoint()
    if start_row > 0:
        print(f"Resuming from row offset: {start_row}")

    # figure out what part index we should be writing next
    start_part_idx = start_row // READ_CHUNKSIZE

    print("Starting CSV read...")
    reader = pd.read_csv(
        IN_CSV,
        usecols=USECOLS,
        chunksize=READ_CHUNKSIZE,
        skiprows=range(1, start_row + 1) if start_row > 0 else None
    )

    processed = start_row
    part_idx = start_part_idx

    for df in reader:
        t0 = time.time()

        # if rerunning and the part already exists, skip it (extra safety)
        out_path = part_path(part_idx)
        if os.path.exists(out_path):
            processed += len(df)
            save_checkpoint(processed)
            print(f"part {part_idx} already exists, skipping; total={processed}")
            part_idx += 1
            continue

        smiles_list = df["canonical_smiles"].astype("string").fillna("").tolist()

        print(f"Processing part {part_idx}, rows={len(df)}...")
        usr_rows = Parallel(
            n_jobs=N_JOBS,
            batch_size="auto",
            backend=JOBLIB_BACKEND,
            verbose=0
        )(
            delayed(compute_usrcat_row)(s)
            for s in tqdm(smiles_list, desc=f"part {part_idx}", unit="mol", leave=False)
        )
        print("Parallel computation done. len=", len(usr_rows), " time=", time.strftime("%H:%M:%S"))
        print("Combining results and writing parquet...")
        
        t1 = time.time()
        usr_df = pd.DataFrame(usr_rows, columns=USR_COLS)
        print(f"Computed USRCAT for {len(usr_df)} molecules in {time.time() - t1:.1f}s")
        
        t2 = time.time()
        out_df = pd.concat([df.reset_index(drop=True), usr_df], axis=1)
        out_df.to_parquet(out_path, engine=PARQUET_ENGINE, index=False)
        print(f"Wrote {out_path}")
        print(f"Combined DataFrame in {time.time() - t2:.1f}s")
        
        processed += len(df)
        save_checkpoint(processed)

        part_idx += 1

    print("USRCAT Descriptor Computation Complete!")
    print(f"Parquet parts in: {OUT_DIR}")

    parts = sorted(glob.glob(f"{OUT_DIR}/part_*.parquet"))
    df = pd.concat([pd.read_parquet(p) for p in parts], ignore_index=True)
    df.to_parquet("/path/to/coconut_with_usrcat.parquet", index=False)
    df.to_csv("/path/to/coconut_with_usrcat.csv", index=False)

if __name__ == "__main__":
    process()