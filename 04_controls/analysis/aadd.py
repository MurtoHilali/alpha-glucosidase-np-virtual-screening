import pandas as pd

def load_and_process_gnina(parquet: str | pd.DataFrame) -> pd.DataFrame:
    if type(parquet) == str:
        df = pd.read_parquet(parquet)

    df["cnn_vs"] = df["cnn_pose_score"] * df["cnn_affinity"]

    df = (
        df
        .sort_values("cnn_vs", ascending=False)
        .drop_duplicates(subset="identifier")
    )

    return df