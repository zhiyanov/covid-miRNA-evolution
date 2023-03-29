import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import tqdm
import sys


if __name__ == "__main__":
    gisaid = sys.argv[1]
    region = sys.argv[2]
    project = sys.argv[3]

    df = pd.read_csv(
        f"./data/{gisaid}_{region}_{project}_distribution.csv",
        sep=","
    )

    condition = \
        (df["weigted stat"] > df.iloc[-1]["weigted stat"]) & \
        (df["stat"] > df.iloc[-1]["stat"])

    print(gisaid, region, project, len(df.loc[condition]) / len(df))


