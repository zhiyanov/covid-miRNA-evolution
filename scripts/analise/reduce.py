import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import sys


def read_fasta(path):
    result = ""
    description = ""

    ff = open(path, "r")

    line = next(ff, None)
    while line:
        line = line.rstrip("\n")
        if not line:
            line = next(ff, None)
            continue

        if line.startswith(">"):
            if result:
                yield description, result
            description = line.lstrip(">")
            result = ""
        else:
            result += line

        line = next(ff, None)

    yield description, result
    ff.close()
    
def percent(seq, chr="N"):
    count = 0
    for char in seq:
        if char == chr:
            count += 1
    return count / len(seq)

def meta_data(meta_path, rna_path):
    meta_df = pd.read_csv(
        meta_path,
        sep=",", index_col=0)

    rna_df = pd.DataFrame(columns=["strain", "length", "polyA", "percN"])
    for covid, rna in read_fasta(rna_path):
        length = len(rna)
        perc = percent(rna, "N")
        poly = len(rna) - len(rna.rstrip("A"))
        append_df = pd.DataFrame(
            columns=rna_df.columns.to_list(),
            data=[[
                covid,
                length,
                poly,
                perc
            ]]
        )
        rna_df = pd.concat([rna_df, append_df], ignore_index=True)

    rna_df.index = rna_df["strain"]
    rna_df = rna_df.drop(columns=["strain"])
    meta_df = meta_df.join(rna_df, how="inner")

    meta_df["date"] = meta_df.index.str.split("|").str[1]
    meta_df["date"] = pd.to_datetime(meta_df["date"], format="%Y-%m-%d")
    meta_df = meta_df.sort_values(by=["date"])

    meta_df["date"] = meta_df.index.str.split("|").str[1]
    meta_df["date"] = pd.to_datetime(meta_df["date"], format="%Y-%m-%d")
    meta_df = meta_df.sort_values(by=["date"])

    meta_df = meta_df.loc[meta_df["date"].astype("str").apply(lambda dt: len(dt.split("-"))) == 3]
    meta_df["year"] = meta_df["date"].astype("str").str.split("-").str[0].astype("int")
    meta_df["month"] = meta_df["date"].astype("str").str.split("-").str[1].astype("int")
    meta_df["day"] = meta_df["date"].astype("str").str.split("-").str[2].astype("int")
    meta_df["days"] = meta_df["day"] + meta_df["month"] * 30 + meta_df["year"] * 365
    start = int(meta_df.iloc[0]["days"])
    meta_df["days"] -= start
    meta_df = meta_df.drop(columns=["year", "month", "day"])

    return meta_df

def reduce_data(meta_path, rna_path, seed_path):
    meta_df = meta_data(meta_path, rna_path)
    meta_df = meta_df.loc[
        (~meta_df["scorpio_call"].isna()) & \
        (meta_df["percN"] == 0.) & \
        (meta_df["length"] > 29000)
    ] 
    
    seed_df = pd.read_csv(seed_path, sep=";", index_col=0)
    seed_df.index = seed_df.index.str.lstrip(">")
    seed_df = seed_df.loc[meta_df.index]

    return seed_df, meta_df

if __name__ == "__main__":
    gisaid = sys.argv[1]

    meta_path = f"../../gisaid/{gisaid}_meta.csv"
    rna_path = f"../../gisaid/{gisaid}.fasta"
    seed_path = f"../../results/{gisaid}/seed.csv"

    seed_df, meta_df = reduce_data(
        meta_path, rna_path,
        seed_path
    )

    print(gisaid)
    seed_df.to_csv(
        f"./data/{gisaid}_seed.csv",
        sep=","
    )
    meta_df.to_csv(
        f"./data/{gisaid}_meta.csv",
        sep=","
    )
