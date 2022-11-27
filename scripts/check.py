import os
import pandas as pd
import tqdm

PATH = "../results/russia_mos_LUAD"

ANNOTATION = {
    "5P": "1-265",
    "NSP1": "266-805",
    "NSP2": "806-2719",
    "NSP3": "2720-8554",
    "NSP4": "8555-10054",
    "NSP5": "10055-10972",
    "NSP6": "10973-11842",
    "NSP7": "11843-12091",
    "NSP8": "12092-12685",
    "NSP9": "12686-13024",
    "NSP10": "13025-13441",
    "NSP11": "13442-13480",
    # "NSP12": "13442-13468",
    "NSP12": "13468-16236",
    "NSP13": "16237-18039",
    "NSP14": "18040-19620",
    "NSP15": "19621-20658",
    "NSP16": "20659-21552",
    "Spike": "21563-25384",
    "NS3": "25393-26220",
    "E": "26245-26472",
    "M": "26523-27191",
    "NS6": "27202-27387",
    "NS7a": "27394-27759",
    "NS7b": "27756-27887",
    "NS8": "27894-28259",
    "N": "28274-29533",
    "NS9b": "28284-28577",
    "NS9c": "28734-28955",
    "3P": "29534-30331"
}


def region_count(reg, start=None, end=None):
    reg = reg.lstrip("[").rstrip("]").split(", ")
    if not (start is None) and not (end is None):
        reg = [s for s in reg if s and (start <= int(s) <= end)]
    else:
        reg = [s for s in reg if s]
    return len(reg)

rna_df = pd.read_csv(f"{PATH}/seed.csv", sep=";", index_col=0)
rna_df["rna proteins"] = 0
for column in tqdm.tqdm(rna_df.drop(columns=["rna proteins"]).columns):
    for protein, region in ANNOTATION.items():
        start, end = map(lambda x: int(x) - 1, region.split("-"))
        rna_df["rna proteins"] += rna_df[column].apply(
            lambda x: region_count(x, start, end)
        )

    rna_df[column] = rna_df[column].apply(lambda x: region_count(x))
rna_df["rna"] = rna_df.drop(columns=["rna proteins"]).sum(axis=1)
rna_df = rna_df[["rna", "rna proteins"]]


flag = False
protein_df = None
for protein, region in tqdm.tqdm(ANNOTATION.items()):
    df = pd.read_csv(f"{PATH}_{protein}/seed.csv", sep=";", index_col=0)
    for column in df.columns:
        df[column] = df[column].apply(lambda x: region_count(x))
    df = df.sum(axis=1).to_frame()
    df = df.rename(columns={0: "sum"})

    if not flag:
        protein_df = df
        protein_df = protein_df.rename(columns={"sum": "protein"})
        flag = True
    else:
        protein_df = protein_df.join(df, how="inner")
        protein_df = protein_df.sum(axis=1).to_frame()
        protein_df = protein_df.rename(columns={0: "protein"})

result_df = rna_df.join(protein_df, how="inner")
