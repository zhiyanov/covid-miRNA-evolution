import numpy as np
import pandas as pd
import tqdm

from scipy.stats import mannwhitneyu

PATH = "/home/dude/huge/dude/long-covid"

PROTEINS = [
    "RNA",
    # "5P",
    # "NSP1",
    # "NSP2",
    # "NSP3",
    # "NSP4",
    # "NSP5",
    # "NSP6",
    # "NSP7",
    # "NSP8",
    # "NSP9",
    # "NSP10",
    # "NSP11",
    # "NSP12",
    # "NSP13",
    # "NSP14",
    # "NSP15",
    # "NSP16",
    # "Spike",
    # "NS3",
    # "E",
    # "M",
    # "NS6",
    # "NS7a",
    # "NS7b",
    # "NS8",
    # "N",
    # "NS9b",
    # "NS9c",
    # "3P"
]

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

def region_count(reg, start=None, end=None):
    reg = reg.lstrip("[").rstrip("]").split(", ")
    if not (start is None) and not (end is None):
        reg = [s for s in reg if s and (start <= int(s) <= end)]
    else:
        reg = [s for s in reg if s]
    return len(reg)


# meta
meta_df = pd.read_csv(
    f"{PATH}/gisaid/russia_mos_meta.csv",
    sep=",", index_col=0)

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

meta_df.loc[:, "group"] = "first"
meta_df.loc[meta_df["days"] > 700, "group"] = "second"

# mirna
miRNA_df = pd.read_csv(f"{PATH}/miRNA/random_LUAD.csv", sep=",")
miRNA_df.index = miRNA_df["MIMAT"]

# TOTAL_EXPRESSION = 1359187
TOTAL_EXPRESSION = miRNA_df["CPM"].sum()
MIMATS = miRNA_df["MIMAT"].to_list()

# count
count_df = None
for protein in PROTEINS:
    if protein != "RNA":
        df = pd.read_csv(f"{PATH}/results/russia_mos_vocs_LUAD_{protein}/seed.csv", sep=";", index_col=0)
    else:
        df = pd.read_csv(f"{PATH}/results/russia_mos_vocs_LUAD/seed.csv", sep=";", index_col=0)

    df.index = df.index.str.lstrip(">")
    df = df[MIMATS]

    df["strain"] = df.index
    df["region"] = [protein] * len(df)

    if count_df is None:
        count_df = df
    else:
        count_df = pd.concat([count_df, df], ignore_index=True)

for mimat in MIMATS:
    count_df[mimat] = count_df[mimat].apply(lambda x: region_count(x))

count_df.loc[:, "count"] = 0
count_df.loc[:, "weighted count"] = 0
# for mimat in tqdm.tqdm(MIMATS):
for mimat in MIMATS:
    expression = (miRNA_df.loc[mimat, "CPM"] / TOTAL_EXPRESSION)
    count_df["count"] += count_df[mimat]
    count_df["weighted count"] += count_df[mimat] * expression

# analisis
count_df = count_df.loc[count_df["region"] == "RNA"]
count_df.index = count_df["strain"]
count_df = count_df.drop(columns=["strain"])

count_df = count_df.join(meta_df[["group"]], how="left")
first_mean = count_df.loc[
    count_df["group"] == "first",
    "weighted count"
].mean()
second_mean = count_df.loc[
    count_df["group"] == "second",
    "weighted count"
].mean()

try:
    _, pvalue = mannwhitneyu(
        count_df.loc[
            count_df["group"] == "first",
            "weighted count"
        ],
        count_df.loc[
            count_df["group"] == "second",
            "weighted count"
        ]
    )
except:
    pvalue = 1

print(first_mean, second_mean, pvalue, end=" ")

first_mean = count_df.loc[
    count_df["group"] == "first",
    "count"
].mean()
second_mean = count_df.loc[
    count_df["group"] == "second",
    "count"
].mean()

try:
    _, pvalue = mannwhitneyu(
        count_df.loc[
            count_df["group"] == "first",
            "count"
        ],
        count_df.loc[
            count_df["group"] == "second",
            "count"
        ]
    )
except:
    pvalue = 1
print(first_mean, second_mean, pvalue)
