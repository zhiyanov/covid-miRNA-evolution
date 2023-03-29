import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import tqdm
import sys

from scipy.stats import linregress
from scipy.stats import mannwhitneyu
from scipy.stats import ks_2samp
from scipy.stats import rankdata

PATH_TO_MIRBASE = "/home/dude/huge/bulk/miRBase/miRBase_22.1.tsv"


def random_mirnas(count, mirna_df):
    miRBase = pd.read_csv(PATH_TO_MIRBASE, sep="\t")
    miRBase = miRBase[["miRNA", "MIMAT"]]
    miRBase = miRBase.loc[
        ~miRBase["MIMAT"].isin(mirna_df["MIMAT"])
    ].drop_duplicates()

    index = len(mirna_df)
    for i in range(count):
        result_df = pd.DataFrame(miRBase.iloc[
            np.random.permutation(
                len(miRBase)
            )[:index]
        ])

        result_df["CPM"] = mirna_df["CPM"].to_list()

        yield result_df

def distribution(count_df, meta_df, mirna_df):
    meta_df.loc[:, "group"] = "first"
    meta_df.loc[
        meta_df["scorpio_call"].str.contains("micron"),
        "group"
    ] = "second"

    distribution = []
    iterator = enumerate(random_mirnas(1001, mirna_df))
    for i, random_df in iterator:
        if i == 1000:
            random_df = mirna_df

        random_df.index = random_df["MIMAT"]
        
        mimats = random_df["MIMAT"].to_list()
        total_expression = random_df["CPM"].sum()

        count_matrix = np.array(count_df[mimats])
        ones_vector = np.ones((len(mimats), 1))
        expression_vector = np.array(
            random_df.loc[mimats]["CPM"]
        ).reshape((-1, 1)) / total_expression

        count_vector = np.matmul(
            count_matrix, ones_vector)
        weighted_count_vector = np.matmul(
            count_matrix, expression_vector)

        df = pd.DataFrame()
        df.index = count_df.index
        df.loc[:, "count"] = count_vector
        df.loc[:, "weighted count"] = weighted_count_vector

        df = df.join(
            meta_df["group"], how="left"
        )

        try:
            wstat, wpv = mannwhitneyu(
                df.loc[df["group"] == "first"]["weighted count"],
                df.loc[df["group"] == "second"]["weighted count"]
            )
        except:
            wstat, wpv = 0, 1
        fwm = df.loc[df["group"] == "first"]["weighted count"].mean()
        swm = df.loc[df["group"] == "second"]["weighted count"].mean()

        try:
            stat, pv = mannwhitneyu(
                df.loc[df["group"] == "first"]["count"],
                df.loc[df["group"] == "second"]["count"]
            )
        except:
            stat, pv = 0, 1
        fm = df.loc[df["group"] == "first"]["count"].mean()
        sm = df.loc[df["group"] == "second"]["count"].mean()

        distribution.append([fwm, swm, wstat, wpv, fm, sm, stat, pv])

    result_df = pd.DataFrame(
        columns=[
            "first weighted mean", "second weighted mean",
            "weigted stat", "weighted pvalue",
            "first mean", "second mean",
            "stat", "pvalue",
        ],
        data=distribution
    )

    return result_df

if __name__ == "__main__":
    gisaid = sys.argv[1]
    region = sys.argv[2]
    project = sys.argv[3]

    count_df = pd.read_csv(
        f"./data/{gisaid}_{region}_count.csv",
        sep=",",
        index_col=0
    )
    
    meta_df = pd.read_csv(
        f"./data/{gisaid}_meta.csv",
        sep=",",
        index_col=0
    )

    mirna_df = pd.read_csv(
        f"../../miRNA/expressed_{project}.csv",
        sep=","
    )
    
    distribution_df = distribution(
        count_df,
        meta_df,
        mirna_df
    )
    
    print(gisaid, region, project)
    distribution_df.to_csv(
        f"./data/{gisaid}_{region}_{project}_distribution.csv",
        index=None,
        sep=","
    )
