import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
import tqdm

from scipy.stats import wilcoxon
from scipy.stats import ttest_ind
from scipy.stats import ttest_rel
from scipy.stats import normaltest
from scipy.stats import mannwhitneyu
from scipy.stats import permutation_test

COVIDS = [
    ">hCoV-19/Russia/MOS-CRIE-7765209290/2022|2022-07-28|2022-08-11",
    ">hCoV-19/Russia/MOS-CRIE-7761820132/2022|2022-07-10|2022-08-11",
    ">hCoV-19/Russia/MOS-CRIE-7762010932/2022|2022-07-11|2022-08-11",
    ">hCoV-19/Russia/MOS-CRIE-7764507299/2022|2022-07-24|2022-08-11",
    ">hCoV-19/Russia/MOS-CRIE-7767361454/2022|2022-08-07|2022-08-30",
    ">hCoV-19/Russia/MOS-CRIE-L106S0285u/2021|2021-10-25|2022-01-10",
    ">hCoV-19/Russia/MOS-CRIE-7762215713/2022|2022-07-12|2022-08-11",
    ">hCoV-19/Russia/MOS-CRIE-7762041931/2022|2022-07-11|2022-08-11",
    ">hCoV-19/Russia/MOS-CRIE-7767400291/2022|2022-08-07|2022-08-30",
    ">hCoV-19/Russia/MOS-CRIE-4709666645/2022|2022-07-24|2022-08-11",
    ">hCoV-19/Russia/MOS-CRIE-7764749034/2022|2022-07-25|2022-08-11",
    ">hCoV-19/Russia/MOS-CRIE-7761979971/2022|2022-07-11|2022-08-11",
    ">hCoV-19/Russia/MOS-CRIE-7762673598/2022|2022-07-14|2022-08-11",
    ">hCoV-19/Russia/MOS-CRIE-4770225711/2022|2022-08-07|2022-08-30",
    ">hCoV-19/Russia/MOS-CRIE-L189S1335/2021|2021-05-25|2021-10-21",
    ">hCoV-19/Russia/MOS-CRIE-7765790035/2022|2022-07-30|2022-08-11",
    ">hCoV-19/Russia/MOS-CRIE-L191V0519_4656/2021|2021-10-29|2022-01-10",
    ">hCoV-19/Russia/MOS-CRIE-L189P3398/2021|2021-10-24|2021-11-30",
    ">hCoV-19/Russia/MOS-CRIE-L106B1105b/2021|2021-04-01|2021-05-12",
    ">hCoV-19/Russia/MOS-CRIE-L189T1645u_5110/2021|2021-11-27|2022-03-24",
    ">hCoV-19/Russia/MOS-CRIE-L189L3945_6653/2021|2021-11-18|2022-03-24",
    ">hCoV-19/Russia/MOS-SAB-1502/2021|2021-02-15|2021-11-19"
]

PROTEINS = [
    "NSP1", "NSP2", "NSP3", "NSP4",
    "NSP5", "NSP6", "NSP7", "NSP8",
    "NSP9", "NSP10", "NSP11", "NSP12",
    "NSP13", "NSP14", "NSP15", "NSP16",
    "Spike", "NS8", ""
]


def statistic(first, second):
    return (first - second).sum()

if __name__ == '__main__':
    coad_miRNA_df = pd.read_csv("../miRNA/COAD_LUAD_DESeq.csv", sep=",")
    luad_miRNA_df = pd.read_csv("../miRNA/LUAD_COAD_DESeq.csv", sep=",")
    for protein in PROTEINS:
        if protein:
            seed_df = pd.read_csv(f"../results/russia_mos_{protein}/result.csv", sep=";", index_col=0)
        else:
            seed_df = pd.read_csv(f"../results/russia_mos/result.csv", sep=";", index_col=0)
        seed_df = seed_df.loc[COVIDS]

        seed_df["date"] = seed_df.index.str.split("|").str[1]
        seed_df["date"] = pd.to_datetime(seed_df["date"], format="%Y-%m-%d")
        seed_df = seed_df.loc[seed_df["date"].astype("str").apply(lambda dt: len(dt.split("-"))) == 3]
        seed_df = seed_df.sort_values(by=["date"])

        seed_df["year"] = seed_df["date"].astype("str").str.split("-").str[0].astype("int")
        seed_df["month"] = seed_df["date"].astype("str").str.split("-").str[1].astype("int")
        seed_df["day"] = seed_df["date"].astype("str").str.split("-").str[2].astype("int")
        seed_df["days"] = seed_df["day"] + seed_df["month"] * 30 + seed_df["year"] * 365
        start = int(seed_df.iloc[0]["days"])
        seed_df["days"] -= start
        seed_df = seed_df.sort_values(by=["days"])
        seed_df = seed_df.drop(columns=["year", "month", "day", "date"])
        seed_df = seed_df.loc[seed_df["days"] > 400]

        seed_df_miRNAs = set(seed_df.columns.to_list()) - set(["date"])

        miRNAs = set(coad_miRNA_df["MIMAT"].to_list()) & seed_df_miRNAs
        coad_miRNA_df.index = coad_miRNA_df["MIMAT"]
        total_exp = coad_miRNA_df["COAD CPM"].sum()
        seed_df["count"] = [0] * len(seed_df)
        for miRNA in miRNAs:
            fraction = 0.
            fraction = (coad_miRNA_df.loc[miRNA, "COAD CPM"] / total_exp)
            seed_df["count"] += seed_df[miRNA].astype("str").apply(lambda ll: len(ll.split(","))) * fraction
        first_sample = np.array(seed_df["count"])
        
        miRNAs = set(luad_miRNA_df["MIMAT"].to_list()) & seed_df_miRNAs
        luad_miRNA_df.index = luad_miRNA_df["MIMAT"]
        total_exp = luad_miRNA_df["LUAD CPM"].sum()
        seed_df["count"] = [0] * len(seed_df)
        for miRNA in miRNAs:
            fraction = 0.
            fraction = (luad_miRNA_df.loc[miRNA, "LUAD CPM"] / total_exp)
            seed_df["count"] += seed_df[miRNA].astype("str").apply(lambda ll: len(ll.split(","))) * fraction
        second_sample = np.array(seed_df["count"])
        
        res = permutation_test(
            (first_sample, second_sample),
            statistic,
            permutation_type="samples"
        )
        stat, pvalue = res.statistic, res.pvalue

        
        if not protein:
            protein = "RNA"

        print(
            protein,
            stat, pvalue,
            sum(first_sample) / len(first_sample),
            sum(second_sample) / len(second_sample)
        )
