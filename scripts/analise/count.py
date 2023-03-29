import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import sys

from scipy.stats import linregress
from scipy.stats import mannwhitneyu
from scipy.stats import ks_2samp
from scipy.stats import rankdata

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

def count_data(seed_df, region):
    if region != "RNA":
        s, e = map(lambda x: int(x) - 1, ANNOTATION[region].split("-"))
    else:
        s, e = None, None

    for mimat in seed_df.columns:
        seed_df[mimat] = seed_df[mimat].apply(lambda x: region_count(x, s, e))

    return seed_df


if __name__ == "__main__":
    gisaid = sys.argv[1]
    region = sys.argv[2]

    seed_df = pd.read_csv(
        f"./data/{gisaid}_seed.csv",
        index_col=0,
        sep=","
    )

    count_df = count_data(seed_df, region)
    
    print(gisaid, region)
    count_df.to_csv(
        f"./data/{gisaid}_{region}_count.csv",
        sep=","
    )
