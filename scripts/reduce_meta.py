import os
import re
import sys
import tqdm
import json
import pandas as pd

from utils import read_fasta


def condition(description):
    if not ("russia" in description.lower()):
        return False
    if not ("mos" in description.lower()):
        return False

    return True

if __name__ == '__main__':
    gisaid_path = sys.argv[1]
    meta_path = sys.argv[2]

    ids = []
    for description, rna in read_fasta(gisaid_path):
        ids.append(description.split("|")[0].lstrip(">"))

    ids = set(ids)
    
    ff = open(meta_path, "r")
    print(next(ff), end="")
    for line in ff:
        meta = line.split("\t")
        id, stam = meta[0], meta[11]
        if id in ids:
            print(line, end="")
    ff.close()

