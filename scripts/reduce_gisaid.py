import os
import re
import sys
import tqdm
import json
import pandas as pd

from lib.utils import read_fasta


def condition(description):
    description = description.lower()
    # if not ("russia" in description):
    #     return False
    # if not ("mos" in description):
    #     return False
    # if not ("england" in description):
    #     return False
    # if not ("lond" in description):
    #     return False
    # if not ("germany" in description):
    #     return False
    # if not ("/be-" in description):
    #     return False
    if not ("taiwan" in description):
        return False
    return True

if __name__ == '__main__':
    gisaid_path = sys.argv[1]

    for description, rna in tqdm.tqdm(read_fasta(gisaid_path), total=12962156):
        if condition(description):
            print(description)
            print(rna)
