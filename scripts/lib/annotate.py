import sys
import os

from Bio.Seq import Seq
from Bio.Align.Applications import MafftCommandline

from .align import *
from .utils import read_fasta

WUHAN_PATH = "/home/dude/huge/dude/long-covid/gisaid/wuhan.fasta"
WUHAN_PROTEINS = {
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

def annotate(seq, path=None):
    _, wuhan = next(read_fasta(WUHAN_PATH))
    ref, aln = align(wuhan, seq, path=path)

    result = {}
    frw = forward(ref) 
    bcw = backward(aln)
    for protein, reg in WUHAN_PROTEINS.items():
        for reg in reg.split("|"):
            s, e = reg.split("-")
            s, e = frw[int(s) - 1], frw[int(e) - 1]
            s, e = bcw[s] + 1, bcw[e] + 1
            
            if protein in result:
                result[protein] += f"|{s}-{e}"
            else:
                result[protein] = f"{s}-{e}"

    return result
