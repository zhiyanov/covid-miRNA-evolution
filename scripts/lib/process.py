import os
import re
import sys
import tqdm
import json
import pandas as pd

import time
from multiprocessing import Process, Queue

from .utils import *
from .annotate import *

PATH_TO_MIRBASE = "/home/dude/huge/bulk/miRBase/miRBase_22.1.tsv"
MAX_QUEUE_SIZE = 1000


def read(qis, path, threshold=None):
    index = 0
    for covid, rna in read_fasta(path):
        if not threshold is None:
            if percent(rna, "N") > threshold:
                continue

        qis[index].put((covid, rna, index))
        index += 1
        if index >= len(qis):
            index = 0

def clean(qo, path):
    ff = open(path, "w")
    while True:
        try:
            output = qo.get(timeout=2).rstrip("\n") + "\n"
            ff.write(output)
        except:
            break

    ff.close()

def find_seeds(path, qi, miRNAs, protein=None):
    _, wuhan = next(read_fasta(WUHAN_PATH))
    ff = open(path, "w")
    alnpath = "/".join(path.split("/")[:-1])
    while True:
        try:
            covid, rna, index = qi.get(timeout=2)
            align = Align(
                wuhan,
                rna,
                path=f"{alnpath}/.{index}.align.fasta"
            )

            if protein:
                # coordinates on wuhan rna
                s, e = map(int, WUHAN_PROTEINS[protein].split("|")[0].split("-"))
                s -= 1
                e -= 1
                # coordinate on aln rna
                # s, e = align.ref_transform((s, e))
            
            result = f"{covid};"
            for miRNA, seed in miRNAs:
                # coordinates on aln rna
                positions = find_all(seed, rna)
                # coordinates on wuhan rna
                positions = align.aln_transform(positions)
                if protein:
                    positions = [pos for pos in positions if (s <= pos <= e)]
                result += str(positions) + ";"
            result = result.rstrip(";")

            ff.write(f"{result}\n")
        except:
            print(path)
            break
    ff.close()

def find_spaces(path, qi, miRNAs, protein=None):
    _, wuhan = next(read_fasta(WUHAN_PATH))
    ff = open(path, "w")
    alnpath = "/".join(path.split("/")[:-1])
    while True:
        try:
            covid, rna, index = qi.get(timeout=2)
            align = Align(
                wuhan,
                rna,
                path=f"{alnpath}/.{index}.align.fasta"
            )

            if protein:
                # coordinates on wuhan rna
                s, e = map(int, WUHAN_PROTEINS[protein].split("|")[0].split("-"))
                s -= 1
                e -= 1
                # coordinates on aln rna
                s, e = align.ref_transform((s, e))

            result = f"{covid};"
            for miRNA, seed in miRNAs:
                # coordinates on wuhan rna
                positions = find_all(seed, wuhan)
                # coordinates on aln rna
                positions = align.ref_transform(positions)
                if protein:
                    positions = [pos for pos in positions if (s <= pos <= e)]
                positions = [pos for pos in positions if ("N" in rna[pos:pos + 6])]
                # coordinates on wuhan rna
                positions = align.aln_transform(positions)
                result += str(positions) + ";"
            result = result.rstrip(";")
            ff.write(f"{result}\n")
        except:
            print(path)
            break
    ff.close()

def merge(path, name):
    output = open(f"{path}/{name}.csv", "w")
    print("header.csv")
    header = open(f"{path}/header.csv", "r")
    line = next(header, None)
    while line:
        output.write(line.rstrip("\n") + "\n")
        line = next(header, None)
    header.close()
    
    for fp in os.listdir(path):
        if "header" in fp:
            continue
        if "png" in fp:
            continue
        if "result" in fp:
            continue
        if "space" in fp:
            continue
        if "seed" in fp:
            continue

        print(fp)
        ff = open(f"{path}/{fp}", "r")
        line = next(ff, None)
        while line:
            output.write(line.rstrip("\n") + "\n")
            line = next(ff, None)
        ff.close()

def process_seeds(
    gisaid_path,
    threshold,
    output_dir,
    process_num,
    protein=None,
    miRNA_path=None
):
    process_num = int(process_num)
    output_dir = output_dir.rstrip("/")

    if miRNA_path:
        miRNAs = pd.read_csv(miRNA_path, sep=",")
        mimats = set(miRNAs["MIMAT"])
    else:
        mimats = None
    
    ff = open(f"{output_dir}/header.csv", "w")
    
    miRBase = pd.read_csv(PATH_TO_MIRBASE, sep="\t")
    if mimats is None:
        miRBase = miRBase.loc[miRBase["MIMAT"].str.lstrip("MIMAT").astype("int") <= 5000]
    else:
        miRBase = miRBase.loc[miRBase["MIMAT"].isin(mimats)]
    miRNAs = []
    for i, miRNA in miRBase.iterrows():
        # if (mimats) and (miRNA['MIMAT'] not in mimats):
        #     continue
        ff.write(f";{miRNA['MIMAT']}")
        seq = uremove(miRNA["Sequence"])
        start = miRNA["Local start"]
        end = miRNA["Local end"]
        seed = reverse(seq[start:start + 6])
        miRNAs.append((miRNA["MIMAT"], seed))

    ff.write("\n")
    ff.close()
     
    qis = []
    workers = []
    for i in range(process_num):
        # qis.append(Queue(MAX_QUEUE_SIZE))
        qis.append(Queue())
        workers.append(Process(target=find_seeds, args=(
            f"{output_dir}/{i}.csv", qis[-1], miRNAs, protein,
        )))

    reader = Process(target=read, args=(
        qis, gisaid_path, threshold,
    ))

    reader.start()
    for i in range(process_num):
        workers[i].start()

    for i in range(process_num):
        workers[i].join()
    reader.join()

    merge(output_dir, "seed")

def process_spaces(
    gisaid_path,
    threshold,
    output_dir,
    process_num,
    protein=None,
    miRNA_path=None
):
    process_num = int(process_num)
    output_dir = output_dir.rstrip("/")

    if miRNA_path:
        miRNAs = pd.read_csv(miRNA_path, sep=",")
        mimats = set(miRNAs["MIMAT"])
    else:
        mimats = None
    
    ff = open(f"{output_dir}/header.csv", "w")
    
    miRBase = pd.read_csv(PATH_TO_MIRBASE, sep="\t")
    if mimats is None:
        miRBase = miRBase.loc[miRBase["MIMAT"].str.lstrip("MIMAT").astype("int") <= 5000]
    else:
        miRBase = miRBase.loc[miRBase["MIMAT"].isin(mimats)]

    miRNAs = []
    for i, miRNA in miRBase.iterrows():
        # if (mimats) and (miRNA['MIMAT'] not in mimats):
        #     continue
        ff.write(f";{miRNA['MIMAT']}")
        seq = uremove(miRNA["Sequence"])
        start = miRNA["Local start"]
        end = miRNA["Local end"]
        seed = reverse(seq[start:start + 6])
        miRNAs.append((miRNA["MIMAT"], seed))

    ff.write("\n")
    ff.close()
     
    qis = []
    workers = []
    for i in range(process_num):
        # qis.append(Queue(MAX_QUEUE_SIZE))
        qis.append(Queue())
        workers.append(Process(target=find_spaces, args=(
            f"{output_dir}/{i}.csv", qis[-1], miRNAs, protein,
        )))

    reader = Process(target=read, args=(
        qis, gisaid_path, threshold,
    ))

    reader.start()
    for i in range(process_num):
        workers[i].start()

    for i in range(process_num):
        workers[i].join()
    reader.join()

    merge(output_dir, "space")
