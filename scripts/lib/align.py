import sys
import os

from Bio.Seq import Seq
from Bio.Align.Applications import MafftCommandline

from multiprocessing import Process, Queue

from .utils import read_fasta

WUHAN_PATH = "/home/dude/huge/dude/long-covid/gisaid/wuhan.fasta"

class Align:
    def __init__(self, ref, aln, path=None):
        self.ref = ref
        self.aln = aln

        _ref, _aln = align(
            ref, aln,
            path=path,
            save=True
        )
        self._ref = _ref
        self._aln = _aln

        self.ref_frw = forward(self._ref)
        self.ref_bcw = backward(self._ref)
        self.aln_frw = forward(self._aln)
        self.aln_bcw = backward(self._aln)
    
    def ref_transform(self, positions):
        frw = self.ref_frw
        bcw = self.aln_bcw
        
        result = []
        for pos in positions:
            result.append(bcw[frw[pos]])
        return result
    
    def aln_transform(self, positions):
        frw = self.aln_frw
        bcw = self.ref_bcw
        
        result = []
        for pos in positions:
            result.append(bcw[frw[pos]])
        return result

def align(first, second, path=None, save=None):
    if not path:
        path = ".align.fasta"
        save = False

    if not os.path.isfile(path):
        alnf = open(path, "w")
        alnf.write(">first\n" + first + "\n")
        alnf.write(">second\n" + second + "\n")
        alnf.close()

        mafft_cline = MafftCommandline(input=path)
        stdout, stderr = mafft_cline()

        outf = open(path, "w")
        outf.write(stdout)
        outf.close()
    
    reader = read_fasta(path)
    _, first = next(reader)
    _, second = next(reader)
    
    print(path)
    if not save:
        os.remove(path)
    
    return first.upper(), second.upper()

def forward(seq):
    # returns coordinates on aln seq
    coord = {}
    i = 0
    for j in range(len(seq)):
        if seq[j] == "-":
            continue
        else:
            coord[i] = j
            i += 1

    return coord

def backward(seq):
    # returns coordinates on ref seq
    coord = {}
    i = 0
    for j in range(len(seq)):
        coord[j] = i
        if seq[j] == "-":
            continue
        else:
            i += 1

    return coord

def read(qis, path):
    index = 0
    for covid, rna in read_fasta(path):
        qis[index].put((covid, rna, index))
        index += 1
        if index >= len(qis):
            index = 0

def work(path, qi):
    _, wuhan = next(read_fasta(WUHAN_PATH))
    while True:
        try:
            covid, rna, index = qi.get(timeout=2)
            align(
                wuhan, rna,
                path=path.rstrip("/") + "/" + covid.lstrip(">") + ".fasta",
                save=True
            )
        except:
            print(path)
            break

def process_align(
    gisaid_path,
    output_dir,
    process_num
):
    process_num = int(process_num)
    output_dir = output_dir.rstrip("/")

    qis = []
    workers = []
    for i in range(process_num):
        # qis.append(Queue(MAX_QUEUE_SIZE))
        qis.append(Queue())
        workers.append(Process(target=work, args=(
            output_dir, qis[-1],
        )))

    reader = Process(target=read, args=(
        qis, gisaid_path,
    ))

    reader.start()
    for i in range(process_num):
        workers[i].start()

    for i in range(process_num):
        workers[i].join()
    reader.join()
