import sys
import os

from Bio.Seq import Seq
from Bio.Align.Applications import MafftCommandline

from multiprocessing import Process, Queue

from .utils import read_fasta


class Align:
    def __init__(self, ref, aln, path=None):
        self.ref = ref
        self.aln = aln

        _ref, _aln = align(
            ref, aln,
            path=path,
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

def align(first, second, path=None):
    if not path:
        path = ".align.fasta"

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
