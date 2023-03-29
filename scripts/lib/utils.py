import re
import io


def find_all(seq, ref):
    return [m.start() for m in re.finditer(seq, ref)]

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
            description = line
            result = ""
        else:
            result += line

        line = next(ff, None)
    
    if description and result:
        yield description, result

    ff.close()

def reverse(seq):
    result = ""
    for char in seq[::-1]:
        if char == "A":
            result += "T"
        if char == "T":
            result += "A"
        if char == "G":
            result += "C"
        if char == "C":
            result += "G"
    return result

def uremove(seq):
    result = ""
    for char in seq:
        if char == "U":
            result += "T"
        else:
            result += char
    return result

def percent(seq, chr="N"):
    count = 0
    for char in seq:
        if char == chr:
            count += 1
    return count / len(seq)

def region_count(regions):
    count = regions.count(",")
    if count > 0:
        return count + 1
    return count

