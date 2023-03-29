from Bio import SeqIO
import tqdm
import sys


if __name__ == "__main__":
    path = sys.argv[1]
    reader = SeqIO.parse(path, "fasta")

    dct = {}
    for inp in tqdm.tqdm(reader, total=28420):
        date = inp.description.split("|")[1]
        if not date in dct:
            dct[date] = [inp]
        else:
            dct[date].append(inp)

    for date in tqdm.tqdm(sorted(dct), total=len(date)):
        for inp in dct[date]:
            print(">" + inp.description)
            print(inp.seq)


