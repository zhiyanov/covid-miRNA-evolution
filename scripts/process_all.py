import sys

from lib.process import *
from lib.utils import *

PATH = "/home/dude/huge/dude/long-covid"


if __name__ == '__main__':
    proteins = [
        "",
        # "5P",
        # "NSP1",
        # "NSP2",
        # "NSP3",
        # "NSP4",
        # "NSP5",
        # "NSP6",
        # "NSP7",
        # "NSP8",
        # "NSP9",
        # "NSP10",
        # "NSP11",
        # "NSP12",
        # "NSP13",
        # "NSP14",
        # "NSP15",
        # "NSP16",
        # "Spike",
        # "NS3",
        # "E",
        # "M",
        # "NS6",
        # "NS7a",
        # "NS7b",
        # "NS8",
        # "N",
        # "NS9b",
        # "NS9c",
        # "3P",
    ]
    
    gisaids = [
        # "germany_ber",
        # "england_lond",
        "wuhan",
        # "taiwan",
    ]

    projects = [
        "",
        # "LUAD",
        # "COAD"
    ]
    
    for gisaid in gisaids:
        for project in projects:
            for protein in proteins:
                print(gisaid, project, protein)
                path = f"{PATH}/results/{gisaid}"
                if project:
                    path += "_" + project
                if protein:
                    path += "_" + protein
                if not os.path.exists(path):
                    os.makedirs(path)
            
                # process_spaces(
                process_seeds(
                    f"{PATH}/gisaid/{gisaid}.fasta",
                    None, # 0.05
                    path,
                    # 60,
                    1,
                    protein,
                    # f"{PATH}/miRNA/expressed_{project}.csv"
                )
