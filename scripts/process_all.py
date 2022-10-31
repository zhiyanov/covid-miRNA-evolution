import sys

from lib.process import *
from lib.utils import *


if __name__ == '__main__':
    proteins = [
        "5P",
        "3P",
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
        # ""
    ]

    projects = [
        "LUAD",
        "COAD"
    ]

    for project in projects:
        for protein in proteins:
            print(project, protein)
            # print(protein)
            if protein:
                if not os.path.exists(f"../results/russia_mos_{project}_{protein}"):
                    os.makedirs(f"../results/russia_mos_{project}_{protein}")
                output_dir = f"../results/russia_mos_{project}_{protein}"
                # if not os.path.exists(f"../results/russia_mos_{protein}"):
                #     os.makedirs(f"../results/russia_mos_{protein}")
                # output_dir = f"../results/russia_mos_{protein}"
            else:
                if not os.path.exists(f"../results/russia_mos_{project}"):
                    os.makedirs(f"../results/russia_mos_{project}")
                output_dir = f"../results/russia_mos_{project}"
                # if not os.path.exists(f"../results/russia_mos"):
                #     os.makedirs(f"../results/russia_mos")
                # output_dir = f"../results/russia_mos"
        
            # process_spaces(
            process_seeds(
                "../gisaid/russia_mos.fasta",
                None, # 0.05
                output_dir,
                60,
                protein,
                f"../miRNA/expressed_{project}.csv"
            )
