import sys
import pandas as pd

from lib.compare import compare

if __name__ == '__main__':
    output_path = sys.argv[1]
    project = sys.argv[2]

    proteins = [
        "NSP1", "NSP2", "NSP3", "NSP4",
        "NSP5", "NSP6", "NSP7", "NSP8",
        "NSP9", "NSP10", "NSP11", "NSP12",
        "NSP13", "NSP14", "NSP15", "NSP16",
        "Spike", "NS8"
    ]

    paths = [f"../results/russia_mos_{project}/result.csv"]
    # paths = [f"../results/russia_mos/result.csv"]
    for protein in proteins:
        paths.append(f"../results/russia_mos_{project}_{protein}/result.csv")
        # paths.append(f"../results/russia_mos_{protein}/result.csv")

    print(project)
    result = compare(*paths)
    
    if len(sys.argv) > 3:
        miRNA_path = sys.argv[3]
        miRNAs = pd.read_csv(miRNA_path, sep=",")
        miRNAs.index = miRNAs["MIMAT"]
        miRNAs = miRNAs.drop(columns=["MIMAT"])
        result = result.join(miRNAs, how="left").fillna(0.)

        result.to_csv(output_path, sep=",")
