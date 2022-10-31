import sys

from lib.process import *
from lib.utils import *


if __name__ == '__main__':
    input_path = sys.argv[1]
    output_dir = sys.argv[2].rstrip("/")
    threshold = float(sys.argv[3])
    
    if len(sys.argv) > 3:
        process_num = int(sys.argv[3])
    else:
        process_num = 1

    if len(sys.argv) > 4:
        protein = sys.argv[4]
    else:
        protein = None

    if len(sys.argv) > 5:
        miRNA_path = sys.argv[4]
    else:
        miRNA_path = None

    process(
        gisaid_path,
        threshold,
        output_dir,
        process_num,
        protein
        miRNA_path
    )
