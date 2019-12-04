import pandas as pd
import sys

def check_file_formats(RNAseq_file, CTCF_file, CTCF_orient_file):
    #check rna-seq file
    RNA_seq_data = pd.read_csv(RNAseq_file, sep="\t")
    if "chr" and "gene" and "FPKM" and "start" and "end" not in set(RNA_seq_data.keys()):
        print("RNA-seq file hasn't nessesary fields as chr,gene, FPKM, start, end", file=sys.stderr)
        sys.exit(1)
    #check CTCF file
    if CTCF_file.endswith(".gz"):  # check gzipped files
        import gzip
        temp_file = gzip.open(CTCF_file)
    else:
        temp_file = open(CTCF_file)
    Nfields = len(temp_file.readline().strip().split())
    temp_file.close()
    if Nfields<6:
        print("CTCF file has less than 6 columns", file=sys.stderr)
        sys.exit(1)
    #check CTCF orient data
    if CTCF_orient_file.endswith(".gz"):  # check gzipped files
        import gzip
        temp_file = gzip.open(CTCF_orient_file)
    else:
        temp_file = open(CTCF_orient_file)
    Nfields = len(temp_file.readline().strip().split())
    temp_file.close()
    if Nfields<5:
        print("CTCF_orient_file has less than 6 columns", file=sys.stderr)
        sys.exit(1)
    #assert CTCF_file contains data for orientation as +, -
    names = list(map(str, list(range(Nfields))))
    data = pd.read_csv(CTCF_orient_file, sep="\t", header=None, names=names)
    # subset and rename
    data = data.iloc[:, [0, 1, 2, 4, 5]]
    data.rename(columns={"0": "chr", "1": "start", "2": "end", "4": "score", "5": "orientation"},
                inplace=True)
    if "+" not in set(data["orientation"]) or "-" not in set(data["orientation"]):
        print("CTCF file orientation column has incorrect values which is inconsistent with + - values", file=sys.stderr)
        sys.exit(1)
