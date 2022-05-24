       # Test for the function "duplicate_region_RNA"
import pandas as pd
import sys
sys.path.insert(0, r'C:\Users\Maria\PycharmProjects\3DPredictor-study\3Dpredictor\source')
from RNASeqReader import RNAseqReader
from ChiPSeqReader import ChiPSeqReader
from shared import Interval

rna_reader = RNAseqReader(r"C:\Users\Maria\Desktop\Полина\duplication_RNASeq\rna-seqPolyA_tsvpre.txt", name="RNA")
rna_reader.read_file(rename={"Gene stable ID": "gene",
                                              "Gene start (bp)": "start",
                                              "Gene end (bp)": "end",
                                              "Chromosome/scaffold name": "chr",
                                              "FPKM": "sigVal"}, encoding='utf-8', sep='\t')
print("Data before duplication:")
print(rna_reader.chr_data['chr11'].iloc[10:20, [rna_reader.chr_data["chr11"].columns.get_loc("chr"),rna_reader.chr_data["chr11"].columns.get_loc("start")
                                                        ,rna_reader.chr_data["chr11"].columns.get_loc("end"),rna_reader.chr_data["chr11"].columns.get_loc("sigVal"),rna_reader.chr_data["chr11"].columns.get_loc("gene")]])
duplication = Interval("chr11", 207500, 253000)     # 1) 207400, 207500  2) 207500, 207600  3) 207600, 209700  4) 207500, 215200
rna_reader.duplicate_region_RNA(duplication, r"C:\Users\Maria\Desktop\mart_export_hg19.txt")   # 5) 207500, 236600  6) 207500, 253000
print(f"Duplication interval: {duplication.start} - {duplication.end}")
print("Data after duplication:")
print(rna_reader.chr_data['chr11'].iloc[10:20, [rna_reader.chr_data["chr11"].columns.get_loc("chr"),rna_reader.chr_data["chr11"].columns.get_loc("start")
                                                        ,rna_reader.chr_data["chr11"].columns.get_loc("end"),rna_reader.chr_data["chr11"].columns.get_loc("sigVal"),rna_reader.chr_data["chr11"].columns.get_loc("gene")]])
