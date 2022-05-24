import sys
sys.path.insert(0, r'C:\Users\admin\3DPredictor-tests\3Dpredictor\source')
from RNASeqReader import RNAseqReader
from ChiPSeqReader import ChiPSeqReader
from shared import Interval

rna_reader = RNAseqReader(r"C:\Users\admin\3DPredictor-tests\3Dpredictor\input\K562\RNA-seq\rna-seqPolyA.tsvpre.txt", name="RNA")
rna_reader.read_file(rename={"Gene stable ID": "gene",
                                              "Gene start (bp)": "start",
                                              "Gene end (bp)": "end",
                                              "Chromosome/scaffold name": "chr",
                                              "FPKM": "sigVal"}, encoding='utf-8', sep='\t')
print("Data before deletion:")
print(rna_reader.chr_data['chr11'].iloc[10:20, [rna_reader.chr_data["chr11"].columns.get_loc("chr"),rna_reader.chr_data["chr11"].columns.get_loc("start")
                                                        ,rna_reader.chr_data["chr11"].columns.get_loc("end"),rna_reader.chr_data["chr11"].columns.get_loc("sigVal"),rna_reader.chr_data["chr11"].columns.get_loc("gene")]])
deletion = Interval("chr11", 207500, 253000)     # 1) 207501, 207505  2) 207500, 207600  3) 207600, 209700  4) 207500, 215200
rna_reader.delete_region_RNA(deletion, r"C:\Users\admin\3DPredictor-tests\tests\deletion_RNASeq_2\mart_export_txt")   # 5) 207500, 236600  6) 207500, 253000
print("Deletion interval:", deletion.start - deletion.end)
print("Data after deletion:")
print(rna_reader.chr_data['chr11'].iloc[10:20, [rna_reader.chr_data["chr11"].columns.get_loc("chr"),rna_reader.chr_data["chr11"].columns.get_loc("start")
                                                        ,rna_reader.chr_data["chr11"].columns.get_loc("end"),rna_reader.chr_data["chr11"].columns.get_loc("sigVal"),rna_reader.chr_data["chr11"].columns.get_loc("gene")]])
