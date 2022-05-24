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
print("Data before inversion:")
print(rna_reader.chr_data['chr11'].iloc[10:20, [rna_reader.chr_data["chr11"].columns.get_loc("chr"),rna_reader.chr_data["chr11"].columns.get_loc("start")
                                                        ,rna_reader.chr_data["chr11"].columns.get_loc("end"),rna_reader.chr_data["chr11"].columns.get_loc("sigVal"),rna_reader.chr_data["chr11"].columns.get_loc("gene")]])

# inversion = Interval("chr11", 207501, 207505) #инверсия не задевает ни одни ген
# inversion = Interval("chr11", 207500, 207600) # конец инверсии задевает tss гена ENSG00000177963
# inversion = Interval("chr11", 207600, 209700) # начало и конец инверсии задевают tss гена ENSG00000177963
# inversion = Interval("chr11", 207500, 215200) # ген ENSG00000177963 полностью лежит внутри инверсии
# inversion = Interval("chr11", 207500, 236600) # инверсия не задевает ген ENSG00000177963, но затрагивает tss генов ENSG...142082 и ENSG...185627
# inversion = Interval("chr11", 207500, 253000) # гены ...177963, ...142082 и ...185627 лежат внутри инверсии

rna_reader.inverse_region_RNA(inversion, r"C:\Users\admin\3DPredictor-tests\tests\deletion_RNASeq_2\mart_export_txt")
print("Inversion interval:", inversion.start - inversion.end)
print("Data after inversion:")
print(rna_reader.chr_data['chr11'].iloc[10:20, [rna_reader.chr_data["chr11"].columns.get_loc("chr"),rna_reader.chr_data["chr11"].columns.get_loc("start")
                                                        ,rna_reader.chr_data["chr11"].columns.get_loc("end"),rna_reader.chr_data["chr11"].columns.get_loc("sigVal"),rna_reader.chr_data["chr11"].columns.get_loc("gene")]])
