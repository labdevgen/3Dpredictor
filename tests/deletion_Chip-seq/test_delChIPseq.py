#загружаем нужные библиотеки
from ChiPSeqReader import ChiPSeqReader
from shared import Interval


#говорим, что переменная ctcf_reader-объект класса ChipSeqReader
ctcf_reader = ChiPSeqReader(r"C:\Users\admin\3DPredictor-tests\tests\GSE96107_NPC_CTCF.IDR0.05.filt.narrowPeak",
                                                            name="CTCF")
#читаем файл
ctcf_reader.read_file()
#устанавливаем ориентация мотивов CTCF, для этого нужен доп файл с ориентациями
ctcf_reader.set_sites_orientation(r"C:\Users\admin\3DPredictor-tests\tests\GSE96107_NPC_CTCF.IDR0.05.filt.narrowPeak-orient.bed")
#напечатаем ту самую табличку, которая получилась в итоге
print(ctcf_reader.chr_data["chr11"].iloc[2043:2049,[ctcf_reader.chr_data["chr11"].columns.get_loc("chr"),ctcf_reader.chr_data["chr11"].columns.get_loc("start")
                                                       ,ctcf_reader.chr_data["chr11"].columns.get_loc("end"),ctcf_reader.chr_data["chr11"].columns.get_loc("sigVal")]])

#заводим переменную deletion и говорим что это объект класса Interval c параметрам, которые указаны в скобках
deletion = Interval("chr11", 108344000, 108414396)


#напечатаем ту самую табличку, которая получилась в итоге
# print(ctcf_reader.chr_data["chr11"])

#заводим переменную deletion и говорим что это объект класса Interval c параметрам, которые указаны в скобках
deletion = Interval("chr11", 108344000, 108414396)
#получим новую табличку которая включает только те строчки которые пересеклись
interval_data = ctcf_reader.get_interval(deletion)

print(interval_data[["chr", "start", "end", "sigVal", "minus_orientation"]])
st, end = ctcf_reader.get_interval(deletion, return_ids=True)
print(st, end)
print(ctcf_reader.chr_data["chr11"].iloc[2043:2049,[ctcf_reader.chr_data["chr11"].columns.get_loc("chr"),ctcf_reader.chr_data["chr11"].columns.get_loc("start")
                                                       ,ctcf_reader.chr_data["chr11"].columns.get_loc("end"),ctcf_reader.chr_data["chr11"].columns.get_loc("sigVal")]])
ctcf_reader.delete_region(deletion, 100)

print(deletion.len)
print(ctcf_reader.chr_data["chr11"].iloc[2043:2049,[ctcf_reader.chr_data["chr11"].columns.get_loc("chr"),ctcf_reader.chr_data["chr11"].columns.get_loc("start")
                                                       ,ctcf_reader.chr_data["chr11"].columns.get_loc("end"),ctcf_reader.chr_data["chr11"].columns.get_loc("sigVal")]])

