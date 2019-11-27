from bigWigFileReader import bigWigReader
from fastaFileReader import fastaReader
from hicFileReader import hicReader
from shared import Interval, Genome

from straw import straw
import numpy as np

from memory_profiler import memory_usage
import datetime
import logging
logging.basicConfig(level=logging.DEBUG)

def make_genome():
    genome = Genome()
    genome.chrmSizes = {"chr1":248956422,
                        "chr2":242193529,
                        "chr3":198295559,
                        "chr4":190214555,
                        "chr5":181538259}
    return genome

genome = make_genome()

def test_bigWig(inMem):
    print ("Loading data")
    now = datetime.datetime.now()
    bwReader = bigWigReader("../input/ENCFF966IHQ.bigWig",name="Test",genome=genome)
    bwReader.readData(inMemory=inMem)
    print ("Time:",datetime.datetime.now() - now)

    print ("Extracting data, inMem=",str(inMem))
    now = datetime.datetime.now()
    start = 10000000
    stop = 101000000
    step = 1000000
    for i in range(start,stop,step):
        res = bwReader.get_interval(Interval("chr1",i,i+step))
    print ("Time:",datetime.datetime.now() - now)
    print (str(len(list(range(start,stop,step))))+" extractions of length "+str(step))

def test_fastaReader():
    print ("Loading data")
    now = datetime.datetime.now()
    path = "../input/hg38/hg38.fa"
    faReader = fastaReader(path,name="hg38",useOnlyChromosomes=["chr1"])
    print (faReader)
    print ("Time:",datetime.datetime.now() - now)

    print ("Extracting data")
    now = datetime.datetime.now()
    start = 10000000
    stop = 101000000
    step = 1000000
    for i in range(start,stop,step):
        res = faReader.get_interval(Interval("chr1",i,i+step))
    print ("Time:",datetime.datetime.now() - now)
    print (str(len(list(range(start,stop,step))))+" extractions of length "+str(step))

def get_hic_data():
    now = datetime.datetime.now()
    #path = "/home/minja/Desktop/hics/mouse_mESC_GSE82185_allValidPairs.hic"
    path = "../input/4DNFI2TK7L2F.hic"
    result = np.array(straw("KR", path,
                            "chr13","chr13","BP",1000))
    print("Time: ", datetime.datetime.now()-now)
    print (result)
    return result

def test_hicReader():
    genome = fastaReader("../input/hg38/test.fa",name="hg38")
    now = datetime.datetime.now()
    hic = hicReader(fname="../input/4DNFI2TK7L2F.hic", genome=genome, resolution = 100000)
    hic.read_data()
    print (hic.norms)
    result = hic.get_contact(Interval("chr1",0,120000000)) # single float value or NaN
    print (result)
    result = hic.get_chr_contact("chr1") # returns sparse matrix of the whole chrm as pandas dataframe

    print (datetime.datetime.now() - now)

now = datetime.datetime.now()
mem = memory_usage((test_hicReader),interval=.5)
print ("Max memory: ",max(mem))
print("Memory log: ",mem)
print("Total time: ",datetime.datetime.now() - now)

#mem = memory_usage((get_hic_data),interval=.5)
#print ("Max memory: ",max(mem))
#print("Memory log: ",mem)


#test_fastaReader()

#mem = memory_usage((test_fastaReader),interval=.5)
#print ("Max memory: ",max(mem))
#print("Memory log: ",mem)

# mem = memory_usage((test_bigWig, (False,)),interval=.5)
# print ("Max memory: ",max(mem))
# print("Memory log: ",mem)

#mem = memory_usage((test_bigWig, (True,)),interval=.5)
#print ("Max memory: ",max(mem))
#print("Memory log: ",mem)

#test_bigWig(inMem=False)

#test_bigWig(inMem=True)