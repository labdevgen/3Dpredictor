from bigWigFileReader import bigWigReader
from fastaFileReader import fastaReader
from shared import Interval, Genome
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
    faReader = fastaReader(path,name="hg38")
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

mem = memory_usage((test_fastaReader),interval=.5)
print ("Max memory: ",max(mem))
print("Memory log: ",mem)

# mem = memory_usage((test_bigWig, (False,)),interval=.5)
# print ("Max memory: ",max(mem))
# print("Memory log: ",mem)

#mem = memory_usage((test_bigWig, (True,)),interval=.5)
#print ("Max memory: ",max(mem))
#print("Memory log: ",mem)

#test_bigWig(inMem=False)

#test_bigWig(inMem=True)