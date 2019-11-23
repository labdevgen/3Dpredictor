from bigWigFileReader import bigWigReader
from shared import Interval, Genome
import datetime
import logging
logging.basicConfig(level=logging.ERROR)

def make_genome():
    genome = Genome()
    genome.chrmSizes = {"chr1":248956422}
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
    stop = 100000000
    step = 1000000
    for i in range(start,stop,step):
        res = bwReader.get_interval(Interval("chr1",i,i+step))
    print ("Time:",datetime.datetime.now() - now)
    print (str(len(list(range(start,stop,step))))+" extractions of length "+str(step))

test_bigWig(inMem=False)
test_bigWig(inMem=True)