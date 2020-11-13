import sys, os
# Add main source directory to import path
root_dir = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
source_dir = os.path.join(root_dir,"source")
sys.path.append(source_dir)

# Add straw source dir



from fastaFileReader import fastaReader
from hicFileReader import hicReader
from shared import Interval, Genome
from bigWigFileReader import bigWigReader

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

def test_straw():
    import straw
    # data = straw("KR", "../input/4DNFI2TK7L2F.hic",
    #                                "19:8000000:10000000", "19:8000000:10000000", "BP", 1000000)
    # print (data)
    #straw_path = root_dir + "/straw/python/straw/straw.py"
    #import importlib.util
    #spec = importlib.util.spec_from_file_location("straw", straw_path)
    #foo = importlib.util.module_from_spec(spec)
    #spec.loader.exec_module(foo)
    #straw8 = foo.straw
    print (1)
    strawObj = straw.straw("../input/4DNFI2TK7L2F.hic")
    print ("19" in strawObj.chromDotSizes.data.keys())
    print (strawObj.chromDotSizes.getLength("19"))
    print (strawObj.chromDotSizes.getLength("19"))

    print (2)
    matrixObj = strawObj.getNormalizedMatrix("19", "19", "KR", "BP", 1000000)
    if matrixObj is None:
        print ("chr not found")
    print (3)
    data = matrixObj.getDataFromBinRegion(8,10,8,10)
    print (4)
    print (data)

def test_bigWig(inMem):
    print ("Loading data")
    now = datetime.datetime.now()
    bwReader = bigWigReader("../input/ENCFF966IHQ.bigWig",name="Test",genome=genome, inMemory=inMem)
    bwReader = bwReader.readData()
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
    path = "../input/hg38/test.fa"
    faReader = fastaReader(path,name="hg38",useOnlyChromosomes=["chr1"])
    faReader = faReader.read_data()
    print (faReader)
    print ("Time:",datetime.datetime.now() - now)

    # print ("Extracting data")
    # now = datetime.datetime.now()
    # start = 10000000
    # stop = 101000000
    # step = 1000000
    # for i in range(start,stop,step):
    #     res = faReader.get_interval(Interval("chr1",i,i+step))
    # print ("Time:",datetime.datetime.now() - now)
    # print (str(len(list(range(start,stop,step))))+" extractions of length "+str(step))

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
    faReader = genome.read_data()
    now = datetime.datetime.now()
    hic = hicReader(fname="../input/4DNFI2TK7L2F.hic", genome=genome, resolution = 100000)
    hic = hic.read_data()
    print (hic.norms)
    result = hic.get_contact(Interval("chr1",0,120000000)) # single float value or NaN
    print (result)
    result = hic.get_chr_contact("chr1") # returns sparse matrix of the whole chrm as pandas dataframe

    print (datetime.datetime.now() - now)

def simple_test():
    logging.basicConfig(level=logging.DEBUG) # set to INFO for less detailed output

    ### load data ###
    # load genome
    input_folder = "../input/"
    faReader = fastaReader(input_folder+"hg38/hg38.fa",useOnlyChromosomes=["chr1"])
    faReader = faReader.read_data()
    # load chipSeq
    bwReader1 = bigWigReader(input_folder+"ENCFF473IZV_H1_CTCF.bigWig", genome = faReader, inMemory=True)
    bwReader1 = bwReader1.readData()

    # load chipSeq
    bwReader2 = bigWigReader(input_folder+"ENCFF966IHQ.bigWig", genome = faReader, inMemory=False)
    bwReader2 = bwReader2.readData()


    #load contacts
    resolution = 5000
    hic = hicReader(input_folder+"4DNFI2TK7L2F.hic", genome=faReader, binsize = resolution,
                    indexedData = True)
    hic = hic.read_data()

    ### run simple check that contact count correlate with ChipSeq signal ###

    ### generate some random samples ####
    # get size of the chr1
    total_length = faReader.get_chr_sizes()["chr1"]

    window_size = 20*resolution # distance between intercting regions in this particular test, in units of resolution

    sample_size = 100000

    # select random points on chr1
    random_points_starts = np.random.random_integers(0,
                                              total_length-window_size,
                                              sample_size)
    random_points_starts = np.array((random_points_starts // resolution)*resolution,
                                    dtype = np.uint64)
    random_points_ends = random_points_starts + window_size

    # for each of selected points get contact between this point and (point + window_size*resolution)
    contacts = []
    chipSignals = []
    seqSignals = []
    now = datetime.datetime.now() # start timer

    logging.info("Starting data generation")
    for start,end in zip(random_points_starts,random_points_ends):
        interval = Interval("chr1",start,end)
        contact = hic.get_contact(interval)
        if contact == None or not np.isfinite(contact):
            continue
        else:
            chipSignal = np.nansum(bwReader1.get_interval(interval))
            if np.isfinite(chipSignal):
                chipSignals.append(chipSignal)
                seqSignal = np.sum(faReader.get_interval(interval))
                seqSignals.append(seqSignal)
                contacts.append(contact)

    logging.info("Time for data generation1: " + str(datetime.datetime.now() - now))
    # now = datetime.datetime.now()
    # chipSignals = []
    # seqSignals = []
    # contacts = []
    # for start,end in zip(random_points_starts,random_points_ends):
    #     interval = Interval("chr1",start,end)
    #     contact = hic.get_contact(interval)
    #     if contact == None:
    #         continue
    #     else:
    #         chipSignal = np.nansum(bwReader2.get_interval(interval))
    #         if np.isfinite(chipSignal):
    #             chipSignals.append(chipSignal)
    #             seqSignal = np.sum(faReader.get_interval(interval))
    #             seqSignals.append(seqSignal)
    #             contacts.append(contact)
    #
    # logging.info("Time for data generation2: " + str(datetime.datetime.now() - now))
    from scipy.stats import spearmanr
    import matplotlib.pyplot as plt

    print(contacts)
    print(chipSignals)

    print (spearmanr(np.array(contacts),np.array(chipSignals)))
    print (np.all(np.isfinite(contacts)))
    print (np.all(np.isfinite(chipSignals)))

    plt.scatter(contacts,chipSignals)
    plt.show()

def test_dump():
    faReader = fastaReader("../input/hg38/test.fa",useOnlyChromosomes=["chr2"])
    faReader = faReader.read_data()
    print (faReader.get_chr_sizes())
    # bwReader = bigWigReader("../input/ENCFF966IHQ.bigWig", genome = faReader, inMemory=True)
    # bwReader = bwReader.readData(debugMode=True)
    # bwReader.get_interval(Interval("chr1",100,10000))
    # print (bwReader.genome.name)
    #
    # hic = hicReader(fname="../input/4DNFI2TK7L2F.hic", genome=faReader, resolution = 100000)
    # hic = hic.read_data(debug_mode = True)
    # print(hic.get_contact(Interval("chr1",300000,500000)))

test_straw()

#simple_test()
#test_dump()

# now = datetime.datetime.now()
# mem = memory_usage((simple_test),interval=.5)
# print ("Max memory: ",max(mem))
# print("Memory log: ",mem)
# print("Total time: ",datetime.datetime.now() - now)
#

# now = datetime.datetime.now()
# mem = memory_usage((test_hicReader),interval=.5)
# print ("Max memory: ",max(mem))
# print("Memory log: ",mem)
# print("Total time: ",datetime.datetime.now() - now)

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