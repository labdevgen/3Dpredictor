from bigWigFileReader import bigWigReader
from fastaFileReader import fastaReader
from hicFileReader import hicReader
from shared import Interval
import numpy as np
import datetime
import logging
logging.basicConfig(level=logging.DEBUG)

import sys, os

# Add main source directory to import path
root_dir = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
source_dir = os.path.join(root_dir,"source")
sys.path.append(source_dir)

def calc_corr(chr, resolution = 5000, window_size = 20):
    logging.basicConfig(level=logging.DEBUG) # set to INFO for less detailed output

    ### load data ###
    # load genome
    faReader = fastaReader("../input/hg38/hg38.fa",useOnlyChromosomes=[chr])
    faReader = faReader.read_data()

    # load chipSeq1
    bwReader1 = bigWigReader("../input/ENCFF473IZV_H1_CTCF.bigWig", genome = faReader, inMemory=True)
    bwReader1 = bwReader1.readData()


    #load contacts

    hic = hicReader("../input/4DNFI2TK7L2F.hic", genome=faReader, resolution = resolution)
    hic = hic.read_data()

    ### run simple check that contact count correlate with ChipSeq signal ###

    ### generate some random samples ####
    # get size of the chr1
    total_length = faReader.get_chr_sizes()[chr]

     # distance between intercting regions in this particular test, in units of resolution
    sample_size = 5000

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
        interval = Interval(chr,start,end)
        assert window_size >=5*resolution
        window = Interval(chr,start+resolution,end)
        contact = hic.get_contact(interval)
        if contact == None:
            contact = 0
        if np.isfinite(contact):
            # chipSignal = np.concatenate((bwReader1.get_interval(Interval(chr,int(start-resolution),int(start+resolution))),
            #                             bwReader1.get_interval(
            #                                 Interval(chr, int(end - resolution), int(end + resolution)))))
            chipSignal = bwReader1.get_interval(window)
            chipSignal = np.nan_to_num(chipSignal)
            chipSignal = np.sum(chipSignal)
            if np.isfinite(chipSignal):
                chipSignals.append(chipSignal)
                seqSignal = np.sum(faReader.get_interval(interval))
                seqSignals.append(seqSignal)
                contacts.append(contact)

    logging.info("Time for data generation: " + str(datetime.datetime.now() - now))
    from scipy.stats import spearmanr, pearsonr

    res = []
    res.append(spearmanr(np.array(contacts),np.array(chipSignals))[0])
    res.append(pearsonr(np.array(contacts),np.array(chipSignals))[0])
    res.append(spearmanr(np.array(contacts),np.array(seqSignals))[0])
    res.append(pearsonr(np.array(contacts),np.array(seqSignals))[0])

    return ("\t".join(list(map(str,res))))
    # import matplotlib.pyplot as plt
    # plt.scatter(contacts,chipSignals)
    # plt.show()

def calc_sparsity():
    logging.basicConfig(level=logging.DEBUG) # set to INFO for less detailed output

    ### load data ###
    # load genome
    chr = "chr2"
    faReader = fastaReader("../input/hg38/hg38.fa",useOnlyChromosomes=[chr])
    faReader = faReader.read_data()

    # load chipSeq
    bwReader1 = bigWigReader("../input/ENCFF473IZV_H1_CTCF.bigWig", genome = faReader, inMemory=True)
    bwReader1 = bwReader1.readData()

    arr = bwReader1.data[chr]
    print(len(arr))
    nonzero = arr[np.nonzero(arr)]
    print(len(nonzero))
    finite = nonzero[np.isfinite(nonzero)]
    print(len(finite))

def calc_insulation_around_CTCF(chr, resolution = 5000, window_size = 20):
    logging.basicConfig(level=logging.DEBUG) # set to INFO for less detailed output

    ### load data ###
    # load genome
    faReader = fastaReader("../input/hg38/hg38.fa",useOnlyChromosomes=[chr])
    faReader = faReader.read_data()

    # load chipSeq1
    bwReader1 = bigWigReader("../input/ENCFF473IZV_H1_CTCF.bigWig", genome = faReader, inMemory=True)
    bwReader1 = bwReader1.readData()


    #load contacts

    hic = hicReader("../input/4DNFI2TK7L2F.hic", genome=faReader, resolution = resolution)
    hic = hic.read_data()

    ### run simple check that contact count correlate with ChipSeq signal ###

    ### generate some random samples ####
    # get size of the chr1
    total_length = faReader.get_chr_sizes()[chr]

    all_CTCF = bwReader1.get_interval(Interval(chr,0,total_length))
    all_CTCF = np.nan_to_num(all_CTCF)
    binsize = 1000
    bins = np.arange(0,total_length-1,binsize)
    sums = [np.sum(all_CTCF[a:a+binsize]) for a in bins]
    peaks = bins[sums>np.percentile(sums,90)]
    with open("../out/test.bed","w") as fout:
        for i in peaks:
            fout.write(chr+"\t"+str(i)+"\t"+str(i+binsize)+"\n")

calc_insulation_around_CTCF(chr="chr3",resolution=5000)

# resolution = 5000
# for chr in ["chr1","chr2","chr3"]:
#     for winsize in [5*resolution,10*resolution,15*resolution,20*resolution]:
#         print(chr,winsize,resolution,calc_corr(chr,resolution=resolution,window_size=winsize))
#calc_corr()
#calc_sparsity()