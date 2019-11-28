from bigWigFileReader import bigWigReader
from fastaFileReader import fastaReader
from hicFileReader import hicReader
from shared import Interval
import numpy as np
import datetime
import logging
logging.basicConfig(level=logging.DEBUG)


def calc_corr():
    logging.basicConfig(level=logging.DEBUG) # set to INFO for less detailed output

    ### load data ###
    # load genome
    chr = "chr2"
    faReader = fastaReader("../input/hg38/hg38.fa",useOnlyChromosomes=[chr])
    faReader = faReader.read_data()

    # load chipSeq
    bwReader1 = bigWigReader("../input/ENCFF966IHQ.bigWig", genome = faReader, inMemory=True)
    bwReader1 = bwReader1.readData()


    #load contacts
    resolution = 5000
    hic = hicReader("../input/4DNFI2TK7L2F.hic", genome=faReader, resolution = resolution)
    hic = hic.read_data()

    ### run simple check that contact count correlate with ChipSeq signal ###

    ### generate some random samples ####
    # get size of the chr1
    total_length = faReader.get_chr_sizes()[chr]

    window_size = 10*resolution # distance between intercting regions in this particular test, in units of resolution
    sample_size = 1000

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
        window = Interval(chr,int(start+2*resolution),int(end-2*resolution))
        contact = hic.get_contact(interval)
        if contact == None:
            contact = 0
        if np.isfinite(contact):
            chipSignal = np.nansum(bwReader1.get_interval(window))
            if np.isfinite(chipSignal):
                chipSignals.append(chipSignal)
                seqSignal = np.sum(faReader.get_interval(interval))
                seqSignals.append(seqSignal)
                contacts.append(contact)

    logging.info("Time for data generation: " + str(datetime.datetime.now() - now))
    from scipy.stats import spearmanr, pearsonr

    print (pearsonr(np.array(contacts),np.array(chipSignals)))
    print (pearsonr(np.array(contacts),np.array(seqSignals)))

    # import matplotlib.pyplot as plt
    # plt.scatter(contacts,chipSignals)
    # plt.show()
