import logging
from ChiPSeqReader import ChiPSeqReader
from Contacts_reader import ContactsReader
from shared import Interval
from matrix_plotter import MatrixPlotter
import matplotlib.pyplot as plt
import numpy as np


logging.basicConfig(level=logging.DEBUG)

def test_ctcf(): #comment
    ctcf_reader = ChiPSeqReader("D:/Lab Archive/ForChrRearrModel/Hepat_WT_MboI_rep1-rep2.IDR0.05.filt.narrowPeak")
    ctcf_reader.read_file()
    d = ctcf_reader.get_interval(Interval("chr1",3448235,3456306))
    logging.info(d)
    d1 = ctcf_reader.get_binned_interval(Interval("chr1",3448200,3457000),binsize=500)
    logging.info(d1)

def test_contacts():
    contacts_reader = ContactsReader()
    contacts_reader.read_files(["C:/Users/FishmanVS/Desktop/RNF3D_beds/chr1.Hepat.contacts"])
    c = contacts_reader.get_contacts(Interval("chr1",5000000,6000000))
    logging.info(c)

def test_matrix_plot():
    contacts_reader = ContactsReader()
    contacts_reader.read_files(["C:/Users/FishmanVS/Desktop/RNF3D_beds/chr1.5MB.Hepat.contacts"])
    #c = contacts_reader.get_contacts(Interval("chr1",5000000,10000000))
    mp = MatrixPlotter()
    chr1contacts = contacts_reader.get_all_chr_contacts("chr1")
    mp.set_data(chr1contacts)
    m = mp.getMatrix4plot(Interval("chr1",5000000,10000000))
    m = np.log(m)
    plt.imshow(m,cmap="OrRd")
    plt.show()

test_matrix_plot()