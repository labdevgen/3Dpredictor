import logging
from ChiPSeqReader import ChiPSeqReader
from Contacts_reader import ContactsReader
from shared import Interval,intersect_intervals
from matrix_plotter import MatrixPlotter
from E1_Reader import E1Reader, fileName2binsize
import matplotlib.pyplot as plt
import numpy as np


logging.basicConfig(level=logging.DEBUG)

def test_ctcf(): #comment
    ctcf_reader = ChiPSeqReader("C:/Users/POLINA/Desktop/lab.ICG/insulatory index/Hepat_WT_MboI_rep1-rep2.IDR0.05.filt.narrowPeak")
    ctcf_reader.read_file()
    d = ctcf_reader.get_interval(Interval("chr1",3448235,3456306))
    logging.info(d)
    d1 = ctcf_reader.get_binned_interval(Interval("chr1",3448200,3457000),binsize=500)
    logging.info(d1)
    d1 = ctcf_reader.get_nearest_peaks(Interval("chr1",3025000,3025000),N=5,side="left")
    logging.info(d1)
    print(ctcf_reader.chr_data['chr1'])

def test_ctcf_orient():
    orient_reader = ChiPSeqReader("D:/Users/Polina/3Dpredictor/input/Hepat_WT_MboI_rep1-rep2_IDR0_05_filt_narrowPeak-orient_N10.bed")
    orient_reader.read_orient_file()
    print(orient_reader.chr_data['chr1'])

def test_contacts():
    contacts_reader = ContactsReader()
    contacts_reader.read_files(["D:/Users/Polina/3Dpredictor/input/chr1.5MB.Hepat.contacts"])
    c = contacts_reader.get_contacts(Interval("chr1",5000000,6000000))
    logging.info(c)

def test_matrix_plot():
    contacts_reader = ContactsReader()
    contacts_reader.read_files(["input/chr1.5MB.Hepat.contacts"])
    c = contacts_reader.get_contacts(Interval("chr1",5000000,10000000))
    mp = MatrixPlotter()
    chr1contacts = contacts_reader.get_all_chr_contacts("chr1")
    logging.debug(chr1contacts.head)
    mp.set_data(chr1contacts)
    m = mp.getMatrix4plot(Interval("chr1",5000000,10000000))
    m = np.log(m)
    plt.imshow(m,cmap="OrRd")
    plt.show()

def test_E1reader():
    files = ["input/chr1.Hepat.E1.50k",
             "input/chr2.Hepat.E1.50k"]
    eig = E1Reader()
    eig.read_files(fnames=files,binSizeFromName = fileName2binsize)
    print(eig.get_E1inInterval(Interval("chr1",1,200000)))
    print("-----------------")
    print(eig.get_E1inInterval(Interval("chr1",194600000,195600000)))
    print("-----------------")
    print(eig.get_E1inInterval(Interval("chr1",189500000,190500000)))
    print("-----------------")
    #print(eig.get_E1inInterval(Interval("chr1",200000000,250000000)))

def test_ChipSeqRemoval():
    ctcf_reader = ChiPSeqReader("input/Hepat_WT_MboI_rep1-rep2.IDR0.05.filt.narrowPeak")
    ctcf_reader.read_file()
    logging.info("------------Before deleting:")

    interval = Interval("chr1",3448235,3700000)
    logging.info(interval)
    d = ctcf_reader.get_interval(interval)
    logging.info(d)

    interval = Interval("chr1", 3448235, 3900000)
    logging.info(interval)
    d = ctcf_reader.get_interval(Interval("chr1",3448235,3900000))
    logging.info(d)

    logging.info("----after deleting----")
    ctcf_reader.delete_region(Interval("chr1",3454000,3611129))

    interval = Interval("chr1",3448235,3700000)
    logging.info(interval)
    d = ctcf_reader.get_interval(interval)
    logging.info(d)

    interval = Interval("chr1", 3448235, 3900000)
    logging.info(interval)
    d = ctcf_reader.get_interval(Interval("chr1",3448235,3900000))
    logging.info(d)


def test_ContactsRemoval():
    contacts_reader = ContactsReader()
    contacts_reader.read_files(["input/chr1.5MB.Hepat.contacts"])
    c = contacts_reader.get_contacts(Interval("chr1",5000000,5150000))
    logging.info(c)
    contacts_reader.delete_region(Interval("chr1",5030000,5100000))
    c = contacts_reader.get_contacts(Interval("chr1",5000000,5150000))
    logging.info(c)

def test_intersect_intervals():
    orient_reader = ChiPSeqReader("D:/Users/Polina/3Dpredictor/input/Hepat_WT_MboI_rep1-rep2_IDR0_05_filt_narrowPeak-orient_N10.bed")
    orient_reader.read_orient_file()
    ctcf_reader = ChiPSeqReader("D:/Users/Polina/3Dpredictor/input/Hepat_WT_MboI_rep1-rep2.IDR0.05.filt.narrowPeak")
    ctcf_reader.read_file()
    print(ctcf_reader.chr_data['chr4'])
    print(orient_reader.chr_data['chr4'])
    result = intersect_intervals(ctcf_reader.chr_data, orient_reader.chr_data)
    print(result['chr4'].iloc[831])
    #print(ctcf_reader.chr_data['chr1'].iloc[75])
    #print(ctcf_reader.chr_data['chr1'].iloc[76])
def test_sites_orientation():
    ctcf_reader = ChiPSeqReader(
        "D:/Users/Polina/3Dpredictor/input/Hepat_WT_MboI_rep1-rep2.IDR0.05.filt.narrowPeak")
    ctcf_reader.read_file()
    orient_fname = "D:/Users/Polina/3Dpredictor/input/Hepat_WT_MboI_rep1-rep2_IDR0_05_filt_narrowPeak-orient_N10.bed"
    result = ChiPSeqReader.set_sites_orientation(ctcf_reader, orient_fname)
    print(result.chr_data['chr1'])
    print(result.chr_data['chr1'].iloc[0])

def test_add_orientation():
    ctcf_reader = ChiPSeqReader(
        "D:/Users/Polina/3Dpredictor/input/Hepat_WT_MboI_rep1-rep2.IDR0.05.filt.narrowPeak")
    ctcf_reader.read_file()
    orient_fname = "D:/Users/Polina/3Dpredictor/input/Hepat_WT_MboI_rep1-rep2_IDR0_05_filt_narrowPeak-orient_N10.bed"
    ctcf_reader.add_orient_to_chr_data(orient_fname)



#test_add_orientation()
#test_sites_orientation()
#test_intersect_intervals()
#test_matrix_plot()
#test_contacts()
#test_E1reader()
#test_ctcf()
#test_ChipSeqRemoval()
#test_ContactsRemoval() #TODO it doesn't throw errors, however the behaviour was not thoroughly tested
test_ctcf_orient()