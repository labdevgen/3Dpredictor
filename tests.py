import logging
from ChiPSeqReader import ChiPSeqReader
from RNASeqReader import RNAseqReader
from Contacts_reader import ContactsReader
from shared import Interval,intersect_intervals
from matrix_plotter import MatrixPlotter
from E1_Reader import E1Reader, fileName2binsize
import matplotlib.pyplot as plt
import numpy as np
from DataGenerators import PredictorGenerator, SitesOrientPredictorGenerator
import pandas as pd
import scipy.stats

logging.basicConfig(level=logging.DEBUG)

def test_ctcf(): #comment
    ctcf_reader = ChiPSeqReader("C:/Users/POLINA/Desktop/lab.ICG/insulatory index/Hepat_WT_MboI_rep1-rep2.IDR0.05.filt.narrowPeak")
    ctcf_reader.read_file()
    d = ctcf_reader.get_interval(Interval("chr1",3448235,3456306))
    #logging.getLogger(__name__).info(d)
    logging.info(d)
    d1 = ctcf_reader.get_binned_interval(Interval("chr1",3448200,3457000),binsize=500)
    logging.getLogger(__name__).info(d1)
    #d1 = ctcf_reader.get_nearest_peaks(Interval("chr1",3025000,3025000),N=5,side="left")
    logging.getLogger(__name__).info(d1)
    logging.info(d1)
    d1 = ctcf_reader.get_nearest_peaks(Interval("chr1",3025000,3025000),N=5,side="left")
    logging.info(d1)
    print(ctcf_reader.chr_data['chr1'])

def test_read_orient():
    ctcf_reader = ChiPSeqReader("D:/Users/Polina/3Dpredictor/input/Hepat_WT_MboI_rep1-rep2.IDR0.05.filt.narrowPeak")
    ctcf_reader.read_file()
    orient_data = ctcf_reader.read_orient_file("D:/Users/Polina/3Dpredictor/input/Hepat_WT_MboI_rep1-rep2_IDR0_05_filt_narrowPeak-orient_N10.bed")
    print(orient_data['chr1'])

def test_contacts():
    contacts_reader = ContactsReader()
    contacts_reader.read_files(["D:/Users/Polina/3Dpredictor/input/chr1.5MB.Hepat.contacts"])
    c = contacts_reader.get_contacts(Interval("chr1",5000000,6000000))
    logging.getLogger(__name__).info(c)

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
    logging.getLogger(__name__).info("------------Before deleting:")

    interval = Interval("chr1",3448235,3700000)
    logging.getLogger(__name__).info(interval)
    d = ctcf_reader.get_interval(interval)
    logging.getLogger(__name__).info(d)

    interval = Interval("chr1", 3448235, 3900000)
    logging.getLogger(__name__).info(interval)
    d = ctcf_reader.get_interval(Interval("chr1",3448235,3900000))
    logging.getLogger(__name__).info(d)

    logging.getLogger(__name__).info("----after deleting----")
    ctcf_reader.delete_region(Interval("chr1",3454000,3611129))

    interval = Interval("chr1",3448235,3700000)
    logging.getLogger(__name__).info(interval)
    d = ctcf_reader.get_interval(interval)
    logging.getLogger(__name__).info(d)

    interval = Interval("chr1", 3448235, 3900000)
    logging.getLogger(__name__).info(interval)
    d = ctcf_reader.get_interval(Interval("chr1",3448235,3900000))
    logging.getLogger(__name__).info(d)


def test_ContactsRemoval():
    contacts_reader = ContactsReader()
    contacts_reader.read_files(["input/chr1.5MB.Hepat.contacts"])
    c = contacts_reader.get_contacts(Interval("chr1",5000000,5150000))
    logging.getLogger(__name__).info(c)
    contacts_reader.delete_region(Interval("chr1",5030000,5100000))
    c = contacts_reader.get_contacts(Interval("chr1",5000000,5150000))
    logging.getLogger(__name__).info(c)

def test_intersect_intervals():
    ctcf_reader = ChiPSeqReader("D:/Users/Polina/3Dpredictor/input/Hepat_WT_MboI_rep1-rep2.IDR0.05.filt.narrowPeak")
    ctcf_reader.read_file()
    orient_data = ctcf_reader.read_orient_file("D:/Users/Polina/3Dpredictor/input/Hepat_WT_MboI_rep1-rep2_IDR0_05_filt_narrowPeak-orient_N10.bed")
    print(ctcf_reader.chr_data['chr4'])
    print(orient_data['chr4'])
    #print(orient_data['chr4'].loc[orient_data['chr4']['start'] == 117163608])
    result = intersect_intervals(ctcf_reader.chr_data, orient_data)
    #print(result['chr4'].loc[result['chr4']['start'] == 117163608])
    print(result['chr4'])
    #print(ctcf_reader.chr_data['chr1'].iloc[75])
    #print(ctcf_reader.chr_data['chr1'].iloc[76])
def test_sites_orientation():
    ctcf_reader = ChiPSeqReader(
        "D:/Users/Polina/3Dpredictor/input/Hepat_WT_MboI_rep1-rep2.IDR0.05.filt.narrowPeak")
    ctcf_reader.read_file()
    ctcf_reader.set_sites_orientation("D:/Users/Polina/3Dpredictor/input/Hepat_WT_MboI_rep1-rep2_IDR0_05_filt_narrowPeak-orient_N10.bed")
    print(ctcf_reader.chr_data['chr1'])
    #print(ctcf_reader.chr_data['chr4'].iloc[23])
    #print(ctcf_reader.chr_data['chr4'])
    #ctcf_reader.export2bed_files_with_orientation("D:/Users/Polina/3Dpredictor/data/")
    result = ctcf_reader.get_only_with_orient_data()
    print(result.chr_data['chr1'])
    print(result.chr_data['chr1'].query("start=='4516413'"))
    #print(result)
def test_N_nearest_peaks_in_interval():
    ctcf_reader = ChiPSeqReader(
        "D:/Users/Polina/3Dpredictor/input/Hepat_WT_MboI_rep1-rep2.IDR0.05.filt.narrowPeak")
    ctcf_reader.read_file()
    ctcf_reader.set_sites_orientation(
        "D:/Users/Polina/3Dpredictor/input/Hepat_WT_MboI_rep1-rep2_IDR0_05_filt_narrowPeak-orient_N10.bed")
    #print(ctcf_reader.chr_data['chr1'])
    result = ctcf_reader.get_N_nearest_peaks_in_interval(interval = Interval("chr1", 100800000, 101125000 ), N=6)
    print('-----------------------------------------')
    #print(result)
    print('sumr', result[0].sigVal.sum())
    print(result[1])
def test_get_nearest_peaks():
    ctcf_reader = ChiPSeqReader(
        "D:/Users/Polina/3Dpredictor/input/Hepat_WT_MboI_rep1-rep2.IDR0.05.filt.narrowPeak")
    ctcf_reader.read_file()
    ctcf_reader.set_sites_orientation(
        "D:/Users/Polina/3Dpredictor/input/Hepat_WT_MboI_rep1-rep2_IDR0_05_filt_narrowPeak-orient_N10.bed")
    print(ctcf_reader.chr_data['chr1'])
    result = ctcf_reader.get_nearest_peaks(Interval("chr1", 3611433, 3611433), N=6, side='left')
    print(result)
def test_get_interval():
    ctcf_reader = ChiPSeqReader(
        "D:/Users/Polina/3Dpredictor/input/Hepat_WT_MboI_rep1-rep2.IDR0.05.filt.narrowPeak")
    ctcf_reader.read_file()
    ctcf_reader.set_sites_orientation(
        "D:/Users/Polina/3Dpredictor/input/Hepat_WT_MboI_rep1-rep2_IDR0_05_filt_narrowPeak-orient_N10.bed")
    interval = ctcf_reader.get_interval(interval = Interval("chr1", 100800000, 101125000))
    print(interval)
def correlation():
    training_data = pd.read_csv("2018-09-17-trainingOrient.RandOnChr1.20000.contacts.3000000.50001.500000.25000.txt",delimiter="\t")
    training_data.fillna(value=0, inplace=True)
    print(training_data['CTCForient_W_sumSigVal'])
    print(training_data['CTCF_W'])
    res = scipy.stats.spearmanr(training_data['CTCForient_W_sumSigVal'], training_data['CTCF_W'])
    print(res)

def test_RNAseqReader():
    RNA = RNAseqReader(fname="input/GSE95111_genes.fpkm_table.txt.pre.txt")
    RNA.read_file(rename={"Gene name":"gene",
                          "Gene start (bp)":"start",
                          "Gene end (bp)":"end",
                          "Chromosome/scaffold name":"chr",
                          "shCtrl-1_0":"sigVal"},
                  sep="\t")
    logging.getLogger(__name__).info(iter(RNA.chr_data.values()).__next__().head())
    for inteval2 in [Interval("chr1",36511867,36528237),
                        Interval("chr1",36511866,36528237),
                        Interval("chr1",36511867, 36528238),
                        Interval("chr1",36511866, 36528238),
                        Interval("chr1",36528238, 36528238),
                        Interval("chr1",36528238, 36528239),
                        Interval("chr1",36528138, 36528149),
                        Interval("chr1",36500000, 36550000)]:
        logging.info("------------------")
        logging.info(inteval2)
        logging.info(str(RNA.get_interval(inteval2)))


#correlation()
#test_get_interval()
#test_ori_predictor_generator()
#test_get_nearest_peaks()
#test_N_nearest_peaks_in_interval()
#test_add_orientation()
#test_sites_orientation()
#test_intersect_intervals()
#test_matrix_plot()
#test_contacts()
#test_E1reader()
#test_ctcf()
#test_ChipSeqRemoval()
#test_ContactsRemoval() #TODO it doesn't throw errors, however the behaviour was not thoroughly tested
#test_read_orient()
test_RNAseqReader()