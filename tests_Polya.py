import logging
from ChiPSeqReader import ChiPSeqReader
from RNASeqReader import RNAseqReader
from Contacts_reader import ContactsReader
from shared import Interval,intersect_intervals
from matrix_plotter import MatrixPlotter
from E1_Reader import E1Reader, fileName2binsize
import matplotlib.pyplot as plt
import numpy as np
from PredictorGenerators import PredictorGenerator, SitesOrientPredictorGenerator
import pandas as pd
import scipy.stats
from Predictor import Predictor
from Weight_funcs_modul import *
from add_loop import add_loop

logging.basicConfig(level=logging.DEBUG)

def test_ctcf(): #comment
    ctcf_reader = ChiPSeqReader("input/GM12878/")
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
    print("!!!!!!1")
    ctcf_reader = ChiPSeqReader("input/Hepat/CTCF/Hepat_WT_MboI_rep1-rep2.IDR0.05.filt.narrowPeak")
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
    contacts_reader.read_files(["input/Hepat/chr1.5MB.Hepat.contacts.gz"], coeff_fname="cofficient_Hepat.txt")
    c = contacts_reader.get_contacts(Interval("chr1",5000000,5150000))
    logging.getLogger(__name__).info(c)
    contacts_reader.delete_region(Interval("chr1",5031234,5100000))
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
    ctcf_reader_2 = ChiPSeqReader("D:/Users/Polina/3Dpredictor/input/NPC/CTCF/GSE96107_NPC_CTCF.IDR0.05.filt.narrowPeak")
    ctcf_reader_2.read_file()
    ctcf_reader_2.set_sites_orientation("D:/Users/Polina/3Dpredictor/input/NPC/CTCF/GSE96107_NPC_CTCF.IDR0.05.filt.narrowPeak-orient.bed")
    print(ctcf_reader_2.chr_data["chr1"])
    ctcf_reader_2.chr_data["chr1"].to_csv("D:/Users/Polina/3Dpredictor/input/NPC/CTCF/testNPC", sep="/t")
    ctcf_reader_1 = ChiPSeqReader(
        "D:/Users/Polina/3Dpredictor/input/Hepat/CTCF/Hepat_WT_MboI_rep1-rep2.IDR0.05.filt.narrowPeak")
    ctcf_reader_1.read_file()
    ctcf_reader_1.set_sites_orientation("D:/Users/Polina/3Dpredictor/input/Hepat/CTCF/Hepat_WT_MboI_rep1-rep2_IDR0_05_filt_narrowPeak-orient_N10.bed")
    print(ctcf_reader_1.chr_data["chr1"])
    ctcf_reader_1.chr_data["chr1"].to_csv("D:/Users/Polina/3Dpredictor/input/Hepat/CTCF/testHepat", sep="\t")

    #print(ctcf_reader.chr_data['chr1'])
    #print(ctcf_reader.chr_data['chr4'].iloc[23])
    #print(ctcf_reader.chr_data['chr4'])
    #ctcf_reader.export2bed_files_with_orientation("D:/Users/Polina/3Dpredictor/data/")
    # ctcf_reader.keep_only_with_orient_data()
    #print(ctcf_reader.chr_data['chr1'])
    #print(ctcf_reader.chr_data['chr1'].query("start=='4516413'"))
    #print(result)
# test_sites_orientation()

def test_N_nearest_peaks_in_interval():
    ctcf_reader = ChiPSeqReader(
        "D:/Users/Polina/3Dpredictor/input/Hepat_WT_MboI_rep1-rep2.IDR0.05.filt.narrowPeak")
    ctcf_reader.read_file()
    ctcf_reader.set_sites_orientation(
        "D:/Users/Polina/3Dpredictor/input/Hepat_WT_MboI_rep1-rep2_IDR0_05_filt_narrowPeak-orient_N10.bed")
    #print(ctcf_reader.chr_data['chr1'])
    result = ctcf_reader.get_N_peaks_near_interval_boundaries(interval = Interval("chr1", 100800000, 101125000 ), N=6)
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
    RNA = RNAseqReader(fname="D:/Users/Polina/3Dpredictor/input/Dactily/RNA-seq/dact_RNA",
                                       name="RNA")
    RNA.read_file(rename={"tr_ex": "gene",
                                          "start": "start",
                                          "end": "end",
                                          "chr": "chr",
                                          "FPKM": "sigVal"},
                                  sep="\t")
    print(RNA.chr_data["chr1"])
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

def test_WeightFunc():
    predictor = Predictor()
    data = predictor.read_file('out/2018-09-25-training.RandOnchr10oe.gz.1000000.50001.500000.25000.txt')
    contacts = np.array(data["contact_count"].values)
    header = predictor.get_avaliable_predictors('out/2018-09-25-training.RandOnchr10oe.gz.1000000.50001.500000.25000.txt')
    predictors = [h for h in header if not h in predictor.constant_nonpredictors]
    predictors_df = data[predictors]
    # predictors_df.columns.get_loc("IsLoop")
    # #print(predictors_df['IsLoop'] != 0)
    # idx_loop = np.flatnonzero(predictors_df['IsLoop'])
    # print(contacts[78])
    # print(predictors_df.loc[78, :])
    res = decorateContactWeither(contactWeitherFunction, power=3, coeff=100)
    result = res(contacts, predictors_df)
    print(result)

def test_RNAseqReader():
    RNA = RNAseqReader(fname="input/GM12878/RNA-seq/rna-seqPolyA.tsvpre.txt")
    RNA.read_file(rename={"Gene name":"gene",
                          "Gene start (bp)":"start",
                          "Gene end (bp)":"end",
                          "Chromosome/scaffold name":"chr",
                          "FPKM":"sigVal"},
                  sep="\t")
    logging.getLogger(__name__).info(iter(RNA.chr_data.values()).__next__().head())
    for inteval2 in [Interval("chr1",36511867,36528237),
                        Interval("chr1",36511866,36528237),
                        Interval("chr1",36511867, 36528238),
                        Interval("chr1",36511866, 36528238),
#                        Interval("chr1",36528238, 36528238),
#                        Interval("chr1",36528238, 36528239),
#                        Interval("chr1",36528138, 36528149),
#                        Interval("chr1",36500000, 36550000),
                        Interval("chr5", 7280119, 7345756),
                     Interval("chr5", 7311493, 7330491)]: #The last gives wrong results
        logging.info("------------------")
        logging.info(inteval2)
        logging.info(str(RNA.get_interval(inteval2))) #New func, based on Polina's intersect_intervals
        logging.info(str(RNA._get_interval(inteval2))) #Original Polina's intersect intervals
    print(RNA.chr_data['chr1'])
    RNA.delete_region(Interval("chr1", 11869, 29000))
    print(RNA.chr_data['chr1'])
def test_nonchip_reader():
    # ctcfreader = ChiPSeqReader("input/GM12878/")

    reader = ChiPSeqReader("input/GM12878/cage/GSM849368_hg19_wgEncodeRikenCageGm12878CellPapClusters.bed")
    reader.read_file(renamer={"0":"chr","1":"start","2":"end","4":"sigVal"})
    print(reader.chr_data['chr1'])
#test_nonchip_reader()

def test_add_loop():
    predictor = Predictor()
    validation_data= predictor.read_file("D:/Users/Polina/3Dpredictor/out/GM12878/2018-10-11-training.RandOnchr1oe.gz.12.1500000.50001.25000.25000.txt")
    add_loop(validation_data, "D:/Users/Polina/3Dpredictor/input/Loops/GM12878/GM12878.25000.loops")
    print(validation_data.query('IsLoop!=0'))
# test_add_loop()
#test_RNAseqReader()
#test_WeightFunc()
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
# test_ChipSeqRemoval()
test_ContactsRemoval() #TODO it doesn't throw errors, however the behaviour was not thoroughly tested
#test_read_orient()
# test_RNAseqReader()