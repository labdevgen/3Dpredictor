import logging
from ChiPSeqReader import ChiPSeqReader
from Contacts_reader import ContactsReader
from RNASeqReader import RNAseqReader
from E1_Reader import E1Reader,fileName2binsize
from shared import Interval, Parameters
from DataGenerator import generate_data
from PredictorGenerators import E1PredictorGenerator,ChipSeqPredictorGenerator, \
                            SmallChipSeqPredictorGenerator,SmallE1PredictorGenerator, SitesOrientPredictorGenerator, OrientBlocksPredictorGenerator, \
                            SitesOnlyOrientPredictorGenerator


logging.basicConfig(format='%(asctime)s %(name)s: %(message)s', datefmt='%I:%M:%S', level=logging.DEBUG)

input_folder ="input/"
#output_folder = "D:/Users/Polina/3Dpredictor/"
output_folder = "out/"
#input_folder =  "input"

params = Parameters()
params.window_size = 25000 #region around contact to be binned for predictors
#params.small_window_size = 12500 #region  around contact ancors to be considered as cis
params.mindist = 50001 #minimum distance between contacting regions
#params.maxdist = params.window_size #max distance between contacting regions
params.maxdist = 1000000
#params.binsize = 20000 #when binning regions with predictors, use this binsize
params.sample_size = 500000 #how many contacts write to file
params.conttype = "contacts.gz"

training_file_name = "2018-09-23-trainingOrient.RandOnChr1."+str(params)+".txt"
validation_file_name = "validatingOrient."+str(params)+".txt"
logging.getLogger(__name__).debug("Using input folder "+input_folder)

#Read contacts data
params.contacts_reader = ContactsReader()
params.contacts_reader.read_files([input_folder + "chr1.5MB.Hepat."+params.conttype,
                           input_folder + "chr2.5MB.Hepat."+params.conttype])
                            #input_folder + "chr10.5MB.Hepat."+params.conttype])
                            #input_folder + "chr6.5MB.Hepat." + params.conttype])

# Read CTCF data
params.ctcf_reader = ChiPSeqReader(input_folder + "Hepat_WT_MboI_rep1-rep2.IDR0.05.filt.narrowPeak",
                                                    name="CTCF")
params.ctcf_reader.read_file()
params.ctcf_reader.set_sites_orientation(
    input_folder + "Hepat_WT_MboI_rep1-rep2_IDR0_05_filt_narrowPeak-orient_N10.bed.gz")


OrientCtcfpg = SitesOrientPredictorGenerator(params.ctcf_reader,
                                             N_closest=4)
NotOrientCTCFpg = SmallChipSeqPredictorGenerator(params.ctcf_reader,
                                                 params.window_size,
                                                 N_closest=4)

# Read CTCF data and drop sites w/o known orientation
#params.ctcf_reader_orintOnly = ChiPSeqReader(input_folder + "Hepat_WT_MboI_rep1-rep2.IDR0.05.filt.narrowPeak",
#                                                    name="CTCF")
#params.ctcf_reader_orintOnly.read_file()
#params.ctcf_reader_orintOnly.set_sites_orientation(
#    input_folder + "Hepat_WT_MboI_rep1-rep2_IDR0_05_filt_narrowPeak-orient_N10.bed.gz")
#params.ctcf_reader_orintOnly.keep_only_with_orient_data()
#onlyOrientCtcfpg = SitesOnlyOrientPredictorGenerator(params.ctcf_reader_orintOnly,
#                                                     N_closest=3)

#Read RNA-Seq data
params.RNAseqReader = RNAseqReader(fname="input/GSE95111_genes.fpkm_table.txt.pre.txt",
                                   name="RNA")
params.RNAseqReader.read_file(rename={"Gene name": "gene",
                      "Gene start (bp)": "start",
                      "Gene end (bp)": "end",
                      "Chromosome/scaffold name": "chr",
                      "shCtrl-1_0": "sigVal"},
              sep="\t")
RNAseqPG = SmallChipSeqPredictorGenerator(params.RNAseqReader,
                                          window_size=params.window_size,
                                          N_closest=3)

#Read E1 data
params.eig_reader = E1Reader()
params.eig_reader.read_files([input_folder + "chr1.Hepat.E1.50k",
                       input_folder + "chr2.Hepat.E1.50k"],
                       #input_folder + "chr10.Hepat.E1.50k"],
                       #input_folder + "chr6.Hepat.E1.50k"],
                      binSizeFromName=fileName2binsize) #infer size of E1 bins from file name using this function

e1pg = SmallE1PredictorGenerator(params.eig_reader,params.window_size)

params.pgs = [e1pg,OrientCtcfpg,NotOrientCTCFpg,RNAseqPG]#,onlyOrientCtcfpg]

#Generate train
trainChrName = "chr1"
params.interval = Interval(trainChrName,
                      params.contacts_reader.get_min_contact_position(trainChrName),
                      params.contacts_reader.get_max_contact_position(trainChrName))
params.out_file = output_folder + training_file_name
generate_data(params,saveFileDescription=True)

#Generate test
for interval in [# Interval("chr10", 59000000, 62000000)]:
                  Interval("chr2", 47900000, 53900000),
                  Interval("chr2", 85000000, 92500000),
                  Interval("chr2",36000000,41000000)]:
                 # Interval("chr1", 100000000, 110000000)]:
    logging.getLogger(__name__).info("Generating validation dataset for interval "+str(interval))
    params.interval = interval
    params.out_file = output_folder + params.interval.toFileName() + validation_file_name
    generate_data(params)

# for object in [params.contacts_reader]+params.pgs:
#     lostInterval = Interval("chr1",103842568,104979840)
#     object.delete_region(lostInterval)
#     params.interval = Interval("chr1",100000000,109000000)
    #logging.getLogger(__name__).info("Saving data to file "+params.interval.toFileName() + "DEL." + lostInterval.toFileName()+validation_file_name)
# params.out_file = params.interval.toFileName() + "DEL." + lostInterval.toFileName()+validation_file_name
#generate_data(params)