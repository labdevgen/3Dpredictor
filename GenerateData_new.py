import logging
from ChiPSeqReader import ChiPSeqReader
from Contacts_reader import ContactsReader
from E1_Reader import E1Reader,fileName2binsize
from shared import Interval, Parameters
from DataGenerators import DataGenerator,E1PredictorGenerator,CTCFPredictorGenerator, \
                            SmallCTCFPredictorGenerator,SmallE1PredictorGenerator

def generate_data(params):
    contacts = params.contacts_reader.get_contacts(params.interval,mindist=params.mindist,maxdist=params.maxdist)
    print( len(contacts) )
    sample_size = min(params.sample_size,len(contacts))
    logging.info("Using sample size "+str(sample_size))
    contacts_sample = contacts.sample(n=sample_size)
    assert len(contacts_sample) == sample_size
    generator = DataGenerator()
    generator.contacts2file(contacts_sample, params.pgs, params.out_file)

logging.basicConfig(level=logging.DEBUG)
input_folder = "D:/Lab Archive/ForChrRearrModel/input/"
#input_folder =  "input"

params = Parameters()
params.window_size = 25000 #region around contact to be binned for predictors
#params.small_window_size = 12500 #region  around contact ancors to be considered as cis
params.mindist = 50001 #minimum distance between contacting regions
#params.maxdist = params.window_size #max distance between contacting regions
params.maxdist = 3000000
params.binsize = 20000 #when binning regions with predictors, use this binsize
params.sample_size = 500000 #how many contacts write to file
params.conttype = "contacts"

training_file_name = "trainingSmall.RandOnChr1."+str(params)+".txt"
validation_file_name = "validatingSmall.noChr2."+str(params)+".txt"
logging.debug("Using input folder "+input_folder)

#Read contacts data
params.contacts_reader = ContactsReader()
params.contacts_reader.read_files([input_folder + "chr1.5MB.Hepat."+params.conttype,
                            input_folder + "chr2.5MB.Hepat."+params.conttype,
                            input_folder + "chr10.5MB.Hepat."+params.conttype])

# Read CTCF data
#params.ctcf_reader = ChiPSeqReader(input_folder + "Hepat_WT_MboI_rep1-rep2.IDR0.05.filt.narrowPeak")
params.ctcf_reader = ChiPSeqReader(input_folder + "Hepat_WT_MboI_rep1-rep2.IDR0.05.filt.narrowPeak_no_chr2")
params.ctcf_reader.read_file()

#Read E1 data
params.eig_reader = E1Reader()
params.eig_reader.read_files([input_folder + "chr1.Hepat.E1.50k",
                       input_folder + "chr2.Hepat.E1.50k",
                       input_folder + "chr10.Hepat.E1.50k"],
                      binSizeFromName=fileName2binsize) #infer size of E1 bins from file name using this function

#e1pg = E1PredictorGenerator(params.eig_reader,params.window_size)
#ctcfpg = CTCFPredictorGenerator(params.ctcf_reader,params.binsize,params.window_size)
#assert params.maxdist <= params.window_size #shouldn't be > window_size
#params.pgs = [e1pg,ctcfpg]

e1pg_small = SmallE1PredictorGenerator(params.eig_reader,params.window_size)
ctcfpg_small = SmallCTCFPredictorGenerator(params.ctcf_reader,params.binsize,params.window_size,N_closest=5)
params.pgs = [e1pg_small,ctcfpg_small]

#Generate train
trainChrName = "chr1"
params.interval = Interval(trainChrName,
                      params.contacts_reader.get_min_contact_position(trainChrName),
                      params.contacts_reader.get_max_contact_position(trainChrName))
params.out_file = training_file_name
#generate_data(params)

#Generate test
#params.interval = Interval("chr10", 59000000, 62000000)
params.interval = Interval("chr2", 47900000, 53900000)
params.out_file = params.interval.toFileName() + validation_file_name
generate_data(params)
