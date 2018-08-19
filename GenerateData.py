import logging
from ChiPSeqReader import ChiPSeqReader
from Contacts_reader import ContactsReader
from E1_Reader import E1Reader,fileName2binsize
from shared import Interval

# Move with 1MB windows though chromosome
# For each window, generate binned CTCF track
# Next get all contacts and put them into file together with contacts count,
# CTCF track and coordinates relative to window


def contacts_to_file(contacts_sample,ctcf_reader,eig_reader,out_file_name,window_size,binsize):
    def contact_to_file(contact, file, ctcf_reader,eig_reader, binsize, window_size, N_fields): #write one contact to file
        window_start = max(0, contact.contact_st - ((window_size - contact.dist) // 2))
        window_end = window_start + window_size
        contacts_relative_start = contact.contact_st - window_start
        assert contacts_relative_start >= 0
        contacts_relative_end = contact.contact_en - window_start
        assert contacts_relative_end >= contacts_relative_start
        line = [contact.chr, contact.contact_st, contact.contact_en, window_start, window_end,
                contacts_relative_start, contacts_relative_end, contact["contact_count"]]
        interval = Interval(contact.chr, window_start, window_end)
        line += ctcf_reader.get_binned_interval(interval, binsize=binsize)
        line += eig_reader.get_E1inInterval(interval)["E1"].tolist()
        assert len(line) == N_fields
        file.write("\t".join(map(str, line)) + "\n")

    #loop over contacts, calculate predictors, write to file
    logging.info("Writing data to file "+out_file_name)
    file = open(out_file_name,"w")
    header = ["chr","contact_st","contact_en","window_start","window_end",
                "contacts_relative_start","contacts_relative_end","contact_count"]
    header += ["CTCF_bin"+str(i) for i in range(window_size // binsize)] #CTCF_values
    header += ["E1_bin" + str(i) for i in range(window_size // eig_reader.get_binsize())]  # E1_values
    file.write("\t".join(header)+"\n")
    contacts_sample.apply(contact_to_file,axis="columns",
                          file=file,ctcf_reader=ctcf_reader,eig_reader=eig_reader,
                          binsize=binsize,window_size=window_size,
                          N_fields = len(header))
    eig_reader.print_varnings()
    file.close()


logging.basicConfig(level=logging.DEBUG)

#constants
window_size = 500000
mindist = 50001
maxdist = window_size
binsize = 20000
sample_size = 500000

training_file_name = "training.RandOnChr1"+".".join(map(str,[window_size,mindist,maxdist,binsize,sample_size]))+".txt"
validation_file_name = "validating.38Mb_58MbOnChr2"+".".join(map(str,[window_size,mindist,maxdist,binsize,sample_size]))+".txt"
#input_folder = "D:/Lab Archive/ForChrRearrModel/"
input_folder = "input/"

#Read contacts data
contacts_reader = ContactsReader()
contacts_reader.read_files([input_folder + "chr1.5MB.Hepat.contacts",
                            input_folder + "chr2.5MB.Hepat.contacts"])

# Read CTCF data
ctcf_reader = ChiPSeqReader(input_folder + "Hepat_WT_MboI_rep1-rep2.IDR0.05.filt.narrowPeak")
ctcf_reader.read_file()

#Read E1 data
eig_reader = E1Reader()
eig_reader.read_files([input_folder + "chr1.Hepat.E1.50k",
                       input_folder + "chr2.Hepat.E1.50k"],
                      binSizeFromName=fileName2binsize)


def generate_training(sample_size):
    #Select random contacts on chromosome chr1 for trainings
    chrName = "chr1"
    region = Interval(chrName,
                      contacts_reader.get_min_contact_position(chrName),
                      contacts_reader.get_max_contact_position(chrName))
    contacts = contacts_reader.get_contacts(region,mindist=mindist,maxdist=maxdist)
    sample_size = min(sample_size,len(contacts))
    logging.info("Using sample size "+str(sample_size))
    contacts_sample = contacts.sample(n=sample_size)
    contacts_to_file(contacts_sample,ctcf_reader,eig_reader,
                     training_file_name,window_size,binsize)


def generate_test(sample_size):
    #Select random contacts on chromosome chr1 for trainings
    region = Interval("chr2",39800000,40000000)
    contacts = contacts_reader.get_contacts(region,mindist=mindist,maxdist=maxdist)
    sample_size = min(sample_size,len(contacts))
    logging.info("Using sample size "+str(sample_size))
    contacts_sample = contacts#.sample(n=sample_size)
    contacts_to_file(contacts_sample,ctcf_reader,eig_reader,
                     validation_file_name,window_size,binsize)


generate_test(sample_size)
#generate_training(sample_size)

from GenerateData2 import DataGenerator,E1PredictorGenerator,CTCFPredictorGenerator, \
                            SmallCTCFPredictorGenerator,SmallE1PredictorGenerator
generator = DataGenerator()
validation_file_name = validation_file_name + ".v1"
e1pg = E1PredictorGenerator(eig_reader,window_size)
ctcfpg = CTCFPredictorGenerator(ctcf_reader,binsize,window_size)
region = Interval("chr2", 39800000, 40000000)
contacts = contacts_reader.get_contacts(region, mindist=mindist, maxdist=maxdist)
print("----------",len(contacts))
generator.contacts2file(contacts,[ctcfpg,e1pg],validation_file_name+".v2")

window_size = 12500
#print(contacts.head())
e1pg_small = SmallE1PredictorGenerator(eig_reader,window_size)
ctcfpg_small = SmallCTCFPredictorGenerator(ctcf_reader,binsize,window_size,N_closest=5)
generator.contacts2file(contacts,[ctcfpg_small,e1pg_small],validation_file_name+".v3")

