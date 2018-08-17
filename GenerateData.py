import logging
from ChiPSeqReader import ChiPSeqReader
from Contacts_reader import ContactsReader
from shared import Interval

# Move with 1MB windows though chromosome
# For each window, generate binned CTCF track
# Next get all contacts and put them into file together with contacts count,
# CTCF track and coordinates relative to window


def contacts_to_file(contacts_sample,ctcf_reader,out_file_name,window_size,binsize,):
    def contact_to_file(contact, file, ctcf_reader, binsize, window_size, N_fields): #write one contact to file
        window_start = max(0, contact.contact_st - ((window_size - contact.dist) // 2))
        window_end = window_start + window_size
        contacts_relative_start = contact.contact_st - window_start
        assert contacts_relative_start >= 0
        contacts_relative_end = contact.contact_en - window_start
        assert contacts_relative_end >= contacts_relative_start
        line = [contact.chr, contact.contact_st, contact.contact_en, window_start, window_end,
                contacts_relative_start, contacts_relative_end, contact["contact_count"]]
        line += ctcf_reader.get_binned_interval(Interval(contact.chr, window_start, window_end), binsize=binsize)
        assert len(line) == N_fields
        file.write("\t".join(map(str, line)) + "\n")

    #loop over contacts, calculate predictors, write to file
    logging.info("Writing data to file "+out_file_name)
    file = open(out_file_name,"w")
    header = ["chr","contact_st","contact_en","window_start","window_end",
                "contacts_relative_start","contacts_relative_end","contact_count"] \
             + [str(i) for i in range(window_size//binsize)]
    file.write("\t".join(header)+"\n")
    contacts_sample.apply(contact_to_file,axis="columns",
                   file=file,ctcf_reader=ctcf_reader,binsize=binsize,window_size=window_size,
                          N_fields = len(header))
    file.close()


logging.basicConfig(level=logging.DEBUG)

#constants
window_size = 1000000
mindist = 50000
maxdist = window_size
binsize = 5000
sample_size = 500000

training_file_name = "training.RandOnChr1"+".".join(map(str,[window_size,mindist,maxdist,binsize,sample_size]))+".txt"
validation_file_name = "validating.38Mb_58MbOnChr2"+".".join(map(str,[window_size,mindist,maxdist,binsize,sample_size]))+".txt"


#Read contacts data
contacts_reader = ContactsReader()
contacts_reader.read_files(["C:/Users/FishmanVS/Desktop/RNF3D_beds/chr1.5MB.Hepat.contacts",
                            "C:/Users/FishmanVS/Desktop/RNF3D_beds/chr2.5MB.Hepat.contacts"])

# Read CTCF data
ctcf_reader = ChiPSeqReader("D:/Lab Archive/ForChrRearrModel/Hepat_WT_MboI_rep1-rep2.IDR0.05.filt.narrowPeak")
ctcf_reader.read_file()


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
    contacts_to_file(contacts_sample,ctcf_reader,training_file_name,window_size,binsize)


def generate_test(sample_size):
    #Select random contacts on chromosome chr1 for trainings
    region = Interval("chr2",38000000,58000000)
    contacts = contacts_reader.get_contacts(region,mindist=mindist,maxdist=maxdist)
    sample_size = min(sample_size,len(contacts))
    logging.info("Using sample size "+str(sample_size))
    contacts_sample = contacts.sample(n=sample_size)
    contacts_to_file(contacts_sample,ctcf_reader,validation_file_name,window_size,binsize)

generate_test(sample_size)
generate_training(sample_size)
