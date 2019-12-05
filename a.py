import logging
from Contacts_reader import ContactsReader
from shared import Parameters
from shared import Interval
import sys
import os

p=os.getcwd()+"\\nn\source"
sys.path.append(p)
from fastaFileReader import fastaReader, rm_chr_from_chrName
from SequencePredictorGenerator import SequencePredictorGenerator
from DataGenerator import generate_data

chr_num = "chr19"  # comma separated number of chromosomes for predictor generation
chr_nums = chr_num.split(",")
conttype = "contacts.gz"  # contacts.gz or oe.gz

# chr_num="12,13,14"
# conttype = "contacts.gz"
logging.basicConfig(format='%(asctime)s %(name)s: %(message)s', datefmt='%I:%M:%S', level=logging.DEBUG)
if __name__ == '__main__':  # Requiered for parallelization, at least on Windows
    for conttype in [conttype]:
        logging.basicConfig(format='%(asctime)s %(name)s: %(message)s', datefmt='%I:%M:%S', level=logging.DEBUG)
        input_folder =  os.getcwd() + "/input/chr19_mm10/"
        output_folder = os.getcwd() + "/output/chr19_mm10/"
        cell_type = "NPC"
        params = Parameters()
        params.window_size = 25000  # region around contact to be binned for predictors
        params.mindist = 50001  # minimum distance between contacting regions
        params.maxdist = 1500000  # maximum distance between contacting regions
        params.sample_size = 1  # how many contacts write to file
        params.conttype = conttype
        params.max_cpus = 11
        params.keep_only_orient = False  # set True if you want use only CTCF with orient
        params.use_only_contacts_with_CTCF = "all_cont"  # "cont_with_CTCF"  #this option use for training to change proportion
        # of contacts with nearest ctcf sites
        write_all_chrms_in_file = False  # set True if you have train with few chromosomes. Need for writing different chromosomes in the same file

        fill_empty_contacts = False
        logging.getLogger(__name__).debug("Using input folder " + input_folder)

        # Read contacts data
        params.contacts_reader = ContactsReader()
        contacts_files = [input_folder + "19.contacts.gz"]
        coeff_fname = input_folder + "coefficient.NPC.5000.txt"
        # set path to the coefficient file and to contacts files
        # contacts file format: bin_start--bin_end--contact_count
        params.contacts_reader.read_files(contacts_files, coeff_fname,
                                          max_cpus=params.max_cpus,
                                          fill_empty_contacts=fill_empty_contacts, maxdist=params.maxdist)
        params.fastaReader = fastaReader(input_folder + "chr19.fa",chrm_names_renamer = rm_chr_from_chrName)
        params.fastaReader.read_data()
        SequencePG = SequencePredictorGenerator(fastaReader=params.fastaReader, binsize=params.contacts_reader.binsize)
        params.pgs = [SequencePG]
        params.out_file = output_folder + "NPC_5000"

        params.sample_size = 100
        params.interval = Interval("19", params.contacts_reader.get_min_contact_position("19"),
                                       params.contacts_reader.get_max_contact_position("19"))
        logging.getLogger(__name__).info("Generating dataset for interval " + str(params.interval))
        generate_data(params)