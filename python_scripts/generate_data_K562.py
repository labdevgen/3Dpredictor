import os
import sys
# add source directory into path to allow import
sourcedir = os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(sys.argv[0]))),
                         "source")
sys.path.append(sourcedir)
import logging
from ChiPSeqReader import ChiPSeqReader
from Contacts_reader import ContactsReader
from RNASeqReader import RNAseqReader
from TssReader import TssReader
from E1_Reader import E1Reader,fileName2binsize
from shared import Interval, Parameters
from DataGenerator import generate_data
from PredictorGenerators import E1PredictorGenerator,ChipSeqPredictorGenerator, \
                SmallChipSeqPredictorGenerator,SmallE1PredictorGenerator, \
                SitesOrientPredictorGenerator, OrientBlocksPredictorGenerator, ConvergentPairPredictorGenerator, Distance_to_TSS_PG
from VectPredictorGenerators import loopsPredictorGenerator
from LoopReader import LoopReader
import pandas as pd
import sys
import pickle
chr_num=sys.argv[1] #comma separated number of chromosomes for predictor generation
chr_nums=chr_num.split(",")
conttype = sys.argv[2] #contacts.gz or oe.gz

# chr_num="12,13,14"
# conttype = "contacts.gz"

if __name__ == '__main__': #Requiered for parallelization, at least on Windows
    for conttype in [conttype]:
        logging.basicConfig(format='%(asctime)s %(name)s: %(message)s', datefmt='%I:%M:%S', level=logging.DEBUG)
        input_folder ="/input/K562/"
        output_folder = "/output/K562/"
        cell_type="K562"
        params = Parameters()
        params.window_size = 5000 #region around contact to be binned for predictors
        params.mindist = 10001 #minimum distance between contacting regions
        params.maxdist = 1500000 #maximum distance between contacting regions
        params.sample_size = 30000 #how many contacts write to file
        params.conttype = conttype
        params.max_cpus = 11
        params.keep_only_orient=False #set True if you want use only CTCF with orient
        params.use_only_contacts_with_CTCF = "all_cont"#"cont_with_CTCF"  #this option use for training to change proportion
                                                                          #of contacts with nearest ctcf sites
        write_all_chrms_in_file=True #set True if you have train with few chromosomes. Need for writing different chromosomes in the same file

        fill_empty_contacts=False
        logging.getLogger(__name__).debug("Using input folder "+input_folder)

        #Read contacts data
        params.contacts_reader = ContactsReader()
        contacts_files = []
        #set path to the coefficient file and to contacts files
        #contacts file format: bin_start--bin_end--contact_count
        [contacts_files.append(input_folder+ "chr"+chr+".5000."+params.conttype) for chr in chr_nums]
        params.contacts_reader.read_files(contacts_files, coeff_fname="coefficient."+cell_type+".5KB.txt", max_cpus=params.max_cpus,
                                          fill_empty_contacts=fill_empty_contacts, maxdist=params.maxdist)

        if params.use_only_contacts_with_CTCF == "cont_with_CTCF":
            params.proportion = 1 #propotion of contacts with ctcf in train data (if 1: contacts with ctcf/random contacts = 1/1)
            #set path to
            params.contacts_reader.use_contacts_with_CTCF(CTCFfile=input_folder + "CTCF/wgEncodeAwgTfbsHaibK562CtcfcPcr1xUniPk.narrowPeak.gz",
                                                          maxdist=params.maxdist, proportion=params.proportion, keep_only_orient=params.keep_only_orient,
                                                          CTCForientfile=input_folder + "CTCF/wgEncodeAwgTfbsHaibK562CtcfcPcr1xUniPk.narrowPeak-orient.bed")
            params.use_only_contacts_with_CTCF += str(params.contacts_reader.conts_with_ctcf)
        # params.contacts_reader.delete_region(Interval("chr22", 16064000, 16075000)) #example of rearrangement

        # Read CTCF data
        # CTCF_file format: ENCODE narrow peak
        # CTCF_orient_file format: chr--start--end--name--score--strand
        logging.info('create CTCF_PG')
        # set path to the CTCF chip-seq file:
        params.ctcf_reader = ChiPSeqReader(input_folder + "CTCF/wgEncodeAwgTfbsHaibK562CtcfcPcr1xUniPk.narrowPeak.gz",
                                                            name="CTCF")
        params.ctcf_reader.read_file()
        # set path to the CTCF_orient file:
        params.ctcf_reader.set_sites_orientation(
            input_folder + "CTCF/wgEncodeAwgTfbsHaibK562CtcfcPcr1xUniPk.narrowPeak-orient.bed")
        if params.keep_only_orient:
            params.ctcf_reader.keep_only_with_orient_data()
        # set corresponding predictor generators and its options:
        OrientCtcfpg = SitesOrientPredictorGenerator(params.ctcf_reader,
                                                     N_closest=4)
        NotOrientCTCFpg = SmallChipSeqPredictorGenerator(params.ctcf_reader,
                                                         params.window_size,
                                                         N_closest=4)

        # Read CTCF data and drop sites w/o known orientation
        params.ctcf_reader_orientOnly = ChiPSeqReader(input_folder + "CTCF/wgEncodeAwgTfbsHaibK562CtcfcPcr1xUniPk.narrowPeak.gz",
                                                            name="CTCF")
        params.ctcf_reader_orientOnly.read_file()
        params.ctcf_reader_orientOnly.set_sites_orientation(
            input_folder + "CTCF/wgEncodeAwgTfbsHaibK562CtcfcPcr1xUniPk.narrowPeak-orient.bed")
        params.ctcf_reader_orientOnly.keep_only_with_orient_data()
        # set corresponding predictor generators and its options:
        OrientBlocksCTCFpg = OrientBlocksPredictorGenerator(params.ctcf_reader_orientOnly,
                                                             params.window_size)
        ConvergentPairPG = ConvergentPairPredictorGenerator(params.ctcf_reader, binsize=params.window_size)

        #Read other chip-seq data. If you have a lot of chip-seq data you can use file with filenames (example in the input folder)
        logging.info('create chipPG')
        chipPG = []
        #filenames format: protein_name--file_name--database--id
        filenames_df = pd.read_csv(input_folder + "peaks/filenames.csv")

        for index, row in filenames_df.iterrows():
                params.chip_reader = ChiPSeqReader(input_folder + 'peaks/' + row["filename"] + '.gz', name=row['name'])
                params.chip_reader.read_file()
                # set corresponding predictor generators and its options:
                chipPG.append(SmallChipSeqPredictorGenerator(params.chip_reader,params.window_size,N_closest=4))
        # # #
        # # Read methylation data
        # logging.info('create metPG')
        # metPG = []
        # filemanes_df = pd.read_csv(input_folder + "methylation/filenames.csv")
        # # assert len(os.listdir(input_folder + 'peaks/')) - 1 == len(filenames_df['name'])
        # for index, row in filemanes_df.iterrows():
        #     #print(row["name"])
        #     params.met_reader = ChiPSeqReader(input_folder + 'methylation/'+ row["filename"], name=row['name'])
        #     params.met_reader.read_file(renamer={"0":"chr","1":"start","2":"end","4":"sigVal"})
        #     metPG.append(SmallChipSeqPredictorGenerator(params.met_reader,params.window_size,N_closest=4))
        #
        # #Read cage data
        # cagePG = []
        # filemanes_df = pd.read_csv(input_folder + "cage/filenames.csv")
        # assert len(os.listdir(input_folder + 'cage/')) - 1 == len(filemanes_df['name'])
        # for index, row in filemanes_df.iterrows():
        #     params.cage_reader = ChiPSeqReader(input_folder+"cage/GSM849365_hg19_wgEncodeRikenCageK562CellPapClusters.bed.gz", name=row['name'])# + "cage/" + row["filename"], name=row['name'])
        #     params.cage_reader.read_file(renamer={"0":"chr","1":"start","2":"end","4":"sigVal"})
        #     cagePG.append(SmallChipSeqPredictorGenerator(params.cage_reader,params.window_size,N_closest=4))

        #Read RNA-Seq data
        #RNA-seq_file format: this file should have fields "gene", "start", "end", "chr","FPKM"
        #you can rename table fields below
        params.RNAseqReader = RNAseqReader(fname=input_folder + "RNA-seq/rna-seqPolyA.tsvpre.txt",
                                           name="RNA")
        #read RNA-seq data and rename table fields
        params.RNAseqReader.read_file(rename={ "Gene name": "gene",
                              "Gene start (bp)": "start",
                              "Gene end (bp)": "end",
                              "Chromosome/scaffold name": "chr",
                              "FPKM": "sigVal"},
                      sep="\t")
        # set corresponding predictor generators and its options:
        RNAseqPG = SmallChipSeqPredictorGenerator(params.RNAseqReader,
                                                  window_size=params.window_size,
                                                  N_closest=3)
        #Read TSS data
        #TSS_file format: "0":"chr","1":"start","2":"end","5":"strand", "6":"TSS_start", "7":"TSS_end". you can rename on your own below
        params.TssReader=TssReader(fname=input_folder + "TSS/NCBI_refSeq_hg19.bed", name="TSS")
        params.TssReader.read_file()
        # set corresponding predictor generators and its options:
        TSSPG=Distance_to_TSS_PG(params.TssReader)


        # #Read E1 data
        # params.eig_reader = E1Reader()
        # e1_files = [input_folder + "E1/"+ "chr" + str(i) + ".GM12878.25K.txt" for i in range(1, 23)]
        # e1_files.append(input_folder +"E1/" + "chrX.GM12878.25K.txt")
        # params.eig_reader.read_files(e1_files, binSizeFromName=fileName2binsize) #infer size of E1 bins from file name using this function
        # e1pg = SmallE1PredictorGenerator(params.eig_reader,params.window_size)

        # write all predictor generators which you want to use:
        params.pgs = [OrientCtcfpg, NotOrientCTCFpg, OrientBlocksCTCFpg ,RNAseqPG, ConvergentPairPG,TSSPG]+chipPG#+cagePG+metPG+chipPG

        # Generate train
        train_chrs=[]
        [train_chrs.append("chr"+chr) for chr in chr_nums]
        if write_all_chrms_in_file:
            train_file_name="training.RandOn"+ str(params)
            params.out_file=output_folder+"_".join(train_chrs)+train_file_name
        for trainChrName in train_chrs:
            training_file_name = "training.RandOn" + trainChrName + str(params) + ".txt"
            # set it if you want to use all contacts of chromosome for training:
            params.sample_size = len(params.contacts_reader.data[trainChrName])
            # if you want to use only an interval of chromosome, set its coordinates:
            params.interval = Interval(trainChrName,
                                  params.contacts_reader.get_min_contact_position(trainChrName),
                                  params.contacts_reader.get_max_contact_position(trainChrName))
            
            if not write_all_chrms_in_file:
                train_file_name = "training.RandOn" + str(params) + ".txt"
                params.out_file = output_folder + params.interval.toFileName() + train_file_name
            generate_data(params,saveFileDescription=True)
            if not write_all_chrms_in_file:
                del(params.out_file)
            del (params.sample_size)


        # Generate test
        validate_chrs=[]
        [validate_chrs.append("chr"+chr) for chr in chr_nums]
        if write_all_chrms_in_file:
            validation_file_name = "validatingOrient." + str(params) + ".txt"
            params.out_file = output_folder + "_".join(validate_chrs) + validation_file_name
        for validateChrName in validate_chrs:

            params.sample_size = len(params.contacts_reader.data[validateChrName])
            # set another interval if you want to predict contacts for specific chromosome interval:
            params.interval = Interval(validateChrName,
                                       params.contacts_reader.get_min_contact_position(validateChrName),
                                       params.contacts_reader.get_max_contact_position(validateChrName))
                              #Interval("chr2", 118000000, 129000000)
            logging.getLogger(__name__).info("Generating validation dataset for interval "+str(params.interval))
            if not write_all_chrms_in_file:
                validation_file_name = "validatingOrient." + str(params) + ".txt"
                params.out_file = output_folder + params.interval.toFileName() + validation_file_name
            generate_data(params)
            if not write_all_chrms_in_file:
                del(params.out_file)
            del (params.sample_size)