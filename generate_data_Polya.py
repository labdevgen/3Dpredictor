# import sys
# sys.path = ["/mnt/storage/home/vsfishman/.local/lib/python3.5/site-packages"] + sys.path
import logging
from ChiPSeqReader import ChiPSeqReader
from Contacts_reader import ContactsReader
from RNASeqReader import RNAseqReader
from E1_Reader import E1Reader,fileName2binsize
from shared import Interval, Parameters
from DataGenerator import generate_data
from PredictorGenerators import E1PredictorGenerator,ChipSeqPredictorGenerator, \
                SmallChipSeqPredictorGenerator,SmallE1PredictorGenerator, \
                SitesOrientPredictorGenerator, OrientBlocksPredictorGenerator
from VectPredictorGenerators import loopsPredictorGenerator
from LoopReader import LoopReader
import pandas as pd
import os


if __name__ == '__main__': #Requered for parallization, at least on Windows
    #,"chr10", "chr1"]:
    for conttype in ["contacts.gz"]:
        logging.basicConfig(format='%(asctime)s %(name)s: %(message)s', datefmt='%I:%M:%S', level=logging.DEBUG)

        input_folder ="input/GM12878/"
        #output_folder = "D:/Users/Polina/3Dpredictor/"
        output_folder = "out/GM12878/validating_chrms_2/"
        cell_type="GM12878"
        #input_folder =  "input"

        params = Parameters()
        params.window_size = 25000 #region around contact to be binned for predictors
        #params.small_window_size = 12500 #region  around contact ancors to be considered as cis
        params.mindist = 50001 #minimum distance between contacting regions
        #params.maxdist = params.window_size #max distance between contacting regions
        params.maxdist = 1500000
        #params.binsize = 20000 #when binning regions with predictors, use this binsize
        params.sample_size = 25000 #how many contacts write to file
        #params.conttype = "oe.gz"
        params.conttype = conttype
        params.max_cpus = 8

        logging.getLogger(__name__).debug("Using input folder "+input_folder)

        #Read contacts data
        params.contacts_reader = ContactsReader()
        contacts_files = []
        contacts_files=[input_folder+ "chr"+str(i)+".5MB.GM12878."+params.conttype for i in (15,16)]
        contacts_files.append(input_folder+ "chrX.5MB.GM12878."+params.conttype)
        params.contacts_reader.read_files(contacts_files, coeff_fname="coefficient."+cell_type+".txt")
        #
        # #Loops predictor
        # loopsReader = LoopReader("input/Loops/GM12878/GM12878.25000.loops")
        # loopsReader.read_loops()
        # loopspg = loopsPredictorGenerator(loopsReader, params.window_size)

        # Read CTCF data
        logging.info('create CTCF_PG')
        params.ctcf_reader = ChiPSeqReader(input_folder + "CTCF/wgEncodeAwgTfbsBroadGm12878CtcfUniPk.narrowPeak.gz",
                                                            name="CTCF")
        params.ctcf_reader.read_file()
        params.ctcf_reader.set_sites_orientation(
            input_folder + "CTCF/wgEncodeAwgTfbsBroadGm12878CtcfUniPk.narrowPeak-orient.bed")


        OrientCtcfpg = SitesOrientPredictorGenerator(params.ctcf_reader,
                                                     N_closest=4)
        NotOrientCTCFpg = SmallChipSeqPredictorGenerator(params.ctcf_reader,
                                                         params.window_size,
                                                         N_closest=4)

        # Read CTCF data and drop sites w/o known orientation
        params.ctcf_reader_orientOnly = ChiPSeqReader(input_folder + "CTCF/wgEncodeAwgTfbsBroadGm12878CtcfUniPk.narrowPeak.gz",
                                                            name="CTCF")
        params.ctcf_reader_orientOnly.read_file()
        params.ctcf_reader_orientOnly.set_sites_orientation(
            input_folder + "CTCF/wgEncodeAwgTfbsBroadGm12878CtcfUniPk.narrowPeak-orient.bed")
        params.ctcf_reader_orientOnly.keep_only_with_orient_data()
        OrientBlocksCTCFpg = OrientBlocksPredictorGenerator(params.ctcf_reader_orientOnly,
                                                             params.window_size)

        # #Read other chip-seq data
        # logging.info('create chipPG')
        # chipPG = []
        # filenames_df = pd.read_csv(input_folder + "peaks/filenames.csv")
        # assert len(os.listdir(input_folder + 'peaks/')) - 1 == len(filenames_df['name'])
        # # print(len(os.listdir(input_folder + 'peaks/')))
        # # print(len(filenames_df['name']))
        # for index, row in filenames_df.iterrows():
        #     params.chip_reader = ChiPSeqReader(input_folder + 'peaks/' + row["filename"] + '.gz', name=row['name'])
        #     params.chip_reader.read_file()
        #     chipPG.append(SmallChipSeqPredictorGenerator(params.chip_reader,params.window_size,N_closest=4))
        # # #
        # #Read methylation data
        # logging.info('create metPG')
        # metPG = []
        # filemanes_df = pd.read_csv(input_folder + "methylation/filenames.csv")
        # assert len(os.listdir(input_folder + 'peaks/')) - 1 == len(filenames_df['name'])
        # for index, row in filemanes_df.iterrows():
        #     #print(row["name"])
        #     params.met_reader = ChiPSeqReader(input_folder + 'methylation/'+ row["filename"], name=row['name'])
        #     params.met_reader.read_file(renamer={"0":"chr","1":"start","2":"end","4":"sigVal"})
        #     metPG.append(SmallChipSeqPredictorGenerator(params.met_reader,params.window_size,N_closest=4))
        #Read cage data
        cagePG = []
        filemanes_df = pd.read_csv(input_folder + "cage/filenames.csv")
        # assert len(os.listdir(input_folder + 'cage/')) - 1 == len(filemanes_df['name'])
        for index, row in filemanes_df.iterrows():
            #print(row["name"])
            params.cage_reader = ChiPSeqReader(input_folder + "cage/" + row["filename"], name=row['name'])
            params.cage_reader.read_file(renamer={"0":"chr","1":"start","2":"end","4":"sigVal"})
            cagePG.append(SmallChipSeqPredictorGenerator(params.cage_reader,params.window_size,N_closest=4))
        #Read RNA-Seq data
        params.RNAseqReader = RNAseqReader(fname=input_folder + "RNA-seq/rna-seqPolyA.tsvpre.txt",
                                           name="RNA")
        params.RNAseqReader.read_file(rename={"Gene name": "gene",
                              "Gene start (bp)": "start",
                              "Gene end (bp)": "end",
                              "Chromosome/scaffold name": "chr",
                              "FPKM": "sigVal"},
                      sep="\t")
        RNAseqPG = SmallChipSeqPredictorGenerator(params.RNAseqReader,
                                                  window_size=params.window_size,
                                                  N_closest=3)

        # #Read E1 data
        # params.eig_reader = E1Reader()
        # e1_files = [input_folder + "E1/"+ "chr" + str(i) + ".GM12878.25K.txt" for i in range(1, 23)]
        # e1_files.append(input_folder +"E1/" + "chrX.GM12878.25K.txt")
        # params.eig_reader.read_files(e1_files, binSizeFromName=fileName2binsize) #infer size of E1 bins from file name using this function
        #
        # e1pg = SmallE1PredictorGenerator(params.eig_reader,params.window_size)
        #
        params.pgs = [OrientCtcfpg, NotOrientCTCFpg, OrientBlocksCTCFpg, RNAseqPG]+cagePG#+metPG+chipPG

        # # Generate train
        # for trainChrName in ["chr1"]:
        #     training_file_name = "2018-10-11-training.RandOn" + trainChrName + str(params) + ".txt"
        #     params.interval = Interval(trainChrName,
        #                           params.contacts_reader.get_min_contact_position(trainChrName),
        #                           params.contacts_reader.get_max_contact_position(trainChrName))
        #     params.out_file = output_folder + training_file_name
        #     generate_data(params,saveFileDescription=True)
        #     del(params.out_file)

        # Generate test
        validate_chrs=["chr20", "chr21","chr22", "chrX"]
        #validate_chrs.append("chrX")
        for validateChrName in validate_chrs:
            params.sample_size = len(params.contacts_reader.data[validateChrName])
            #print(params.sample_size)
            validation_file_name = "validatingOrient." + str(params) + ".txt"
            params.interval = Interval(validateChrName,
                                       params.contacts_reader.get_min_contact_position(validateChrName),
                                       params.contacts_reader.get_max_contact_position(validateChrName))
            logging.getLogger(__name__).info("Generating validation dataset for interval "+str(params.interval))
            params.out_file = output_folder + params.interval.toFileName() + validation_file_name
            generate_data(params)
            del(params.out_file)
            del (params.sample_size)

        # for interval in [Interval("chr5", 11500000, 17500000)]:
        # #                  Interval("chr10", 47900000, 53900000),
        # #                  Interval("chr10", 15000000, 20000000),
        # #                  Interval("chr10",36000000,41000000)]:
        # # Interval("chr1", 100000000, 110000000)]:
        #    logging.getLogger(__name__).info("Generating validation dataset for interval "+str(interval))
        #    validation_file_name = "validatingOrient." + str(params) + ".txt"
        #    params.interval = interval
        #    params.out_file = output_folder + params.interval.toFileName() + validation_file_name
        #    generate_data(params)

        # for object in [params.contacts_reader]+params.pgs:
        #     lostInterval = Interval("chr1",103842568,104979840)
        #     object.delete_region(lostInterval)
        #     params.interval = Interval("chr1",100000000,109000000)
        #     logging.getLogger(__name__).info("Saving data to file "+params.interval.toFileName() + "DEL." + lostInterval.toFileName()+validation_file_name)
        # params.out_file = params.interval.toFileName() + "DEL." + lostInterval.toFileName()+validation_file_name
        # generate_data(params)