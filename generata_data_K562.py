# import sys
# sys.path = ["/mnt/storage/home/vsfishman/.local/lib/python3.5/site-packages"] + sys.path
import logging
from ChiPSeqReader import ChiPSeqReader
from Contacts_reader import ContactsReader
from RNASeqReader import RNAseqReader
from TssReader import TssReader
from E1_Reader import E1Reader,fileName2binsize
from shared import Interval, Parameters
from DataGenerator import generate_data
from PredictoGenerators_new_edition import E1PredictorGenerator,ChipSeqPredictorGenerator, \
                SmallChipSeqPredictorGenerator,SmallE1PredictorGenerator, \
                SitesOrientPredictorGenerator, OrientBlocksPredictorGenerator, ConvergentPairPredictorGenerator, Distance_to_TSS_PG
from VectPredictorGenerators import loopsPredictorGenerator
from LoopReader import LoopReader
import pandas as pd
import os
import sys
import pickle
chr_num=sys.argv[1]
chr_nums=chr_num.split(",")
conttype = sys.argv[2]

# chr_num="12"
# conttype = "contacts.gz"

if __name__ == '__main__': #Requered for parallization, at least on Windows
    #,"chr10", "chr1"]:
    for conttype in [conttype]:
        print("hello")
        logging.basicConfig(format='%(asctime)s %(name)s: %(message)s', datefmt='%I:%M:%S', level=logging.DEBUG)

        input_folder ="/mnt/scratch/ws/psbelokopytova/201905031108polinaB/3DPredictor/input/K562/"
        #output_folder = "D:/Users/Polina/3Dpredictor/"
        output_folder = "/mnt/scratch/ws/psbelokopytova/201905031108polinaB/3DPredictor/out/K562/5KB/new/"
        cell_type="K562"
        lengths_dict = {'chr1': 1494930, 'chr3': 609806, 'chr5': 518646, 'chr7': 682860, 'chr11': 726290, 'chr13': 115324}
        params = Parameters()
        params.window_size = 5000 #region around contact to be binned for predictors
        params.mindist = 10001 #minimum distance between contacting regions
        params.maxdist = 1500000
        params.sample_size = 691309 #how many contacts write to file
        params.conttype = conttype
        params.max_cpus = 11
        params.keep_only_orient=False
        params.use_only_contacts_with_CTCF = "all_cont"#"cont_with_CTCF"#"all_cont"#"cont_with_CTCF "

        write_all_chrms_in_file=True
        fill_empty_contacts = False


        logging.getLogger(__name__).debug("Using input folder "+input_folder)

        #Read contacts data
        params.contacts_reader = ContactsReader()
        contacts_files = []
        # contacts_files=[input_folder+ "5KBcontacts/chr"+chr_num+".5KB."+params.conttype ]
        [contacts_files.append(input_folder+ "5KBcontacts/chr"+chr+".5KB."+params.conttype) for chr in chr_nums]
        params.contacts_reader.read_files(contacts_files, coeff_fname="coefficient."+cell_type+".5KB.txt", max_cpus=params.max_cpus,
                                          fill_empty_contacts=fill_empty_contacts, maxdist=params.maxdist)

        if params.use_only_contacts_with_CTCF == "cont_with_CTCF":
            params.proportion = 1
            params.contacts_reader.use_contacts_with_CTCF(CTCFfile=input_folder + "CTCF/wgEncodeAwgTfbsHaibK562CtcfcPcr1xUniPk.narrowPeak.gz",
                                                          maxdist=params.maxdist, proportion=params.proportion, keep_only_orient=params.keep_only_orient,
                                                          CTCForientfile=input_folder + "CTCF/wgEncodeAwgTfbsHaibK562CtcfcPcr1xUniPk.narrowPeak-orient.bed")
            params.use_only_contacts_with_CTCF += str(params.contacts_reader.conts_with_ctcf)
        # params.contacts_reader.duplicate_region(Interval("chr22", 16064000, 16075000))



        #Loops predictor
        # loopsReader = LoopReader("input/Loops/GM12878/GM12878.25000.loops")
        # loopsReader.read_loops()
        # loopspg = loopsPredictorGenerator(loopsReader, params.window_size)

        # Read CTCF data
        logging.info('create CTCF_PG')
        params.ctcf_reader = ChiPSeqReader(input_folder + "CTCF/wgEncodeAwgTfbsHaibK562CtcfcPcr1xUniPk.narrowPeak.gz",
                                                            name="CTCF")
        params.ctcf_reader.read_file()
        params.ctcf_reader.set_sites_orientation(
            input_folder + "CTCF/wgEncodeAwgTfbsHaibK562CtcfcPcr1xUniPk.narrowPeak-orient.bed")
        if params.keep_only_orient:
            params.ctcf_reader.keep_only_with_orient_data()

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
        OrientBlocksCTCFpg = OrientBlocksPredictorGenerator(params.ctcf_reader_orientOnly,
                                                             params.window_size)
        ConvergentPairPG = ConvergentPairPredictorGenerator(params.ctcf_reader, binsize=params.window_size)

        # #Read other chip-seq data
        logging.info('create chipPG')
        chipPG = []
        filenames_df = pd.read_csv(input_folder + "peaks/filenames.csv")
        # assert len(os.listdir(input_folder + 'peaks/')) - 1 == len(filenames_df['name'])
        # print(len(os.listdir(input_folder + 'peaks/')))
        # print(len(filenames_df['name']))
        proteins=set(["RAD21", "SMC3", "POLR2A", "H3K27ac", "H3K27me3", "DNase-seq", "H3K9me3", "H3K4me1", "H3K4me2", "H3K4me3", "YY1"])
        for index, row in filenames_df.iterrows():
            if row["name"] in proteins:
                params.chip_reader = ChiPSeqReader(input_folder + 'peaks/' + row["filename"] + '.gz', name=row['name'])
                params.chip_reader.read_file()
                chipPG.append(SmallChipSeqPredictorGenerator(params.chip_reader,params.window_size,N_closest=4))
        assert len(chipPG)==len(proteins)
        # # # #
        # # # #Read methylation data
        # # logging.info('create metPG')
        # # metPG = []
        # # filemanes_df = pd.read_csv(input_folder + "methylation/filenames.csv")
        # # # assert len(os.listdir(input_folder + 'peaks/')) - 1 == len(filenames_df['name'])
        # # for index, row in filemanes_df.iterrows():
        # #     #print(row["name"])
        # #     params.met_reader = ChiPSeqReader(input_folder + 'methylation/'+ row["filename"], name=row['name'])
        # #     params.met_reader.read_file(renamer={"0":"chr","1":"start","2":"end","4":"sigVal"})
        # #     metPG.append(SmallChipSeqPredictorGenerator(params.met_reader,params.window_size,N_closest=4))
        # # #Read cage data
        # # cagePG = []
        # # filemanes_df = pd.read_csv(input_folder + "cage/filenames.csv")
        # # # assert len(os.listdir(input_folder + 'cage/')) - 1 == len(filemanes_df['name'])
        # # for index, row in filemanes_df.iterrows():
        # #     #print(row["name"])
        # #     params.cage_reader = ChiPSeqReader(input_folder+"cage/GSM849365_hg19_wgEncodeRikenCageK562CellPapClusters.bed.gz", name=row['name'])# + "cage/" + row["filename"], name=row['name'])
        # #     params.cage_reader.read_file(renamer={"0":"chr","1":"start","2":"end","4":"sigVal"})
        # #     cagePG.append(SmallChipSeqPredictorGenerator(params.cage_reader,params.window_size,N_closest=4))
        #Read RNA-Seq data
        params.RNAseqReader = RNAseqReader(fname=input_folder + "RNA-seq/rna-seqPolyA.tsvpre.txt",
                                           name="RNA")
        params.RNAseqReader.read_file(rename={ "Gene name": "gene",
                              "Gene start (bp)": "start",
                              "Gene end (bp)": "end",
                              "Chromosome/scaffold name": "chr",
                              "FPKM": "sigVal"},
                      sep="\t")
        RNAseqPG = SmallChipSeqPredictorGenerator(params.RNAseqReader,
                                                  window_size=params.window_size,
                                                  N_closest=3)
        #Read TSS data
        params.TssReader=TssReader(fname=input_folder + "TSS/NCBI_refSeq_hg19.bed", name="TSS")
        params.TssReader.read_file()
        TSSPG=Distance_to_TSS_PG(params.TssReader)


        # #Read E1 data
        # params.eig_reader = E1Reader()
        # e1_files = [input_folder + "E1/"+ "chr" + str(i) + ".GM12878.25K.txt" for i in range(1, 23)]
        # e1_files.append(input_folder +"E1/" + "chrX.GM12878.25K.txt")
        # params.eig_reader.read_files(e1_files, binSizeFromName=fileName2binsize) #infer size of E1 bins from file name using this function
        #
        # e1pg = SmallE1PredictorGenerator(params.eig_reader,params.window_size)
        #
        # params.pgs = [OrientCtcfpg]
        params.pgs = [OrientCtcfpg, NotOrientCTCFpg, OrientBlocksCTCFpg ,RNAseqPG, ConvergentPairPG,TSSPG]+chipPG#+cagePG+metPG+chipPG

        # Generate train
        train_chrs=[]
        [train_chrs.append("chr"+chr) for chr in chr_nums]
        if write_all_chrms_in_file:
            train_file_name="training.RandOn"+ str(params)
            params.out_file=output_folder+"_".join(train_chrs)+train_file_name
        for trainChrName in train_chrs:
            print(trainChrName)
            # training_file_name = "training.RandOn" + trainChrName + str(params) + ".txt"

            params.sample_size = lengths_dict[trainChrName]
            params.interval = Interval(trainChrName,
                                  params.contacts_reader.get_min_contact_position(trainChrName),
                                  params.contacts_reader.get_max_contact_position(trainChrName))
            
            # params.out_file = output_folder + training_file_name
            if not write_all_chrms_in_file:
                train_file_name = "training.RandOn" + str(params) + ".txt"
                params.out_file = output_folder + params.interval.toFileName() + train_file_name
            generate_data(params,saveFileDescription=True)
            if not write_all_chrms_in_file:
                del(params.out_file)
            del (params.sample_size)


        # # Generate test
        # validate_chrs=[]
        # [validate_chrs.append("chr"+chr) for chr in chr_nums]#,"chr16", "chr17"]#, "chr18"]#, "chr18", "chr19", "chr20"]#,"chr14", "chr15"]
        # if write_all_chrms_in_file:
        #     validation_file_name = "validatingOrient." + str(params) + ".txt"
        #     params.out_file = output_folder + "_".join(validate_chrs) + validation_file_name
        # for validateChrName in validate_chrs:
        #     print("chromosome", validateChrName)
        #     # print(params.contacts_reader.data[validateChrName])
        #     # interval=Interval("chr2", 118000000, 129000000)
        #     params.sample_size = len(params.contacts_reader.data[validateChrName])
        #
        #     params.interval = Interval(validateChrName,
        #                                params.contacts_reader.get_min_contact_position(validateChrName),
        #                                params.contacts_reader.get_max_contact_position(validateChrName))
        #     # params.interval = interval
        #     logging.getLogger(__name__).info("Generating validation dataset for interval "+str(params.interval))
        #     if not write_all_chrms_in_file:
        #         validation_file_name = "validatingOrient." + str(params) + ".txt"
        #         params.out_file = output_folder + params.interval.toFileName() + validation_file_name
        #     generate_data(params)
        #     if not write_all_chrms_in_file:
        #         del(params.out_file)
        #     del (params.sample_size)

        # for interval in [Interval("chr2", 118000000, 129000000)]:
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