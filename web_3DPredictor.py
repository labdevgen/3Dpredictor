import logging
import pandas as pd
from ChiPSeqReader import ChiPSeqReader
import pickle
from PredictorGenerators import SmallChipSeqPredictorGenerator, SitesOrientPredictorGenerator, OrientBlocksPredictorGenerator,ConvergentPairPredictorGenerator
from shared import Interval, Parameters
from RNASeqReader import RNAseqReader
from DataGenerator import DataGenerator

import sys
import argparse

def createParser():
    parser = argparse.ArgumentParser(formatter_class=argparse.RawDescriptionHelpFormatter, description="3DPredictor: A tool for 3D chromatin structure prediction.", epilog="GitHub: https://github.com/labdevgen/3Dpredictor")
    parser.add_argument('--version', action='version', version="beta")
    subparsers = parser.add_subparsers (title='Commands', dest='command')
    Default_parser = subparsers.add_parser ('Predictor', help='3D Predictor')
    Default_parser.add_argument ('-N', '--rna-seq', required=True, dest="RNA_seq_file", help='RNA-Seq File')
    Default_parser.add_argument ('-C', '--ctcf', required=True, dest="CTCF_file", help='CTCF File')
    Default_parser.add_argument ('-O', '--orient', required=True, dest="CTCF_orient_file", help='CTCF Orientation File')
    Default_parser.add_argument ('-c', '--chrom', required=True, dest="chr", help='Chromosome')
    Default_parser.add_argument ('-s', '--start', required=True, dest="interval_start", help='Start Position')
    Default_parser.add_argument ('-e', '--end', required=True, dest="interval_end", help='End Position')
    Default_parser.add_argument ('-r', '--resolution', required=True, choices=("5000", "25000"), dest="resolution", help='Resolution')
    Default_parser.add_argument ('-d', '--destination', required=True, dest="out_dir", help='Output Folder')
    Default_parser.add_argument ('-m', '--model', required=True, dest="model_path", help='Model Path')
    return parser

logging.basicConfig(format='%(asctime)s %(message)s', datefmt='%I:%M:%S', level=logging.DEBUG)

if __name__ == '__main__':
    parser = createParser()
    namespace = parser.parse_args(sys.argv[1:])
    
    # здесь должна быть проверка
    
    RNA_seq_file = namespace.RNA_seq_file
    CTCF_file = namespace.CTCF_file
    CTCF_orient_file = namespace.CTCF_orient_file
    chr= namespace.chr
    interval_start = namespace.interval_start
    interval_end = namespace.interval_end
    resolution = namespace.resolution
    model_path = namespace.model_path
    out_dir=namespace.out_dir

    #RNA_seq_file = "/mnt/scratch/ws/psbelokopytova/202001051010polina_data/3DPredictor/input/K562/RNA-seq/rna-seqPolyA.tsvpre.txt"
    #CTCF_file = "/mnt/scratch/ws/psbelokopytova/202001051010polina_data/3DPredictor/input/K562/CTCF/wgEncodeAwgTfbsHaibK562CtcfcPcr1xUniPk.narrowPeak.gz"
    #CTCF_orient_file = "/mnt/scratch/ws/psbelokopytova/202001051010polina_data/3DPredictor/input/K562/CTCF/wgEncodeAwgTfbsHaibK562CtcfcPcr1xUniPk.narrowPeak-orient.bed"
    #chr= "chr14"
    #interval_start = "101000000"
    #interval_end = "101500000"
    #resolution = "5000"
    #model_path = "/mnt/scratch/ws/psbelokopytova/202001051010polina_data/3DPredictor/out/models/models up to 15.11.19/5662021667"
    #out_dir="/mnt/scratch/ws/psbelokopytova/202001051010polina_data/3DPredictor/out/"

    params = Parameters()
    params.window_size = int(resolution) #region around contact to be binned for predictors
    params.mindist = int(resolution)*2+1 #minimum distance between contacting regions
    params.maxdist = 1500000 #maximum distance between contacting regions
    params.max_cpus = 1
    # params.keep_only_orient=False #set True if you want use only CTCF with orient
    params.multiprocessing=False
    # params.write_to_file = False
    # Read CTCF data
    # CTCF_file format: ENCODE narrow peak
    # CTCF_orient_file format: chr--start--end--name--score--strand
    logging.info('create CTCF_PG')
    # set path to the CTCF chip-seq file:
    params.ctcf_reader = ChiPSeqReader(CTCF_file, name="CTCF")
    params.ctcf_reader.read_file()
    # set path to the CTCF_orient file:
    params.ctcf_reader.set_sites_orientation(CTCF_orient_file)
    # set corresponding predictor generators and its options:
    OrientCtcfpg = SitesOrientPredictorGenerator(params.ctcf_reader,N_closest=4)
    NotOrientCTCFpg = SmallChipSeqPredictorGenerator(params.ctcf_reader,params.window_size,N_closest=4)
    ctcf_reader_orientOnly = ChiPSeqReader(CTCF_file, name="CTCF")
    ctcf_reader_orientOnly.read_file()
    ctcf_reader_orientOnly.set_sites_orientation(CTCF_orient_file)
    ctcf_reader_orientOnly.keep_only_with_orient_data()
    # set corresponding predictor generators and its options:
    OrientBlocksCTCFpg = OrientBlocksPredictorGenerator(ctcf_reader_orientOnly,params.window_size)
    ConvergentPairPG = ConvergentPairPredictorGenerator(params.ctcf_reader, binsize=params.window_size)

    #Read RNA-Seq data
    #RNA-seq_file format: this file should have fields "gene", "start", "end", "chr","FPKM"
    #you can rename table fields below
    params.RNAseqReader = RNAseqReader(RNA_seq_file,name="RNA")
    #read RNA-seq data and rename table fields
    params.RNAseqReader.read_file(rename={"Gene name": "gene",
                                  "Gene start (bp)": "start",
                                  "Gene end (bp)": "end",
                                  "Chromosome/scaffold name": "chr",
                                  "FPKM": "sigVal"},
                                          sep="\t")
    # set corresponding predictor generators and its options:
    RNAseqPG = SmallChipSeqPredictorGenerator(params.RNAseqReader,window_size=params.window_size,N_closest=3)
    params.pgs = [OrientCtcfpg, NotOrientCTCFpg, OrientBlocksCTCFpg, RNAseqPG, ConvergentPairPG]

    interval_start_bin = int(interval_start)//int(resolution)*int(resolution)
    interval_end_bin = int(interval_end)//int(resolution)*int(resolution)

    logging.info("dump trained model")
    trained_predictor = pickle.load(open(model_path, "rb"))
    outfile = open(out_dir+"predicted_contacts_"+chr+"_"+interval_start+"_"+"interval_end"+"_"+resolution, "w")
    outfile.write("chr"+"\t"+"contact_st"+"\t"+"contact_en"+"\t"+"predicted_contact_count"+"\n")

    for contact_st in range(interval_start_bin, interval_end_bin + params.maxdist+1, int(resolution)):
        for contact_en in range (contact_st, interval_end_bin + params.maxdist+1, int(resolution)):
            if (contact_en - contact_st) <= params.maxdist and (contact_en - contact_st) >= params.mindist:
                contact_df = pd.DataFrame([[chr, contact_st, contact_en, contact_en-contact_st]],columns=['chr', 'contact_st', 'contact_en', 'dist'])
                generator = DataGenerator()
                result = generator.contact2predictors(contact_df, params)
                result_df = pd.DataFrame([[chr, contact_st, contact_en, contact_en-contact_st, 0]+result[1]],
                                         columns=['chr', 'contact_st', 'contact_en', 'contact_dist', "contact_count"]+result[0])

                trained_predictor.validate(result_df, show_plot=False,out_dir=out_dir, df_input=True)
                # print("predicted")
                # print(trained_predictor.predicted)
                outfile.write(str(chr)+"\t"+str(contact_st)+"\t"+str(contact_en)+"\t"+str(trained_predictor.predicted[0])+"\n")
                #     break
        # break
    outfile.close()
