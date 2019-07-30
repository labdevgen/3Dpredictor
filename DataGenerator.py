# Written by Minja, 08 2018
# A general module purpose is generate predictors
# Class DataGenerator does it by applying contact2file to a dataframe
# which contains contacts information in a format contact_st -- contact_en -- contact_count
# To generate predictors for each contact, we use instances of PredictorGenerators
# Each instance should be able to generate value(s) of specific predictor, e.g. ChiPSeq signal,
# RNAseq or 1st Eigenvecor
# Function contact2file aggregates these predictors and writes data to file
#
# Generate data is high-level function that generates contacts sample,
# Calls DataGenerator instance and, finally, writes xml description about model

import logging
import numpy as np
import pandas as pd
import datetime
import multiprocessing
from collections import OrderedDict
from shared import write_XML,makedirs
import os

from multiprocessing.reduction import ForkingPickler
import struct


processed = 0

def get_split_array_indexes(contacts,n_splits):
    # this is copy-pasted from np.split_array
    Nsections = n_splits
    Neach_section, extras = divmod(len(contacts), Nsections)
    section_sizes = ([0] +
                     extras * [Neach_section + 1] +
                     (Nsections - extras) * [Neach_section])
    div_points = np.core.numeric.array(section_sizes).cumsum()
    start_points = div_points[0:-1]
    end_points = div_points[1:]
    assert div_points[0] == 0
    assert div_points[-1] == len(contacts)
    return start_points,end_points

def initializer(_main_data, _DataGeneratorObj):
    global main_data
    global DataGeneratorObj
    main_data = _main_data # share _main_data so every child process has access to _main_data under "main_data" name
    DataGeneratorObj = _DataGeneratorObj # same for _DataGeneratorObj

def _apply_df(args):
    #df, DataGeneratorObj = args

    # main data is a big df with all contacts
    # args is a tuple with 2 numbers (from,two), which indicates which part of main data to process in current child process
    df = main_data.iloc[args[0]:args[1]]
    df.reset_index(inplace=True,drop=True) # To allow concat
    # Get not-vectorizable predicors
    basic_info = [df[["chr", "contact_st", "contact_en"]], df["contact_en"] - df["contact_st"],
                    df["contact_count"]]
    not_vect_result = df.apply(contact2file,DataGeneratorObj=DataGeneratorObj,
                    axis="columns")

    # Get vecrorizable predictors
    vect_results = []
    logging.getLogger(__name__).info("Running vectorized predictors")
    for pg in DataGeneratorObj.vect_predictor_generators:
        current_result = pg.get_predictors(df)
        assert len(current_result) == len(df) == len(not_vect_result)
        vect_results.append(current_result)
    return pd.concat(basic_info + [not_vect_result] + vect_results, axis = 1)


def generate_data(params, saveFileDescription = True):
    contacts = params.contacts_reader.get_contacts(params.interval,mindist=params.mindist,maxdist=params.maxdist)
    sample_size = min(params.sample_size,len(contacts))
    logging.getLogger(__name__).info("Using sample size "+str(sample_size))
    contacts_sample = contacts.sample(n=sample_size)
    assert len(contacts_sample) == sample_size
    generator = DataGenerator()
    result = generator.contacts2file(contacts_sample, params)
    if saveFileDescription:
        XML_report = generator.toXMLDict()
        write_XML(XML_report,
                  header = params.out_file,
                  fname = params.out_file+".xml")
    return result

def contact2file(contact,DataGeneratorObj,report = 5000):
        global processed
        processed += 1
        if (processed > report) and (processed % report == 0):
            print(str(datetime.datetime.now())+"Processed: "+str(processed))

        line=[]
        for pg in DataGeneratorObj.not_vect_predictor_generators:
            line += pg.get_predictors(contact)
        if len(line) != DataGeneratorObj.N_notVect_fields:
            logging.error(str(len(line))+" "+str(DataGeneratorObj.N_notVect_fields))
            logging.error(line)
            raise Exception("Length of predictors does not match header")
        return "\t".join(map(str, line))

class DataGenerator():
    def __init__(self,**kwargs):
        for (k,v) in kwargs:
            self.__setattr__(self,k,v)

    def contacts2file(self,contacts,params):
        #contacts - dataframe with contact counts
        #predictor_generators - list, each item is an instance of PredictorGenerator
        if len(contacts) == 0:
            logging.error("Empty contacts dataset")
            raise Exception("Empty contacs dataset")
        logging.getLogger(__name__).info("Writing data to file " + params.out_file)

        #Save some variables if we would like to have stats later on
        self.not_vect_predictor_generators = [p for p in params.pgs if not p.vectorizable]
        self.vect_predictor_generators = [p for p in params.pgs if p.vectorizable]

        self.predictor_generators = self.not_vect_predictor_generators + self.vect_predictor_generators

        # example of contact df row
        # needed for toXML_dict funct
        self.contact_example = contacts.iloc[0,:]
        self.N_contacs = len(contacts)

        contacts.reset_index(drop=True,inplace=True) #This is requered to solve problems with concat later on

        self.params = params
        if os.path.exists(params.out_file):
            out_file = open(params.out_file, "a")
            write_header=False
        else:
            # create directory if it does not exist
            if not os.path.exists(os.path.dirname(params.out_file)):
                makedirs(os.path.dirname(params.out_file))
            out_file = open(params.out_file, "w")
            write_header=True

        #Check that predictor names are unique
        pg_names = [pg.name for pg in self.predictor_generators]
        assert len(pg_names) == len(set(pg_names))

        #Get header row and calculate number of fields
        header = []
        for pg in self.not_vect_predictor_generators:
            header += pg.get_header(contacts.iloc[0,:])
        self.N_notVect_fields = len(header)
        header = ["chr", "contact_st", "contact_en", "contact_dist", "contact_count"] + header
        for pg in self.vect_predictor_generators:
            header += pg.get_header(contacts.iloc[0,:])
        assert len(header) == len(set(header))
        if write_header:
            out_file.write("\t".join(header) + "\n")

        logging.getLogger(__name__).debug("Going to generate predictors for "+ \
                                          str(len(contacts))+" contacts")
        # Calculate CPU's number
        n_cpus = multiprocessing.cpu_count()
        n_cpus = min(len(contacts) // 100 + 1, n_cpus) # don't use many CPUs
                                                # if we have <100 contacts

        if hasattr(params,"max_cpus"): #allow user to limit cpus used
            n_cpus = min(n_cpus,params.max_cpus)
        logging.getLogger(__name__).debug("Number of CPUs set to " + str(n_cpus ))
        logging.getLogger(__name__).debug("Generating contact predictors")

        # Now get predictors
        pool = multiprocessing.Pool(processes=n_cpus,initializer=initializer,initargs=(contacts,self))
        start_points,end_points = get_split_array_indexes(contacts,n_cpus)
        result = pool.map(_apply_df, [(st, end) for st,end in zip(start_points,end_points)])
        #result = pool.map(_apply_df, [(d, self) for d in np.array_split(contacts, n_cpus+10)])
        #result = list(map(_apply_df, [(d, self) for d in np.array_split(contacts, n_cpus)]))
        pool.close()

        logging.getLogger(__name__).debug("Writing to file")
        for i in result:
                i.apply(lambda x: out_file.write("\t".join(map(str,x))+"\n"), axis="columns")
        for pg in self.predictor_generators:
            pg.print_warnings_occured_during_predGeneration()
        out_file.close()
        logging.getLogger(__name__).info("Done!")


    def toXMLDict(self):
        if self.N_contacs == 0:
            raise Exception("Trying to get stats on empty data")

        res = OrderedDict()
        res["date"] = str(datetime.datetime.now())
        res["class"] = self.__class__.__name__
        res["output_file_name"] = self.params.out_file
        res["N_contacts"] = self.N_contacs
        res["Global paramteres"] = self.params.toXMLDict()

        for pg in self.predictor_generators:
            res["Predictor generator " + pg.name] = pg.toXMLDict(self.contact_example)
        return res
