# Written by Minja, 08 2018
# A general module purpose is generate predictors
# Class DataGenerator does it by applying contact2file to a dataframe
# which contains contacts information in a format contact_st -- contact_en -- contact_count
# To generate predictors for each contact, we use instances of PredictorGenerators
# Each instance should be able to generate value(s) of specific predicor, e.g. ChiPSeq signal,
# RNAseq of 1st Eigenvecor
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
from shared import write_XML


processed = 0

def _apply_df(args):
    df, DataGeneratorObj = args
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
    generator.contacts2file(contacts_sample, params)
    if saveFileDescription:
        XML_report = generator.toXMLDict()
        write_XML(XML_report,
                  header = params.out_file,
                  fname = params.out_file+".xml")

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
        self.contacts = contacts
        self.contacts.reset_index(drop=True,inplace=True) #This is requered to solve probles with concat later on
        self.params = params

        out_file = open(params.out_file, "w")

        #Check that predictor names are unique
        pg_names = [pg.name for pg in self.predictor_generators]
        #print(pg_names)
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
        pool = multiprocessing.Pool(processes=n_cpus)
        result = pool.map(_apply_df, [(d, self) for d in np.array_split(contacts, n_cpus)])
        #result = list(map(_apply_df, [(d, self) for d in np.array_split(contacts, n_cpus)]))
        pool.close()

        logging.getLogger(__name__).debug("Writing to file")
        for i in result:
            i.apply(lambda x: out_file.write("\t".join(map(str,x))+"\n"), axis="columns")
        for pg in self.predictor_generators:
            pg.print_warnings_occured_during_predGeneration()
        out_file.close()

    def toXMLDict(self):
        if len(self.contacts) == 0:
            raise Exception("Trying to get stats on empty data")

        res = OrderedDict()
        res["date"] = str(datetime.datetime.now())
        res["class"] = self.__class__.__name__
        res["output_file_name"] = self.params.out_file
        res["N_contacts"] = len(self.contacts)
        res["Global paramteres"] = self.params.toXMLDict()
        for pg in self.predictor_generators:
            res["Predictor generator " + pg.name] = pg.toXMLDict(self.contacts.iloc[0,:])
        return res
