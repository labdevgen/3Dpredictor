# Written by Minja, 08 2018
# A general module purpose is generate predictors
# Class DataGenerator does it by applying contact2file to a dataframe
# which contains contacts information in a format contact_st -- contact_en -- contact_count
# To generate predictors for each contact, we use instances of PredictorGenerators
# Each instance should be able to generate value(s) of specific predicor, e.g. ChiPSeq signal,
# RNAseq of 1st Eigenvecor
# Function contact2file aggregates these predictors and writes data to file

import logging
import datetime

processed = 0

def contact2file(contact,DataGeneratorObj,report = 5000):
        global processed
        processed += 1
        if (processed > report) and (processed % report == 0):
            print(str(datetime.datetime.now())+"Processed: "+str(processed))

        line = [contact.chr, contact.contact_st, contact.contact_en,
                contact.contact_en - contact.contact_st, contact["contact_count"]]

        for pg in DataGeneratorObj.predictor_generators:
            line += pg.get_predictors(contact)
        if len(line) != DataGeneratorObj.N_fields:
            logging.error(str(len(line))+" "+str(DataGeneratorObj.N_fields))
            logging.error(line)
            raise Exception("Length of predictors does not match header")
        DataGeneratorObj.out_file.write("\t".join(map(str, line)) + "\n")

class DataGenerator():
    def __init__(self,**kwargs):
        for (k,v) in kwargs:
            self.__setattr__(self,k,v)

    def contacts2file(self,contacts,predictor_generators,out_file_name):
        #contacts - dataframe with contact counts
        #predictor_generators - list, each item is an instance of PredictorGenerator
        if len(contacts) == 0:
            logging.error("Empty contacts dataset")
            raise Exception("Empty contacs dataset")
        logging.getLogger(__name__).info("Writing data to file " + out_file_name)
        self.out_file = open(out_file_name, "w")
        self.predictor_generators = predictor_generators

        #Check that predictor names are unique
        pg_names = [pg.name for pg in predictor_generators]
        #print(pg_names)
        assert len(pg_names) == len(set(pg_names))

        #Get header row and calculate number of fields
        header = ["chr", "contact_st", "contact_en", "contact_dist", "contact_count"]
        for pg in predictor_generators:
            header += pg.get_header(contacts.iloc[0,:])
        #print(header)
        assert len(header) == len(set(header))
        self.N_fields = len(header)
        self.out_file.write("\t".join(header) + "\n")

        logging.getLogger(__name__).debug("Going to generate predictors for "+str(len(contacts))+" contacts")
        contacts.apply(contact2file, DataGeneratorObj=self, axis="columns")
        for pg in predictor_generators:
            pg.print_warnings_occured_during_predGeneration()
        self.out_file.close()
