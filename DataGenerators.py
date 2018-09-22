import logging
from shared import Interval
import numpy as np
import datetime

#A general module purpose is generate prictors
#Class DataGenerator does it by applying contact2file to a dataframe
#which contains contacts information in a format contact_st -- contact_en -- contact_count
#To generate predictors for each contact, it uses instances of PredictorGenerator
#Each instance should be able to generate value(s) of specific predicor, e.g. ChiPSeq signal of 1st Eigenvecor
#Function contact2file aggregates these predictors and writes data to file

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
        if len(contacts) == 0:
            logging.error("Empty contacts dataset")
            raise
        logging.getLogger(__name__).info("Writing data to file " + out_file_name)
        self.out_file = open(out_file_name, "w")
        self.predictor_generators = predictor_generators

        #Check that predictor names are unique
        pg_names = [pg.name for pg in predictor_generators]
        assert len(pg_names) == len(set(pg_names))

        #Get header row and calculate number of fields
        header = ["chr", "contact_st", "contact_en", "contact_dist", "contact_count"]
        for pg in predictor_generators:
            header += pg.get_header(contacts.iloc[0,:])
        self.N_fields = len(header)
        self.out_file.write("\t".join(header) + "\n")

        logging.getLogger(__name__).debug("Going to generate predictors for "+str(len(contacts))+" contacts")
        contacts.apply(contact2file, DataGeneratorObj=self, axis="columns")
        for pg in predictor_generators:
            pg.print_warnings_occured_during_predGeneration()
        self.out_file.close()

class PredictorGenerator(object):
    def __init__(self,**kwargs):
        for key,val in kwargs:
            self.__setattr__(key,val)

    def get_header(self,contact): #Returns predictor header line, i.e. name(s) of predictors
        pass

    def get_predictors(self,contact): #predictor value(s)
        pass

    def symmetric_window_around_contact(self,contact):
        #Returns interval with a size self.window_size center around region spanning contact
        window_start = max(0, contact.contact_st - ((self.window_size - contact.dist) // 2))
        window_end = window_start + self.window_size
        contacts_relative_start = contact.contact_st - window_start
        assert contacts_relative_start >= 0
        contacts_relative_end = contact.contact_en - window_start
        assert contacts_relative_end >= contacts_relative_start
        return window_start,window_end,contacts_relative_start,contacts_relative_end

    def intevals_around_ancor(self,contact):
        half = self.window_size // 2
        assert contact.contact_en - contact.contact_st > half
        return (Interval(contact.chr,
                         contact.contact_st - half,
                         contact.contact_st + half),
                Interval(contact.chr,
                         contact.contact_st + half,
                         contact.contact_en - half),
                Interval(contact.chr,
                         contact.contact_en - half,
                         contact.contact_en + half ),
                )
    def print_warnings_occured_during_predGeneration(self):
        pass


class ChipSeqPredictorGenerator(PredictorGenerator):
    def __init__(self, chipSeq_reader, binsize, window_size, **kwargs):
        self.name =  chipSeq_reader.proteinName
        self.chipSeq_reader = chipSeq_reader
        self.binsize = binsize
        self.window_size = window_size
        super(ChipSeqPredictorGenerator, self).__init__(**kwargs)

    def get_header(self,contact):
        self.header = [self.name + "win_st", self.name + "win_en",
                       self.name + "win_rel_st", self.name + "win_rel_en"] \
                      + [self.name + "_bin" + str(i) for i in range(self.window_size // self.binsize)]
        return self.header

    def get_predictors(self,contact):
        window_start, window_end, contacts_relative_start, contacts_relative_end = \
                                        self.symmetric_window_around_contact(contact)
        interval = Interval(contact.chr, window_start, window_end)
        return [window_start,window_end,contacts_relative_start,contacts_relative_end] \
               +self.chipSeq_reader.get_binned_interval(interval, binsize=self.binsize)

    def delete_region(self,interval):
        self.chipSeq_reader.delete_region(interval)

class SmallChipSeqPredictorGenerator(ChipSeqPredictorGenerator):
    # Less predictors:
    # 1. Total CTCF in 3 regions: L_ancor + winsize/2 ... R_ancor - winsize/2
    #                      -winsize/2 Lancor +winsize/2
    #                      -winsize/2 Rancor +winsize/2
    # 2. Distances and sigVral for K nearest CTCF peaks from L_ancor and R_ancor
    def __init__(self, chipSeq_reader, window_size, N_closest, **kwargs):
        self.name = chipSeq_reader.proteinName
        self.N_closest = N_closest
        super(SmallChipSeqPredictorGenerator, self).__init__(chipSeq_reader, 0, window_size, **kwargs)

    def get_header(self,contact):
        self.header = [self.name + "_L", self.name + "_W", self.name + "_R"]
        for side in ["L", "R"]:
            for metric in ["Sig", "Dist"]:
                for i in range(self.N_closest):
                    self.header += [self.name + "_" + side + metric + "_" + str(i)]
        return self.header

    def get_predictors(self,contact):
        intL,intM,intR = self.intevals_around_ancor(contact)
        sig_L = self.chipSeq_reader.get_interval(intL).sigVal.sum()
        sig_R = self.chipSeq_reader.get_interval(intR).sigVal.sum()
        sig_mid = self.chipSeq_reader.get_interval(intM).sigVal.sum()
        Left_top = self.chipSeq_reader.get_nearest_peaks(Interval(contact.chr,
                                                                  contact.contact_st - (self.window_size // 2),
                                                                  contact.contact_st - (self.window_size // 2)),
                                                         N=self.N_closest, side="left")

        Left_top = Left_top["sigVal"].values.tolist() + \
                   (contact.contact_st-Left_top["mids"]).values.tolist()

        Right_top = self.chipSeq_reader.get_nearest_peaks(Interval(contact.chr,
                                                                   contact.contact_en + (self.window_size // 2),
                                                                   contact.contact_en + (self.window_size // 2)),
                                                          N=self.N_closest, side="right")
        Right_top = Right_top["sigVal"].values.tolist() + \
                    (Right_top["mids"]-contact.contact_en).values.tolist()

        return [sig_L,sig_mid,sig_R]+Left_top+Right_top


class E1PredictorGenerator(PredictorGenerator):
    def __init__(self,eig_reader,window_size,name="E1"):
        self.eig_reader = eig_reader
        self.window_size = window_size
        self.name = name

    def get_header(self,contact):
        self.header = [self.name + "win_st",self.name + "win_en",self.name + "win_rel_st",self.name + "win_rel_en"]
        self.header += [self.name + "_bin" + str(i) \
                         for i in range(self.window_size // self.eig_reader.get_binsize())]
        return self.header

    def get_predictors(self,contact):
        window_start, window_end, contacts_relative_start, contacts_relative_end = \
                                        self.symmetric_window_around_contact(contact)
        interval = Interval(contact.chr, window_start, window_end)
        return [window_start, window_end, contacts_relative_start, contacts_relative_end] + \
                self.eig_reader.get_E1inInterval(interval)["E1"].tolist()

    def print_warnings_occured_during_predGeneration(self):
        self.eig_reader.print_warnings()

    def delete_region(self,interval):
        self.eig_reader.delete_region(interval)

class SmallE1PredictorGenerator(E1PredictorGenerator):
    # Less predictors:
    # 1. E1 in +/- winsize region around contact ancors
    def get_header(self,*args):
        self.header = [self.name + "_L",self.name + "_R"]
        return self.header

    def get_predictors(self,contact):
        intL,intM,intR = self.intevals_around_ancor(contact)
        return list(map(np.average,[self.eig_reader.get_E1inInterval(intL)["E1"].tolist(),
               self.eig_reader.get_E1inInterval(intR)["E1"].tolist()]))