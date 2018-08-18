import logging
from shared import Interval


def contact2file(contact,DataGeneratorObj):
        line = [contact.chr, contact.contact_st, contact.contact_en, contact["contact_count"]]
        for pg in DataGeneratorObj.predictor_generators:
            line += pg.get_predictors(contact)
        assert len(line) == DataGeneratorObj.N_fields
        DataGeneratorObj.out_file.write("\t".join(map(str, line)) + "\n")

class DataGenerator():
    def __init__(self,**kwargs):
        for (k,v) in kwargs:
            self.__setattr__(self,k,v)

    def contacts2file(self,contacts,predictor_generators,out_file_name):
        if len(contacts) == 0:
            logging.error("Empty contacts dataset")
            raise
        logging.info("Writing data to file " + out_file_name)
        self.out_file = open(out_file_name, "w")
        self.predictor_generators = predictor_generators

        header = ["chr", "contact_st", "contact_en", "contact_count"]
        for pg in predictor_generators:
            header += pg.get_header(contacts.iloc[0,:])
        self.N_fields = len(header)
        self.out_file.write("\t".join(header) + "\n")
        contacts.apply(contact2file, DataGeneratorObj=self, axis="columns")
        self.out_file.close()

class PredictorGenerator(object):
    def get_header(self,contact):
        pass
    def get_predictors(self,contact):
        pass
    def symmetric_window_around_contact(self,contact):
        window_start = max(0, contact.contact_st - ((self.window_size - contact.dist) // 2))
        window_end = window_start + self.window_size
        contacts_relative_start = contact.contact_st - window_start
        assert contacts_relative_start >= 0
        contacts_relative_end = contact.contact_en - window_start
        assert contacts_relative_end >= contacts_relative_start
        return window_start,window_end,contacts_relative_start,contacts_relative_end


class CTCFPredictorGenerator(PredictorGenerator):
    def __init__(self,ctcf_reader,binsize,window_size,**kwargs):
        self.ctcf_reader = ctcf_reader
        self.binsize = binsize
        self.window_size = window_size
        super(CTCFPredictorGenerator, self).__init__(**kwargs)

    def get_header(self,contact):
        self.header = ["CTCFwin_st","CTCFwin_en","CTCFwin_rel_st","CTCFwin_rel_en"] \
                      + ["CTCF_bin" + str(i) for i in range(window_size // binsize)]
        return self.header

    def get_predictors(self,contact):
        window_start, window_end, contacts_relative_start, contacts_relative_end = \
                                        self.symmetric_window_around_contact(contact)
        interval = Interval(contact.chr, window_start, window_end)
        return [window_start,window_end,contacts_relative_start,contacts_relative_end] \
               +self.ctcf_reader.get_binned_interval(interval, binsize=self.binsize)

class SmallCTCFPredictorGenerator(CTCFPredictorGenerator):
    # Less predictors:
    # 1. Total CTCF in 3 regions: L_ancor + winsize/2 ... R_ancor - winsize/2
    #                      -winsize/2 Lancor +winsize/2
    #                      -winsize/2 Rancor +winsize/2
    # 2. Distances and sigVral for K nearest CTCF peaks from L_ancor and R_ancor
    def __init__(self,ctcf_reader,binsize,window_size,N_closest,**kwargs):
        self.N_closest = N_closest
        super(SmallCTCFPredictorGenerator, self).__init__(ctcf_reader,binsize,window_size,**kwargs)

    def get_header(self,contact):
        self.header = ["CTCF_L","CTCF_R","CTCF_W"]
        self.header += ["CTCF_"+side+"Top_"+metric+str(i) for i in range(self.N_closest) for side in ["L","R"] for metric in ["D","V"]]
    def get_predictors(self,contact):
        CTCF_L = self.ctcf_reader.get_interval(Interval(contact.chr,
                                                        contact.st-self.window_size//2,
                                                        contact.st+self.window_size//2)
                                               ).sigVal.sum
        CTCF_R = self.ctcf_reader.get_interval(Interval(contact.chr,
                                                        contact.en-self.window_size//2,
                                                        contact.en+self.window_size//2)
                                               ).sigVal.sum
        CTCF_mid = self.ctcf_reader.get_interval(Interval(contact.chr,
                                                        contact.st+self.window_size//2,
                                                        contact.en-self.window_size//2)
                                               ).sigVal.sum
        Left_top = self.ctcf_reader.get_nearest_peaks(Interval(contact.chr,
                                                               contact.st-self.window_size,
                                                               contact.st-self.window_size),
                                                      self.N_closest)
        Left_top = Left_top["sigVal"].values.tolist() + \
                   ((contact.st-self.window_size)-Left_top["mids"]).values.tolist()

        Right_top = self.ctcf_reader.get_nearest_peaks(Interval(contact.chr,
                                                               contact.en+self.window_size,
                                                               contact.en+self.window_size),
                                                       self.N_closest)
        Right_top = Right_top["sigVal"].values.tolist() + \
                    (Left_top["mids"]-(contact.st+self.window_size)).values.tolist()


class E1PredictorGenerator(PredictorGenerator):
    def __init__(self,eig_reader,window_size):
        self.eig_reader = eig_reader
        self.window_size = window_size
        self.header = ["E1win_st","E1win_en","E1win_rel_st","E1win_rel_en"]
        self.header += ["E1_bin" + str(i) for i in range(window_size // eig_reader.get_binsize())]

    def get_header(self,contact):m
        return self.header


    # Less predictors:
    # 1. E1 in +/- winsize region aroung contact ancors
