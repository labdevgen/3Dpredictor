import logging, os
import pandas as pd
import numpy as np
from shared import Interval
import datetime

MAX_CHR_DIST = 3000000000

class ContactsReader():
    def __init__(self):
        self.data = {}
        self.binsize = -1

    def read_file(self,chr,fname,coeff_fname):
        logging.getLogger(__name__).info("Reading file "+fname)
        if chr in self.data:
            logging.getLogger(__name__).warning("Chromosome "+chr+" will be rewritten")
        # print(coeff_data.keys())
        if fname.split(".")[-2]=="oe":
            coeff = 1
        elif fname.split(".")[-2]=="contacts":
            coeff_data = pd.read_csv(coeff_fname, delimiter="\t")
            coeff=coeff_data["coeff"]
        # print(coeff)
        data = pd.read_csv(fname, delimiter="\t", names=["contact_st", "contact_en", "contact_count"])
        data.dropna(inplace=True)
        # print(data["contact_count"])
        logging.info("get normalized contacts")
        print(datetime.datetime.now())
        data["contact_count"] = data["contact_count"].apply(lambda x: x/coeff)
        # print("!!!!")
        print(datetime.datetime.now())
        # print(data["contact_count"])
        data["chr"] = [chr] * len(data)
        data["dist"] = data["contact_en"] - data["contact_st"]
        assert np.all(data["dist"]) >= 0
        assert 0 <= np.all(data["contact_count"]) <= 1
        binsize = min(data["dist"][data["dist"] > 0])
        if self.binsize != -1 and binsize != self.binsize:
            logging.error("Binsize in file "+str(fname)
                          +"("+binsize+") does not match binsize "
                          +str(self.binsize))
        elif self.binsize == -1:
            self.binsize = binsize
            logging.getLogger(__name__).info("Bin Size set to "+str(binsize))

        #data.sort_values(by=["st","en"],inplace=True)
        self.data[chr] = data

    def read_files(self,fnames, coeff_fname):
        for f in fnames:
            self.read_file(os.path.basename(f).split(".")[0],f, coeff_fname=coeff_fname)

    def get_contacts(self,interval,mindist=0,maxdist=MAX_CHR_DIST):
        return self.data[interval.chr].query(
              "@interval.start < contact_st < @interval.end & "
            + "@interval.start < contact_en < @interval.end & "
            + "dist <=@maxdist & "
            + "dist >=@mindist")

    def get_min_contact_position(self,chr):
        return min(self.data[chr]["contact_st"].values)

    def get_max_contact_position(self,chr):
        return max(self.data[chr]["contact_en"].values)

    def get_chrms(self):
        return list(self.data.keys())

    def get_all_chr_contacts(self,chr):
        return self.data[chr]

    def delete_region(self,interval):
        data = self.data[interval.chr]
        #Drop contacts withing interval
        bad_ids = data.query("@interval.start < contact_st < @interval.end | "
            + "@interval.start < contact_en < @interval.end").index #either start or end in region to be removed
        #logging.getLogger(__name__).info (bad_ids)
        data.drop(bad_ids,inplace=True)
        #logging.getLogger(__name__).info(data.head())

        self.data[interval.chr] = data
        #change coordinates
        data["st_red"] = data.contact_st.apply(lambda x: x >= interval.start)
        data["en_red"] = data.contact_en.apply(lambda x: x >= interval.start)

        new_starts = data.contact_st.apply(lambda x: (x - interval.len) if (x >= interval.start) else x).values
        new_ends = data.contact_en.apply(lambda x: (x - interval.len) if (x >= interval.start) else x).values
        new_dist = new_ends - new_starts
        #logging.getLogger(__name__).debug(data.iloc[new_dist < 0,:].head())
        #logging.getLogger(__name__).debug(new_starts[new_dist < 0])
        #logging.getLogger(__name__).debug(new_ends[new_dist < 0])

        assert np.all(new_dist >= 0)

        self.data[interval.chr].loc[:,"contact_st"] = new_starts
        self.data[interval.chr].loc[:,"contact_en"] = new_ends
        self.data[interval.chr].loc[:,"dist"] = new_dist
