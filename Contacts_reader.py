import logging, os
import pandas as pd
import numpy as np
from shared import Interval

MAX_CHR_DIST = 3000000000

class ContactsReader():
    def __init__(self):
        self.data = {}
        self.binsize = -1

    def read_file(self,chr,fname):
        logging.info("Reading file "+fname)
        if chr in self.data:
            logging.warning("Chromosome "+chr+" will be rewritten")

        data = pd.read_csv(fname, delimiter="\t", names=["contact_st", "contact_en", "contact_count"])
        data.dropna(inplace=True)
        data["chr"] = [chr] * len(data)
        data["dist"] = data["contact_en"] - data["contact_st"]
        assert np.all(data["dist"]) >= 0
        binsize = min(data["dist"][data["dist"] > 0])
        if self.binsize != -1 and binsize != self.binsize:
            logging.error("Binsize in file "+str(fname)
                          +"("+binsize+") does not match binsize "
                          +str(self.binsize))
        elif self.binsize == -1:
            self.binsize = binsize
            logging.info("Bin Size set to "+str(binsize))

        #data.sort_values(by=["st","en"],inplace=True)
        self.data[chr] = data

    def read_files(self,fnames):
        for f in fnames:
            self.read_file(os.path.basename(f).split(".")[0],f)

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