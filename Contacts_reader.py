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

        data = pd.read_csv(fname, delimiter="\t", names=["st", "en", "count"])
        data["chr"] = [chr] * len(data)
        data["dist"] = data["en"] - data["st"]
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

    def get_contacts(self,interval,maxdist=MAX_CHR_DIST):
        logging.debug(self.data[interval.chr].head())
        return self.data[interval.chr].query(
              "@interval.start < st < @interval.end & "
            + "@interval.start < en < @interval.end & "
            + "dist <=@maxdist")

    def get_min_contact_position(self,chr):
        return min(self.data[chr]["st"].values)

    def get_max_contact_position(self,chr):
        return max(self.data[chr]["en"].values)

    def get_chrms(self):
        return self.data.keys()