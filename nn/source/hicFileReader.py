# Code by Minja Fishman
# Nov 2019, ICG SB RAS

import sys
import os
import straw
import pandas as pd
import numpy as np
import logging

# Add main directory to import path
import_path = os.path.dirname(os.path.dirname(os.getcwd()))
logging.getLogger(__name__).info("Appending import path: "+import_path)
sys.path.append(import_path)

from shared import FileReader


class hicReader(FileReader):
    def __init__(self, fname, genome, resolution, maxdist=1000000, normalization = "KR", name = None):
        self.genome = genome
        if name == None:
            self.name = fname
        self.fname = fname
        self.maxdist = maxdist
        self.resolution = resolution
        self.normalization = normalization
        self.data = {}
        self.norms = {}
        super(hicReader, self).__init__(fname)

    def read_data(self):
        for chr in self.genome.chrmSizes.keys():
            logging.getLogger(__name__).info("Processing chrm "+chr)
            import datetime
            now = datetime.datetime.now()
            try:
                result = straw.straw(self.normalization, self.fname,
                                    chr, chr, "BP", self.resolution)
            except TypeError:
                if "chr" in chr:
                    new_chr = chr.replace("chr","",1)
                    logging.warning("Failed to find chr "+chr+"; trying to find "+new_chr)
                    result = straw.straw(self.normalization, self.fname,
                                    new_chr, new_chr, "BP", self.resolution)
                else:
                    raise TypeError
            logging.debug("Load time: " + str(datetime.datetime.now() - now))
            now = datetime.datetime.now()

            result = np.array(result).T
            logging.debug("Transpose time: " + str(datetime.datetime.now() - now))
            now = datetime.datetime.now()

            result = pd.DataFrame(result, columns = ["st", "en", "count"], copy=False)
            logging.debug("DF conversion time: " + str(datetime.datetime.now() - now))
            now = datetime.datetime.now()

            #result = result.loc[pd.notnull(result["count"])]
            #logging.debug("N/A filtering time: " + str(datetime.datetime.now() - now))
            now = datetime.datetime.now()

            # Let's normalize data to have sum over each contact within chrm ~1.0
            subsample_size = 1000
            subsample = np.random.choice(result.st.values,size=subsample_size,replace=False)
            subsample = np.unique(subsample)
            assert len(subsample) >= subsample_size / 10
            s = []
            for i in subsample:
                local_count = result.query("st==@i | en==@i")["count"].sum()
                if local_count == 0:
                    # these are probably NAN samples
                    continue
                    #logging.error("zero count for region ", i)
                    #logging.debug(str(result.query("st==@i | en==@i")))
                else:
                    s.append(local_count)
            assert len(s) >= len(subsample) / 2
            if np.std(s) / np.average(s) >= 0.2:
                logging.warning("Sums of contacs for loci are very different. Examples: ")
                logging.warning(str(s))
                logging.warning("Using average for 'magic normalization coefficient' might be not correct")
            logging.debug("Magic coefficient calc time: " + str(datetime.datetime.now() - now))
            now = datetime.datetime.now()

            result.query("(en-st)<@self.maxdist",inplace=True)
            assert len(result) > 0
            #TODO uncomment this after debug tests! It will check consistency with genome
            #assert max(result.en.values) <= self.genome.chrmSizes[chr] + self.resolution
            result["count"] = result["count"] / np.average(s)
            self.data[chr] = result
            self.norms[chr] = np.average(s)

        assert len(self.data.keys()) > 0

    def get_contacts(self, interval):
        raise NotImplementedError

    def get_contact(self, interval):
        if interval.start % self.resolution != 0 or interval.end % self.resolution != 0:
            logging.error("Start or end of the contact does not match resolution")
            raise Exception()
        contact = self.data[interval.chr].query("st==@interval.start & en==@interval.end")
        if len(contact) == 1:
            return contact["count"]
        if len(contact) == 0:
            return None
        if len(contact) > 1:
            raise Exception()

    def get_chr_contact(self, chr):
        return self.data[chr]