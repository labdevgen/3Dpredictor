# Code by Minja Fishman
# Nov 2019, ICG SB RAS
#

# TODO: consider scipy.sparse for future usage

import sys
import os
import straw
import pandas as pd
import numpy as np
import logging
import datetime

# Add main directory to import path
import_path = os.path.dirname(os.path.dirname(os.getcwd()))
logging.getLogger(__name__).info("Appending import path: "+import_path)
sys.path.append(import_path)

from shared import FileReader


class hicReader(FileReader):
    def __init__(self, fname, genome, resolution, maxdist=1000000, normalization = "KR", name = None):
        self.genome = genome
        self.genomeName = genome.name
        self.fname = fname
        self.maxdist = maxdist
        self.resolution = resolution
        self.normalization = normalization
        self.data = {}
        self.norms = {}
        if name == None:
            self.name = os.path.basename(fname)
        self.full_name = str(self.toXMLDict(exludedMembers=("data","norms")))
        super(hicReader, self).__init__(fname)

    def dump(self):
        descriptiveXML = {"XML": self.toXMLDict(exludedMembers=("data", "dataPointer")),
                          "header": self.name}
        temp = self.genome
        del self.genome
        super(hicReader, self).dump(descriptiveXML=descriptiveXML)
        self.genome = temp

    def load(self, genome):
        result = super(hicReader, self).load()
        result.genome = genome
        return result

    def read_data(self, debug_mode = False, noDump = False):
        # if noDump, will not try load data from dump file
        #  debug_mode = False will skip some assertions, useful only for intentionally incorrect input data

        # first try to load data from dump
        if os.path.isfile(self.get_dump_path()) and not noDump:
            return self.load(self.genome)

        # if we found no dump, lets read data and dump file
        for chr in self.genome.chrmSizes.keys():
            logging.getLogger(__name__).info("Processing chrm "+chr)
            load_start_time = datetime.datetime.now()
            try:
                result = straw.straw(self.normalization, self.fname,
                                    chr, chr, "BP", self.resolution)
            except TypeError:
                if "chr" in chr:
                    new_chr = chr.replace("chr","",1)
                    logging.getLogger(__name__).warning("Failed to find chr "+chr+"; trying to find "+new_chr)
                    result = straw.straw(self.normalization, self.fname,
                                    new_chr, new_chr, "BP", self.resolution)
                else:
                    raise TypeError
            logging.getLogger(__name__).debug("Load time: " + str(datetime.datetime.now() - load_start_time))
            now = datetime.datetime.now()

            result = np.array(result).T
            logging.getLogger(__name__).debug("Transpose time: " + str(datetime.datetime.now() - now))
            now = datetime.datetime.now()

            result = pd.DataFrame(result, columns = ["st", "en", "count"], copy=False)
            logging.getLogger(__name__).debug("DF conversion time: " + str(datetime.datetime.now() - now))
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
                logging.getLogger(__name__).warning("Sums of contacs for loci are very different. Examples: ")
                logging.getLogger(__name__).warning(str(s))
                logging.getLogger(__name__).warning("Using average for 'magic normalization coefficient' might be not correct")
            logging.getLogger(__name__).debug("Magic coefficient calc time: " + str(datetime.datetime.now() - now))

            result.query("(en-st)<@self.maxdist",inplace=True)
            assert len(result) > 0

            if not debug_mode:
                assert max(result.en.values) <= self.genome.chrmSizes[chr] + self.resolution
            result["count"] = result["count"] / np.average(s)
            result = result.set_index(["st","en"])
            self.data[chr] = result
            self.norms[chr] = np.average(s)
            logging.getLogger(__name__).info("Total hic load time: "+str(datetime.datetime.now()-load_start_time))

        assert len(self.data.keys()) > 0

        self.dump()
        return self

    def get_contacts(self, interval):
        raise NotImplementedError

    def get_contact(self, interval):
        if interval.start % self.resolution != 0 or interval.end % self.resolution != 0:
            logging.getLogger(__name__).error("Start or end of the contact does not match resolution")
            raise Exception()
        if interval.len >= self.maxdist:
            logging.getLogger(__name__).error("Provided interval: " + str(interval) + "\n" + \
                          " length >= than maximal distance between contacts specified in paramteres ("+\
                          str(self.maxdist)+")")
            raise Exception()
        chr_contacts = self.data[interval.chr]
        try:
            contacts = chr_contacts.loc(axis=0)[(interval.start,interval.end)]
        except KeyError:
            return None # TODO maybe should change this to 0
        if len(contacts) == 1:
            return contacts.iloc[0]
        elif len(contacts) > 1:
            logging.getLogger(__name__).error("More than 1 contact for region "+str(interval))
            raise Exception()
        else:
            raise NotImplementedError

    def get_chr_contact(self, chr):
        return self.data[chr]