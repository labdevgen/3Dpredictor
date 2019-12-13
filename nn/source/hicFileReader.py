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

# Add main source directory to import path
root_dir = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
source_dir = os.path.join(root_dir,"source")
sys.path.append(source_dir)

from shared import FileReader


class hicReader(FileReader):
    def __init__(self, fname, genome, binsize, maxdist=1000000, normalization = "KR", name = None,
                 indexedData = False):
        # indexedData - store data in new indexed format
        # allows fast extraction of single contact, but not compatible with most of contact_reader function
        self.genome = genome
        self.genomeName = genome.name
        self.fname = fname
        self.maxdist = maxdist
        self.binsize = binsize
        self.normalization = normalization
        self.indexedData = indexedData
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

    def read_data(self, debug_mode = False, noDump = False, fill_empty_contacts = False):
        # if noDump, will not try load data from dump file
        #  debug_mode = False will skip some assertions, useful only for intentionally incorrect input data

        if fill_empty_contacts:
            raise NotImplemented

        # first try to load data from dump
        if os.path.isfile(self.get_dump_path()) and not noDump:
            return self.load(self.genome)

        # if we found no dump, lets read data and dump file
        for chr in self.genome.chrmSizes.keys():
            logging.getLogger(__name__).info("Processing chrm "+chr)
            load_start_time = datetime.datetime.now()
            try:
                result = straw.straw(self.normalization, self.fname,
                                    chr, chr, "BP", self.binsize)
            except TypeError:
                if "chr" in chr:
                    new_chr = chr.replace("chr","",1)
                    logging.getLogger(__name__).warning("Failed to find chr "+chr+"; trying to find "+new_chr)
                    result = straw.straw(self.normalization, self.fname,
                                    new_chr, new_chr, "BP", self.binsize)
                else:
                    raise TypeError
            logging.getLogger(__name__).debug("Load time: " + str(datetime.datetime.now() - load_start_time))
            now = datetime.datetime.now()

            result = np.array(result).T
            logging.getLogger(__name__).debug("Transpose time: " + str(datetime.datetime.now() - now))
            now = datetime.datetime.now()
            result = pd.DataFrame(result, columns = ["contact_st", "contact_en", "contact_count"], copy=False)
            logging.getLogger(__name__).debug("DF conversion time: " + str(datetime.datetime.now() - now))
            now = datetime.datetime.now()

            #result = result.loc[pd.notnull(result["count"])]
            #logging.debug("N/A filtering time: " + str(datetime.datetime.now() - now))
            now = datetime.datetime.now()

            # Let's normalize data to have sum over each contact within chrm ~1.0
            subsample_size = 1000
            subsample = np.random.choice(result.contact_st.values,size=subsample_size,replace=False)
            subsample = np.unique(subsample)
            assert len(subsample) >= subsample_size / 10
            s = []
            for i in subsample:
                local_count = result.query("contact_st==@i | contact_en==@i")["contact_count"].sum()
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

            result.query("(contact_en-contact_st)<@self.maxdist",inplace=True)
            assert len(result) > 0

            if not debug_mode:
                assert max(result.en.values) <= self.genome.chrmSizes[chr] + self.binsize
            result["count"] = result["contact_count"] / np.average(s)
            result["dist"] = result["contact_en"] - result["contact_st"]
            assert np.all(result["dist"].values>=0)
            if self.indexedData:
                result = result.set_index(["contact_st","contact_en"])
            self.data[chr] = result
            self.norms[chr] = np.average(s)
            logging.getLogger(__name__).info("Total hic load time: "+str(datetime.datetime.now()-load_start_time))

        assert len(self.data.keys()) > 0

        self.dump()
        return self

    def get_contacts(self, interval):
        raise NotImplementedError

    def get_contact(self, interval):
        # note: this won't work for not-indexed data
        assert self.indexedData
        if interval.start % self.binsize != 0 or interval.end % self.binsize != 0:
            logging.getLogger(__name__).error("Start or end of the contact does not match binsize")
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

    def read_file(self):
        logging.getLogger(__name__).error("This function is not implemented")
        raise Exception()

    def get_binsize(self):
        return self.binsize