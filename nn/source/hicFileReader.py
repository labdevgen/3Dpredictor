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
from distutils.version import LooseVersion

# Add main source directory to import path
root_dir = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
source_dir = os.path.join(root_dir,"source")
sys.path.append(source_dir)

from shared import FileReader
from ChiPSeqReader import ChiPSeqReader
from Contacts_reader import ContactsReader

MAX_CHR_DIST = 3000000000

class hicReader(FileReader):
    def __init__(self, fname, genome, binsize, maxdist=1000000, normalization = "KR",
                 name = None, indexedData = False, removeNans = True):
        # indexedData - store data in new indexed format
        # allows fast extraction of single contact, but not compatible with most of contact_reader function
        # dropNans - drop contacts with contact frequency = Nan
        # Note that contacts are stored in sparse matrix format. This means that all missing contacts = 0
        # When we drop Nans, we are unable to distinguish between contacts equal to 0
        # and contacts with KR norm values equal to Nan (i.e. contacts within white lines on hic-maps)

        self.genome = genome
        self.genomeName = genome.name
        self.fname = fname
        self.maxdist = maxdist
        self.binsize = binsize
        self.normalization = normalization
        self.indexedData = indexedData
        self.norms = {}
        self.data = {}
        self.droppedNans = removeNans
        if name == None:
            self.name = os.path.basename(fname)
        self.full_name = str(sorted(self.toXMLDict(exludedMembers=("data","norms")).items()))
        super(hicReader, self).__init__(fname)

    def dropNans(self):
        # drop all contacts with Nan values
        # IMPORTANT: see comments in init func explaining why this leads to confusion of 0 and Nan contacts
        for chr in self.data.keys():
            self.data[chr] = self.data[chr][pd.notna(self.data[chr]["contact_count"])]
            assert len(self.data[chr]) > 0
        self.droppedNans = True

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
        # debug_mode = False will skip some assertions, useful only for intentionally incorrect input data
        # returns object with "self.data" containig data

        def get_data_straw006():
            try:
                return straw.straw(self.normalization, self.fname,
                                     chr, chr, "BP", self.binsize)
            except TypeError:
                if "chr" in chr:
                    new_chr = chr.replace("chr", "", 1)
                    logging.getLogger(__name__).warning("Failed to find chr " + chr + "; trying to find " + new_chr)
                    if new_chr in chr:
                        return straw.straw(self.normalization, self.fname,
                                             new_chr, new_chr, "BP", self.binsize)
                    else:
                        return None

        def get_data_straw008():
            strawObj = straw.straw(self.fname)

            # check chr exists in hic file
            if not chr in strawObj.chromDotSizes.data.keys():
                if "chr" in chr:
                    new_chr = chr.replace("chr", "", 1)
                    logging.getLogger(__name__).warning("Failed to find chr " + chr + "; trying to find " + new_chr)
                    if not new_chr in chr:
                        return None
                else:
                    return None
                hic_chr = new_chr
            else:
                hic_chr = chr

            # check chromosome sizes
            if self.genome.chrmSizes[chr] != strawObj.chromDotSizes.getLength(hic_chr):
                logging.getLogger(__name__).error("Genome version mismatch!")
                raise BaseException

            chr1, X1, X2 = strawObj.chromDotSizes.figureOutEndpoints(hic_chr) #get start & end of chromosome
            matrxObj = strawObj.getNormalizedMatrix(chr1, chr, self.normalization,  # chr1 = chr = hic_chr
                                                    "BP", self.binsize)

            assert matrxObj is not None
            return matrxObj.getDataFromGenomeRegion(X1, X2, X1, X2) #return data for the chromosome

        if fill_empty_contacts:
            raise NotImplemented

        # first try to load data from dump
        if os.path.isfile(self.get_dump_path()) and (not noDump):
            return self.load(self.genome)

        # if we found no dump, lets read data and dump file

        # first define straw version
        if LooseVersion(straw.__version__) == "0.0.6":
            get_data_straw = get_data_straw006
        elif LooseVersion(straw.__version__) >= "0.0.8":
            get_data_straw = get_data_straw008
        else:
            logging.error("Unsupported straw version. "+ straw.__version__ +\
                          "Please use straw 0.0.6 or >=0.0.8")
            raise ValueError

        # now get and process data using straw
        for chr in self.genome.chrmSizes.keys():
            logging.getLogger(__name__).info("Processing chrm "+chr)
            load_start_time = datetime.datetime.now()
            result = get_data_straw()
            print(str(result)) #A
            if result is None:
                logging.getLogger(__name__).warning("Failed to find chr " + chr + " in hic file!")
                continue
            logging.getLogger(__name__).warning("Failed to find chr " + chr + " in hic file!")

            logging.getLogger(__name__).debug("Load time: " + str(datetime.datetime.now() - load_start_time))
            now = datetime.datetime.now()

            result = np.array(result).T
            logging.getLogger(__name__).debug("Transpose time: " + str(datetime.datetime.now() - now))
            now = datetime.datetime.now()
            result = pd.DataFrame(result, columns = ["contact_st", "contact_en", "count"], copy=False)
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
                local_count = result.query("contact_st==@i | contact_en==@i")["count"].sum()
                if local_count == 0:
                    # these are probably NAN samples
                    continue
                    #logging.error("zero count for region ", i)
                    #logging.debug(str(result.query("st==@i | en==@i")))
                else:
                    s.append(local_count)
            #assert len(s) >= len(subsample) / 2
            if np.std(s) / np.average(s) >= 0.2:
                logging.getLogger(__name__).warning("Sums of contacs for loci are very different. Examples: ")
                logging.getLogger(__name__).warning(str(s))
                logging.getLogger(__name__).warning("Using average for 'magic normalization coefficient' might be not correct")
            logging.getLogger(__name__).debug("Magic coefficient calc time: " + str(datetime.datetime.now() - now))

            result.query("(contact_en-contact_st)<@self.maxdist",inplace=True)
            assert len(result) > 0

            if not debug_mode:
                assert max(result.contact_en.values) <= self.genome.chrmSizes[chr] + self.binsize
            result["contact_count"] = result["count"] / np.average(s)
            result["chr"] = [str(chr)] * len(result)
            result["dist"] = result["contact_en"] - result["contact_st"]
            assert np.all(result["dist"].values>=0)
            if self.indexedData:
                result = result.set_index(["contact_st","contact_en"])
                assert result.index.is_unique
            self.data[chr] = result

            if fill_empty_contacts:
                logging.getLogger(__name__).info("going to fill empty contacts")
                print("i'm here!!!")
                self.data[chr]["contact_st"] = self.data[chr]["contact_st"].apply(lambda x: int(x))
                self.data[chr]["contact_st"] = self.data[chr]["contact_en"].apply(lambda x: int(x))
                min_bin = min(np.min(self.data[chr]["contact_st"].values), np.min(self.data[chr]["contact_en"].values))
                max_bin = max(np.max(self.data[chr]["contact_st"].values), np.min(self.data[chr]["contact_en"].values))
                binsize = self.binsize
                contacts = []
                for i in range(min_bin, max_bin + 1, binsize):
                    for j in range(i + binsize, max_bin + 1, binsize):
                        if binsize * 2 + 1 <= abs(j - i) <= self.maxdist:
                            contacts.append((i, j))
                all_contacts = pd.DataFrame(contacts, columns=['contact_st', 'contact_en'])
                result_merge = pd.merge(all_contacts, self.data[chr], how='left', on=['contact_st', 'contact_en'])
                result_merge["chr"].fillna(chr, inplace=True)
                result_merge["contact_count"].fillna(0, inplace=True)
                result_merge.loc[pd.isna(result_merge["dist"]), "dist"] = result_merge.loc[pd.isna(result_merge["dist"]), "contact_en"] - \
                                                              result_merge.loc[pd.isna(result_merge["dist"]), "contact_st"]
                self.data[chr] = result_merge
            self.norms[chr] = np.average(s)
            logging.getLogger(__name__).info("Total hic load time: "+str(datetime.datetime.now()-load_start_time))
        assert len(self.data.keys()) > 0
        self.dropNans()
        self.dump()
        return self

    def get_contacts(self, interval, mindist=0, maxdist=MAX_CHR_DIST):
        assert not self.indexedData
        return self.data[interval.chr].query(
            "@interval.start <= contact_st < @interval.end & "
            + "@interval.start < contact_en <= @interval.end & "
            + "dist <=@maxdist & "
            + "dist >=@mindist")

    def get_contact(self, interval):
        # interval - should contain interval, start (left anchor) and end (right anchor) of two loci
        # returns number of contacts between the left and right anchor
        #
        # note: this function won't work for not-indexed data
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
            return chr_contacts.loc()[(interval.start,interval.end),"contact_count"]
            # Uncomment these if you would like to return dataframe with contact(s) information
            # contacts = chr_contacts.loc()[[(interval.start,interval.end)]]
        except KeyError:
            return None # TODO maybe should change this to 0
        # Uncomment these if you would like to return dataframe with contact(s) information
        #if len(contacts) == 1:
        #    return contacts.iloc[0]
        # elif len(contacts) > 1:
        #     logging.getLogger(__name__).error(chr_contacts.head())
        #     logging.getLogger(__name__).error("More than 1 contact for region "+str(interval))
        #     logging.getLogger(__name__).error(str(contacts))
        #     logging.getLogger(__name__).error(str(len(contacts)))
        #     raise Exception()
        # else:
        #     raise NotImplementedError

    def read_file(self):
        logging.getLogger(__name__).error("This function is not implemented")
        raise Exception()

    def get_binsize(self):
        return self.binsize

    def get_min_contact_position(self,chr):
        print(self.data.keys())
        return min(self.data[chr]["contact_st"].values)

    def get_max_contact_position(self,chr):
        return max(self.data[chr]["contact_en"].values)

    def get_chrms(self):
        return list(self.data.keys())

    def get_all_chr_contacts(self,chr):
        return self.data[chr]

    def use_contacts_with_CTCF(self, CTCFfile, maxdist, proportion, keep_only_orient, CTCForientfile):
        #keep_only_orient [True or False] use only CTCF binding sites with known orientation,
        # required CTCForientfile from GimmeMotifs

        mindist = self.binsize * 2 + 1
        ctcf_reader = ChiPSeqReader(CTCFfile)
        ctcf_reader.read_file()
        if keep_only_orient:
            ctcf_reader.set_sites_orientation(CTCForientfile)
            ctcf_reader.keep_only_with_orient_data()
        conts_with_ctcf = []
        for chr in self.data.keys():
            contacts_data=self.data[chr]
            ctcf_data=ctcf_reader.chr_data[chr]
            ctcf_bins=[] #list of bins which contain CTCF
            ctcf_data["mids"].apply(lambda x: ctcf_bins.extend([x//self.binsize*self.binsize, x//self.binsize*self.binsize+self.binsize,
                                                            x//self.binsize*self.binsize-self.binsize]))
            assert len(ctcf_data)*3==len(ctcf_bins)
            ctcf_bins=sorted(list(set(ctcf_bins)))
            contacts_with_ctcf=[]
            for i in range(0, len(ctcf_bins)):
                for j in range(i + 1, len(ctcf_bins)):
                    if self.binsize * 2 + 1 <= abs(ctcf_bins[j] - ctcf_bins[i]) <= maxdist:
                        contacts_with_ctcf.append((ctcf_bins[i], ctcf_bins[j]))
            contacts_with_ctcf_df = pd.DataFrame(contacts_with_ctcf, columns=['contact_st', 'contact_en'])
            merging_dfs = pd.merge(contacts_data, contacts_with_ctcf_df, how='outer', on=['contact_st', 'contact_en'], indicator=True)
            df_with_CTCF=merging_dfs[merging_dfs["_merge"]=="both"].query("dist <=@maxdist & dist >=@mindist")
            df_wo_CTCF=merging_dfs[merging_dfs["_merge"]=="left_only"].query("dist <=@maxdist & dist >=@mindist")
            df_wo_CTCF =df_wo_CTCF.sample(n=len(df_with_CTCF)*proportion)
            result=pd.concat([df_with_CTCF, df_wo_CTCF])
            self.data[chr]=result
            conts_with_ctcf.append(len(df_with_CTCF))
        self.conts_with_ctcf= np.sum(conts_with_ctcf)

    def delete_region(self,interval):
        data = self.data[interval.chr]
        #Drop contacts withing interval
        #len_bins = int(round(float(interval.len)/self.binsize,0)) * self.binsize

        bad_ids = data.query("@interval.start < contact_st < @interval.end | "
            + "@interval.start < contact_en < @interval.end").index #either start or end in region to be removed
        #logging.getLogger(__name__).info (bad_ids)
        data.drop(bad_ids,inplace=True)
        logging.getLogger(__name__).info("dropping contacts with index", bad_ids)
        #logging.getLogger(__name__).info(data.head())

        self.data[interval.chr] = data
        #change coordinates
        # TODO remove following 2 commands
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