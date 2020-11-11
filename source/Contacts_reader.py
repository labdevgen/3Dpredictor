import logging, os
import pandas as pd
import numpy as np
from shared import Interval
import datetime
#import swifter
from ChiPSeqReader import ChiPSeqReader

MAX_CHR_DIST = 3000000000

class ContactsReader(): #Class process files with contacts
    def __init__(self):
        self.data = {}
        self.binsize = -1

    def read_file(self,chr,fname,coeff_fname,fill_empty_contacts, max_cpus, maxdist):
        logging.getLogger(__name__).info("Reading file "+fname)
        if chr in self.data:
            logging.getLogger(__name__).warning("Chromosome "+chr+" will be rewritten")
        if fname.split(".")[-2]=="oe":
            coeff = 1
        elif fname.split(".")[-2]=="contacts":
            coeff_data = pd.read_csv(coeff_fname, delimiter="\t")
            coeff=coeff_data["coeff"]
        data = pd.read_csv(fname, delimiter="\t", names=["contact_st", "contact_en", "contact_count"])
        data.dropna(inplace=True)
        logging.info("get normalized contacts")
        logging.info(datetime.datetime.now())
        data["contact_count"]=data["contact_count"]/int(coeff)
        if fname.split(".")[-2]=="contacts":
            assert 0<=np.all(data["contact_count"])<=1
        logging.info(datetime.datetime.now())
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
            logging.getLogger(__name__).info("Bin Size set to "+str(binsize))

            # get dataframe with all contacts
        if not fill_empty_contacts:
            self.data[chr] = data
        elif fill_empty_contacts:
            min_bin = min(np.min(data["contact_st"].values), np.min(data["contact_en"].values))
            max_bin = max(np.max(data["contact_st"].values), np.min(data["contact_en"].values))
            binsize = self.binsize
            contacts = []
            for i in range(min_bin, max_bin + 1, binsize):
                for j in range(i + binsize, max_bin + 1, binsize):
                    if binsize*2+1<=abs(j-i)<=maxdist:
                        contacts.append((i, j))
            all_contacts = pd.DataFrame(contacts, columns=['contact_st', 'contact_en'])
            result = pd.merge(all_contacts, data, how='left', on=['contact_st', 'contact_en'])
            result["chr"].fillna(chr, inplace=True)
            result["contact_count"].fillna(0, inplace=True)
            result.loc[pd.isna(result["dist"]), "dist"] = result.loc[pd.isna(result["dist"]), "contact_en"] - \
                                                          result.loc[pd.isna(result["dist"]), "contact_st"]
            self.data[chr] = result



    def read_files(self,fnames, coeff_fname, max_cpus, fill_empty_contacts, maxdist):
        for f in fnames:
            self.read_file(os.path.basename(f).split(".")[0],f, coeff_fname=coeff_fname, max_cpus=max_cpus, fill_empty_contacts=fill_empty_contacts,
                           maxdist=maxdist)

    def get_contacts(self,interval,mindist=0,maxdist=MAX_CHR_DIST):
        return self.data[interval.chr].query(
              "@interval.start <= contact_st < @interval.end & "
            + "@interval.start < contact_en <= @interval.end & "
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
        #len_bins = int(round(float(interval.len)/self.binsize,0)) * self.binsize

        bad_ids = data.query("@interval.start < contact_st < @interval.end | "
            + "@interval.start < contact_en < @interval.end").index #either start or end in region to be removed
        #logging.getLogger(__name__).info (bad_ids)
        data.drop(bad_ids,inplace=True)
        logging.getLogger(__name__).info("dropping contacts with index", bad_ids)
        #logging.getLogger(__name__).info(data.head())

        self.data[interval.chr] = data
        #change coordinates

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


    def duplicate_region(self, interval):
        data = self.data[interval.chr]
        old_length=len(self.data[interval.chr])
        bad_ids = data.query("@interval.start < contact_st < @interval.end | "
                             + "@interval.start < contact_en < @interval.end").index
        dup_data= data.loc[bad_ids,:]
        dup_data["contact_st"]+=interval.len
        dup_data["contact_en"] += interval.len
        dup_data["dist"] += interval.len
        # change coordinates
        new_starts = data.contact_st.apply(lambda x: (x + interval.len) if (x >= interval.start) else x).values
        new_ends = data.contact_en.apply(lambda x: (x + interval.len) if (x >= interval.start) else x).values
        new_dist = new_ends - new_starts
        assert np.all(new_dist >= 0)

        self.data[interval.chr].loc[:, "contact_st"] = new_starts
        self.data[interval.chr].loc[:, "contact_en"] = new_ends
        self.data[interval.chr].loc[:, "dist"] = new_dist
        self.data[interval.chr] = pd.concat([self.data[interval.chr], dup_data])
        assert len(self.data[interval.chr]) - len(dup_data)==old_length

    def use_contacts_with_CTCF(self, CTCFfile, maxdist, proportion, keep_only_orient, CTCForientfile):
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
    def generate_contacts_for_region(self, interval, binsize, maxdist):
        self.binsize = binsize
        interval_start_bin = int(interval.start) // int(binsize) * int(binsize)
        interval_end_bin = int(interval.end) // int(binsize) * int(binsize)
        contact_starts = []
        contact_ends = []
        dists = []
        for contact_st in range(interval_start_bin, interval_end_bin + binsize, int(binsize)):
            logging.info(str(datetime.datetime.now()) + " " + str(contact_st))
            for contact_en in range(contact_st, interval_end_bin + binsize, int(binsize)):
                if (contact_en - contact_st) <= maxdist and (contact_en - contact_st) >= maxdist:
                    contact_starts.append(contact_st)
                    contact_ends.append(contact_en)
                    dists.append(contact_en - contact_st)
        assert np.all(dists) >= 0
        dict={'chr':[interval.chr]*len(dists), 'contact_st':contact_starts, 'contact_en':contact_ends, 'contact_count':[1]*len(dists)}
        self.data[interval.chr] = pd.DataFrame(dict)
