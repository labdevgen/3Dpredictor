import logging, os
import pandas as pd
import numpy as np
from shared import Interval
import datetime
import swifter
from ChiPSeqReader import ChiPSeqReader

MAX_CHR_DIST = 3000000000

class ContactsReader():
    def __init__(self):
        self.data = {}
        self.binsize = -1

    def read_file(self,chr,fname,coeff_fname,fill_empty_contacts, max_cpus, maxdist):
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
        # print(data["contact_count"])
        data["contact_count"]=data["contact_count"]/int(coeff)
        # print(data["contact_count"])
        # data["contact_count"] = data["contact_count"].swifter.set_npartitions(max_cpus).apply(lambda x: x/coeff )
        if fname.split(".")[-2]=="contacts":
            assert 0<=np.all(data["contact_count"])<=1
        print(datetime.datetime.now())
        # print("done")
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
            print("min_bin", min_bin)
            print("max_bin", max_bin)
            print("binsize", binsize)
            contacts = []
            for i in range(min_bin, max_bin + 1, binsize):
                for j in range(i + binsize, max_bin + 1, binsize):
                    if binsize*2+1<=abs(j-i)<=maxdist:
                        contacts.append((i, j))
            print("get result")
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
        # print("get contacts !!!!!!!!!!!!!!!")
        # print(len(self.data[interval.chr]))
        # print(len(self.data[interval.chr].query(
        #       "dist <=@maxdist & "
        #     + "dist >=@mindist")))
        # print(interval.start, interval.end)
        # print(len(self.data[interval.chr].query("19025000 <= contact_st < 115105000 & "
        #     + "19025000 < contact_en <= 115105000")))
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
        print("dropping contacts with index", bad_ids)
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


    def duplicate_region(self, interval):
        data = self.data[interval.chr]
        old_length=len(self.data[interval.chr])
        bad_ids = data.query("@interval.start < contact_st < @interval.end | "
                             + "@interval.start < contact_en < @interval.end").index
        print(bad_ids)
        dup_data= data.loc[bad_ids,:]
        print("dup", len(dup_data))
        dup_data["contact_st"]+=interval.len
        dup_data["contact_en"] += interval.len
        dup_data["dist"] += interval.len
        # change coordinates
        new_starts = data.contact_st.apply(lambda x: (x + interval.len) if (x >= interval.start) else x).values
        new_ends = data.contact_en.apply(lambda x: (x + interval.len) if (x >= interval.start) else x).values
        new_dist = new_ends - new_starts
        # logging.getLogger(__name__).debug(data.iloc[new_dist < 0,:].head())
        # logging.getLogger(__name__).debug(new_starts[new_dist < 0])
        # logging.getLogger(__name__).debug(new_ends[new_dist < 0])
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
            # print("ctcf_reader", ctcf_data)
            ctcf_bins=[] #list of bins which contain CTCF
            ctcf_data["mids"].apply(lambda x: ctcf_bins.extend([x//self.binsize*self.binsize, x//self.binsize*self.binsize+self.binsize,
                                                            x//self.binsize*self.binsize-self.binsize]))
            assert len(ctcf_data)*3==len(ctcf_bins)
            ctcf_bins=sorted(list(set(ctcf_bins)))
            # print("ctcf_bins",ctcf_bins)
            contacts_with_ctcf=[]
            for i in range(0, len(ctcf_bins)):
                for j in range(i + 1, len(ctcf_bins)):
                    if self.binsize * 2 + 1 <= abs(ctcf_bins[j] - ctcf_bins[i]) <= maxdist:
                        contacts_with_ctcf.append((ctcf_bins[i], ctcf_bins[j]))
            contacts_with_ctcf_df = pd.DataFrame(contacts_with_ctcf, columns=['contact_st', 'contact_en'])
            merging_dfs = pd.merge(contacts_data, contacts_with_ctcf_df, how='outer', on=['contact_st', 'contact_en'], indicator=True)
            print(len(merging_dfs[merging_dfs["_merge"]=="both"]))
            df_with_CTCF=merging_dfs[merging_dfs["_merge"]=="both"].query("dist <=@maxdist & dist >=@mindist")
            print("len(df_with_CTCF))", len(df_with_CTCF))
            df_wo_CTCF=merging_dfs[merging_dfs["_merge"]=="left_only"].query("dist <=@maxdist & dist >=@mindist")
            df_wo_CTCF =df_wo_CTCF.sample(n=len(df_with_CTCF)*proportion)
            print("len(df_wo_CTCF))", len(df_wo_CTCF))
            result=pd.concat([df_with_CTCF, df_wo_CTCF])
            self.data[chr]=result
            conts_with_ctcf.append(len(df_with_CTCF))
        self.conts_with_ctcf= np.sum(conts_with_ctcf)
            # print(self.data)


