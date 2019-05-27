
import pandas as pd
import numpy as np
import logging
from ChiPSeqReader import ChiPSeqReader
from shared import intersect_intervals, intersect_with_interval
class TssReader(ChiPSeqReader):
     #Init is inhereted from ChiPSeqReader
     #It will set self.data to None
     #And take care about self.fname and self.name

     #read file
     #rename, if not None, passed to dataframe.rename(columns=rename)
     #As a results, dataframe should get column names
     #chr start end  sigVal gene
     #If rename==None dataframe.rename won't be called
     #if start_chr_name_with_chr = True, check wheather chromosome names starts with chr
     #if not, add chr at the begningn
     #args and kwargs are passed to pandas readcsv
     def read_file(self,renamer = {"0":"chr","1":"start","2":"end","5":"strand", "6":"TSS_start", "7":"TSS_end"},
                   start_chr_name_with_chr=True,*args,**kwargs):
         # set random temporary labels
         if self.fname.endswith(".gz"):  # check gzipped files
             import gzip
             temp_file = gzip.open(self.fname)
         else:
             temp_file = open(self.fname)
         Nfields = len(temp_file.readline().strip().split())
         temp_file.close()

         names = list(map(str, list(range(Nfields))))
         data = pd.read_csv(self.fname, sep="\t", header=None, names=names)
         # subset and rename
         data_fields = list(map(int, renamer.keys()))
         data = data.iloc[:, data_fields]
         data.rename(columns=renamer,inplace=True)
         data = data[["chr","start","end","strand","TSS_start", "TSS_end"]]

         if start_chr_name_with_chr:
             #check wheather chr names include "chr"
             already_correct = np.any(data["chr"].apply(lambda x: x.startswith("chr")).values)
             if not already_correct:
                 data["chr"] = data["chr"].apply(lambda x: "chr"+str(x))
         data.loc[data["strand"]=="+", "TSS"] = data["TSS_start"]
         data.loc[data["strand"]=="-", "TSS"] = data["TSS_end"]
         data['TSS'] = data['TSS'].astype('int64')
         #Check end > start for all genes
         assert np.all(data["end"].values - data["start"] > 0)

        # check duplicates, set mids, and split by chromosomes and sort
         chr_data = self.process_data(data)
         del data
         self.chr_data=chr_data


     def get_interval(self, interval, return_ids=False): #Return all genes that intersect interval
                                       #Also counts partial intersections
        return intersect_with_interval(self.chr_data,interval, return_ids=return_ids)

     def delete_region(self, interval):
         debug = len(self.get_interval(interval))
         data = self.chr_data[interval.chr]
         st, en = self.get_interval(interval, return_ids=True)
         self.chr_data[interval.chr].iloc[en:, data.columns.get_loc("start")] -= interval.len
         self.chr_data[interval.chr].iloc[en:, data.columns.get_loc("end")] -= interval.len
         self.chr_data[interval.chr].iloc[en:, data.columns.get_loc("mids")] -= interval.len
         old_length = len(self.chr_data[interval.chr])
         self.chr_data[interval.chr].drop(data.index[st:en], inplace=True)
         assert len(self.chr_data[interval.chr]) + debug == old_length