#Written by Minja, 2018-09
#A reader for RNA seq data
#It should store RNAseq data as Pandas DataFrame:
#Chr -- Start -- End -- sigVal -- gene
#The difference with ChipSeq data is that
#Genes can be long, so we cannot convert whole gene
#To one point called "mid", as we did with ChipSeq

import pandas as pd
import numpy as np
import logging
from ChiPSeqReader import ChiPSeqReader
from shared import intersect_intervals, intersect_with_interval

# TODO :
# This class is inhereted from ChiPSeqReader
# Some functions may work without rewriting, some are already rewritten
# However some are specific to ChiPseq signals. They won't work correct now, and it makes
# no scence to rewrite them because they are not planned to be used
# in future, it's probably better to split ChiPSeqReader to 2 classes
# one that would be parent for both ChiPSeqReader and RNAseqReader with shared functions
# and one for specific ChiPseq funcs, mainly related to oriented data


class RNAseqReader(ChiPSeqReader):
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
     def read_file(self,rename=None,start_chr_name_with_chr=True,*args,**kwargs):
        #Read raw data
        data = pd.read_csv(self.fname,*args,**kwargs)

        #Rename if needed
        if rename != None:
            data.rename(columns=rename,inplace=True)
        #Drop unused fields
        data = data[["chr","start","end","sigVal","gene"]]

        if start_chr_name_with_chr:
            #check wheather chr names include "chr"
            already_correct = np.any(data["chr"].apply(lambda x: x.startswith("chr")).values)
            if not already_correct:
                data["chr"] = data["chr"].apply(lambda x: "chr"+str(x))

        #Check end > start for all genes
        assert np.all(data["end"].values - data["start"] > 0)

        # check duplicates, set mids, and split by chromosomes and sort
        self.chr_data = self.process_data(data)
        del data


     # def get_nearest_peaks(...): - Should work well when inhereted from ChiPSeq reader

     #Depricated, use get_interval
     def _get_interval(self, interval): #Return all genes that intersect interval
                                       #Also counts partial intersections
        search_query = pd.DataFrame({"start":[interval.start],"end":[interval.end]})
        result = intersect_intervals(self.chr_data,
                                     {interval.chr:search_query},
                                     suppreseChrNumberCheck=True)
        gene_idxs = result[interval.chr]["intersection"]
        return self.chr_data[interval.chr].iloc[gene_idxs,:]

     # This function was overloaded from ChiPseq readr
     # Because original ChiPseq func only return those items
     # which mid-pos overlap input interval
     # For large genes midpos may not overalap interval, yet we want to return these gene
     def get_interval(self, interval, return_ids=False): #Return all genes that intersect interval
                                       #Also counts partial intersections
        return intersect_with_interval(self.chr_data,interval, return_ids=return_ids)


     def get_binned_interval(self):
        logging.getLogger(__name__).error("Function not yet ready")
        raise Exception("Not ready")

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