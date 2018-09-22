#Written by Minja, 2018-09
#A reader for RNA seq data
#It should store RNAseq data as Pandas DataFrame:
#Chr -- Start -- End -- Value
#The difference with ChipSeq data is that
#Genes can be long, so we cannot convert whole gene
#To one point called "mid", as we did with ChipSeq

import pandas as pd
import numpy as np
import logging
from ChiPSeqReader import ChiPSeqReader
from shared import intersect_intervals

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


     def get_nearest_peaks(self):
        logging.getLogger(__name__).error("Function not yet ready")
        raise Exception("Not ready")

     def get_interval(self,interval): #Return all genes that intersect interval
                        #Also counts partiall intersections
        search_query = pd.DataFrame({"start":[interval.start],"end":[interval.end]})
        result = intersect_intervals(self.chr_data,
                                     {interval.chr:search_query},
                                     suppreseChrNumberCheck=True)
        gene_idxs = result[interval.chr]["intersection"]
        return self.chr_data[interval.chr].iloc[gene_idxs,:]

     def get_binned_interval(self):
        logging.getLogger(__name__).error("Function not yet ready")
        raise Exception("Not ready")

     def delete_region(self):
        logging.getLogger(__name__).error("Function not yet ready")
        raise Exception("Not ready")