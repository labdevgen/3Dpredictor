import ChiPSeqReader
import pandas as pd
import numpy as np
import logging

#A reader for RNA seq data
#It should store RNAseq data as Pandas DataFrame:
#Chr -- Start -- End -- Value
#The difference with ChipSeq data is that
#Genes can be long, so we cannot convert whole gene
#To one point called "mid", as we did with ChipSeq

class RNAseqReader(ChiPSeqReader):
     #Init is inhereted from ChiPSeqReader
     #It will set self.data to None
     #And take care about self.fname and self.name

     #read file
     #rename, if not None, passed to dataframe.rename(columns=rename)
     #As a results, dataframe should get column names
     #chr start end  sigVal gene
     #If rename==None dataframe.rename won't be called
     #args and kwargs are passed to pandas readcsv
     def read_file(self,rename=None,*args,**kwargs):
        #Read raw data
        data = pd.read_csv(self.fname,args,kwargs)

        #Rename if needed
        if rename != None:
            data.rename(columns=rename,inplace=True)

        #Drop unused fields
        data = data[["chr","start","end","sigVal","gene"]]

        #Check end > start for all genes
        assert np.all(data["end"].values - data["start"] > 0)

        # check duplicates, set mids, and split by chromosomes
        self.chr_data = self.process_data(data)
        del data

   def get_nearest_peaks():
        logging.getLogger(__name__).error("Function not yet ready")
        raise Exception("Not ready")

    def get_interval():
        logging.getLogger(__name__).error("Function not yet ready")
        raise Exception("Not ready")

    def get_binned_interval():
        logging.getLogger(__name__).error("Function not yet ready")
        raise Exception("Not ready")

    def delete_region():
        logging.getLogger(__name__).error("Function not yet ready")
        raise Exception("Not ready")