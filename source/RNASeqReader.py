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

     def process_data(self,data):
        sorted_data = super(RNAseqReader, self).process_data(data)
        # set IntervalIndex used by intersect_intervals function
        for chr in sorted_data.keys():
            df = sorted_data[chr]
            sorted_data[chr] = df.set_index(df.apply(
                    lambda x: pd.Interval(x.start, x.end, closed="both"),
                    axis="columns"))
        return sorted_data

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

     def get_interval(self, interval, return_ids=False):
                                       # Return all genes that intersect interval
                                       # Also counts partial intersections
        return intersect_with_interval(self.chr_data,interval, return_ids=return_ids)


     def get_binned_interval(self):
        logging.getLogger(__name__).error("Function not yet ready")
        raise Exception("Not ready")

     def delete_region(self, interval):
         # TODO
         # this won't work for long genes. Imagine, we are deleting regions about 40 base pairs from inside the gene
         # then, we delete all this gene however usually in cells transcription is remained in this region.
         # We should use TSS snd TSE for accurate prediction!
         # raise NotImplementedError
         data = self.chr_data[interval.chr]
         st, en = self.get_interval(interval, return_ids=True)
         debug = len(self.get_interval(interval))
         self.chr_data[interval.chr].iloc[en+1:, data.columns.get_loc("start")] -= interval.len
         self.chr_data[interval.chr].iloc[en+1:, data.columns.get_loc("end")] -= interval.len
         self.chr_data[interval.chr].iloc[en+1:, data.columns.get_loc("mids")] -= interval.len
         old_length = len(self.chr_data[interval.chr])
         if st != -1:
            data.index = data.index.map(str)
            data.drop(index=data.index[st:en+1], inplace=True)
            self.chr_data[interval.chr] = data.set_index(data.apply(
                lambda x: pd.Interval(x.start, x.end, closed="both"),
                axis="columns"))
         assert len(self.chr_data[interval.chr]) + debug == old_length

     def duplicate_region_RNA(self, interval):                  # Modifies data according to the relationship
         st, en = self.get_interval(interval, return_ids=True)  # between the "start" and "end" of genes, transcription direction and
         old_length = len(self.chr_data[interval.chr])          # the "start" and "end" of the duplication interval
         tss_strand = pd.read_csv(r"C:\Users\Maria\PycharmProjects\example\tss_strand.csv", encoding='utf-8', sep='\t')
         condition = np.where(                                                     # Search for genes affected by
             (((tss_strand.position + 2000 * tss_strand.strand) > interval.end) &  # duplication in the file with
              (tss_strand.position < interval.end) & (tss_strand.strand == 1)) |   # the transcription direction
             (((tss_strand.position + 2000 * tss_strand.strand) < interval.end) &
              (tss_strand.position > interval.end) & (tss_strand.strand == -1)))
         dup_indices = list(self.chr_data[interval.chr].index[(self.chr_data[interval.chr].start > interval.start) & # List of gene indices
                                                              (self.chr_data[interval.chr].end < interval.end)])     # to be duplicated
         drop_indices = list(self.chr_data[interval.chr].index[np.in1d(self.chr_data[interval.chr].gene,   # Finding the intersection of affected genes in files with
                                                                       tss_strand.gene.iloc[condition])])  # the direction of transcription and with RNAseq data.
         debug = len(dup_indices) - len(drop_indices)                                                      # Creating a list from their indices.
         dup_data = self.chr_data[interval.chr].loc[dup_indices]    # Duplicated genes as df
         dup_data[["mids", "start", "end"]] += interval.len         # Change coordinates of duplicated genes
         self.chr_data[interval.chr].iloc[en + 1:, [self.chr_data[interval.chr].columns.get_loc("start"),                # Change coordinates
                                                    self.chr_data[interval.chr].columns.get_loc("end"),                  # of other genes
                                                    self.chr_data[interval.chr].columns.get_loc("mids")]] += interval.len
         if len(drop_indices) > 0:
             self.chr_data[interval.chr].drop(drop_indices, inplace=True)   # Delete genes affected by duplication
         self.chr_data[interval.chr] = pd.concat([self.chr_data[interval.chr], dup_data]) # # Adds duplicated genes
         self.chr_data[interval.chr].sort_values(by="start", inplace=True)
         self.chr_data[interval.chr].set_index(self.chr_data[interval.chr].apply(
             lambda x: pd.Interval(x.start, x.end, closed="both"), axis="columns"), inplace=True)  # Set new indices according new "start" and "end"
         assert len(self.chr_data[interval.chr]) - debug == old_length
