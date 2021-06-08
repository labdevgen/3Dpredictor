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

     # Modifies data according to the relationship between the "start" and "end" of genes, transcription direction
     # (tss file) and the "start" and "end" of the deletion interval.
     def delete_region_RNA(self, interval, tss_file):
         st, en = self.get_interval(interval, return_ids=True)
         old_length = len(self.chr_data[interval.chr])
         tss_file = pd.read_csv(tss_file, encoding='utf-8', sep='\t')
         # Search for genes affected by deletion in the file with the transcription direction:
         condition = np.where(((tss_file.strand == 1) & (interval.end > tss_file.position) &
                                (interval.start < (tss_file.position + 2000 * tss_file.strand))) |
                               ((tss_file.strand == -1) & (interval.start < tss_file.position) &
                                (interval.end > (tss_file.position + 2000 * tss_file.strand))))
         # Finding the intersection of affected genes by gene name in file with the direction of transcription and
         # file with RNAseq data. Creating a list of their indices:
         drop_indices = list(self.chr_data[interval.chr].index[np.in1d(self.chr_data[interval.chr].gene,
                                                                       tss_file.gene.iloc[condition])])
         debug = -len(drop_indices)
         # Change coordinates of other genes:
         self.chr_data[interval.chr].iloc[en + 1:,[self.chr_data[interval.chr].columns.get_loc("start"),
                                                   self.chr_data[interval.chr].columns.get_loc("end"),
                                                   self.chr_data[interval.chr].columns.get_loc("mids")]] -= interval.len
         # Delete genes affected by deletion:
         if len(drop_indices) > 0:
             self.chr_data[interval.chr].drop(drop_indices, inplace=True)
         # Set new indices according new "start" and "end":
         self.chr_data[interval.chr].set_index(self.chr_data[interval.chr].apply(
             lambda x: pd.Interval(x.start, x.end, closed="both"), axis="columns"), inplace=True)
         assert len(self.chr_data[interval.chr]) - debug == old_length

     # Modifies data according to the relationship between the "start" and "end" of genes, transcription direction
     # (tss file) and the "start" and "end" of the duplication interval.
     def duplicate_region_RNA(self, interval, tss_file):
         st, en = self.get_interval(interval, return_ids=True)
         old_length = len(self.chr_data[interval.chr])
         tss_file = pd.read_csv(tss_file, encoding='utf-8', sep='\t')
         # Search for genes affected by duplication in the file with the transcription direction:
         condition = np.where(
             (((tss_file.position + 2000 * tss_file.strand) > interval.end) &
              (tss_file.position < interval.end) & (tss_file.strand == 1)) |
             (((tss_file.position + 2000 * tss_file.strand) < interval.end) &
              (tss_file.position > interval.end) & (tss_file.strand == -1)))
         # List of genes indices to be duplicated:
         dup_indices = list(self.chr_data[interval.chr].index[(self.chr_data[interval.chr].start > interval.start) &
                                                              (self.chr_data[interval.chr].end < interval.end)])
         # Finding the intersection of affected genes by gene name in file with the direction of transcription and
         # file with RNAseq data. Creating a list of their indices:
         drop_indices = list(self.chr_data[interval.chr].index[np.in1d(self.chr_data[interval.chr].gene,
                                                                       tss_file.gene.iloc[condition])])
         debug = len(dup_indices) - len(drop_indices)
         dup_data = self.chr_data[interval.chr].loc[dup_indices]    # Duplicated genes as df
         dup_data[["mids", "start", "end"]] += interval.len         # Change coordinates of duplicated genes
         # Change coordinates of other genes:
         self.chr_data[interval.chr].iloc[en + 1:, [self.chr_data[interval.chr].columns.get_loc("start"),
                                                    self.chr_data[interval.chr].columns.get_loc("end"),
                                                    self.chr_data[interval.chr].columns.get_loc("mids")]] += interval.len
         # Delete genes affected by duplication:
         if len(drop_indices) > 0:
             self.chr_data[interval.chr].drop(drop_indices, inplace=True)
         # Add duplicated genes:
         self.chr_data[interval.chr] = pd.concat([self.chr_data[interval.chr], dup_data])
         self.chr_data[interval.chr].sort_values(by="start", inplace=True)
         # Set new indices according new "start" and "end":
         self.chr_data[interval.chr].set_index(self.chr_data[interval.chr].apply(
             lambda x: pd.Interval(x.start, x.end, closed="both"), axis="columns"), inplace=True)
         assert len(self.chr_data[interval.chr]) - debug == old_length
