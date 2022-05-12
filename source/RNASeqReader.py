# Written by Minja, 2018-09
# A reader for RNA seq data
# It should store RNAseq data as Pandas DataFrame:
# Chr -- Start -- End -- sigVal -- gene
# The difference with ChipSeq data is that
# Genes can be long, so we cannot convert whole gene
# To one point called "mid", as we did with ChipSeq

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
    # Init is inhereted from ChiPSeqReader
    # It will set self.data to None
    # And take care about self.fname and self.name

    def process_data(self, data):
        sorted_data = super(RNAseqReader, self).process_data(data)
        # set IntervalIndex used by intersect_intervals function
        for chr in sorted_data.keys():
            df = sorted_data[chr]
            sorted_data[chr] = df.set_index(df.apply(
                lambda x: pd.Interval(x.start, x.end, closed="both"),
                axis="columns"))
        return sorted_data

    # read file
    # rename, if not None, passed to dataframe.rename(columns=rename)
    # As a results, dataframe should get column names
    # chr start end  sigVal gene
    # If rename==None dataframe.rename won't be called
    # if start_chr_name_with_chr = True, check wheather chromosome names starts with chr
    # if not, add chr at the begningn
    # args and kwargs are passed to pandas readcsv
    def read_file(self, rename=None, start_chr_name_with_chr=True, *args, **kwargs):
        # Read raw data
        data = pd.read_csv(self.fname, *args, **kwargs)

        # Rename if needed
        if rename != None:
            data.rename(columns=rename, inplace=True)
        # Drop unused fields
        data = data[["chr", "start", "end", "sigVal", "gene"]]

        if start_chr_name_with_chr:
            # check wheather chr names include "chr"
            already_correct = np.any(data["chr"].apply(lambda x: x.startswith("chr")).values)
            if not already_correct:
                data["chr"] = data["chr"].apply(lambda x: "chr" + str(x))

        # Check end > start for all genes
        assert np.all(data["end"].values - data["start"] > 0)

        # check duplicates, set mids, and split by chromosomes and sort
        self.chr_data = self.process_data(data)
        del data

    # def get_nearest_peaks(...): - Should work well when inhereted from ChiPSeq reader

    # Depricated, use get_interval
    def _get_interval(self, interval):  # Return all genes that intersect interval
        # Also counts partial intersections
        search_query = pd.DataFrame({"start": [interval.start], "end": [interval.end]})
        result = intersect_intervals(self.chr_data,
                                     {interval.chr: search_query},
                                     suppreseChrNumberCheck=True)
        gene_idxs = result[interval.chr]["intersection"]
        return self.chr_data[interval.chr].iloc[gene_idxs, :]

    def get_interval(self, interval, return_ids=False):
        # Return all genes that intersect interval
        # Also counts partial intersections
        return intersect_with_interval(self.chr_data, interval, return_ids=return_ids)

    def get_binned_interval(self):
        logging.getLogger(__name__).error("Function not yet ready")
        raise Exception("Not ready")

    # Modifies data according to the relationship between the "start" and "end" of genes, transcription direction
    # (tss file) and the "start" and "end" of the deletion interval.
    def delete_region_RNA(self, interval, tss_file):
        st, en = self.get_interval(interval, return_ids=True)
        old_length = len(self.chr_data[interval.chr])
        tss_file = pd.read_csv(tss_file, encoding='utf-8', sep='\t')
        # Add a column with the direction of transcription by merging RNAseq data and tss file:
        self.chr_data[interval.chr] = pd.merge(self.chr_data[interval.chr],
                                               tss_file.loc[:, ['Strand', 'Gene stable ID']], how="left",
                                               left_on="gene", right_on="Gene stable ID")
        # Search for genes affected by deletion. Creating a list of their indices to be deleted:
        drop_indices = list(
            np.where(((self.chr_data[interval.chr].Strand == 1) & (interval.end > self.chr_data[interval.chr].start) &
                      (interval.start < (
                                  self.chr_data[interval.chr].start + 2000 * self.chr_data[interval.chr].Strand))) |
                     ((self.chr_data[interval.chr].Strand == -1) & (interval.start < self.chr_data[interval.chr].end) &
                      (interval.end > (self.chr_data[interval.chr].end + 2000 * self.chr_data[interval.chr].Strand))))[
                0])
        debug = -len(drop_indices)
        # Change coordinates of other genes:
        self.chr_data[interval.chr].iloc[en + 1:, [self.chr_data[interval.chr].columns.get_loc("start"),
                                                   self.chr_data[interval.chr].columns.get_loc("end"),
                                                   self.chr_data[interval.chr].columns.get_loc("mids")]] -= interval.len
        # Delete genes affected by deletion:
        if len(drop_indices) > 0:
            self.chr_data[interval.chr].drop(drop_indices, inplace=True)
        # Set new indices according new "start" and "end":
        self.chr_data[interval.chr].set_index(self.chr_data[interval.chr].apply(
            lambda x: pd.Interval(x.start, x.end, closed="both"), axis="columns"), inplace=True)
        self.chr_data[interval.chr].drop(['Strand', 'Gene stable ID'], axis=1, inplace=True)
        assert len(self.chr_data[interval.chr]) - debug == old_length

    # Modifies data according to the relationship between the "start" and "end" of genes, transcription direction
    # (tss file) and the "start" and "end" of the duplication interval.
    def duplicate_region_RNA(self, interval, tss_file):
        st, en = self.get_interval(interval, return_ids=True)
        old_length = len(self.chr_data[interval.chr])
        tss_file = pd.read_csv(tss_file, encoding='utf-8', sep='\t')
        # Add a column with the direction of transcription by merging RNAseq data and tss file:
        self.chr_data[interval.chr] = pd.merge(self.chr_data[interval.chr],
                                               tss_file.loc[:, ['Strand', 'Gene stable ID']], how="left",
                                               left_on="gene", right_on="Gene stable ID")
        # Search for genes affected by duplication. Creating a list of their indices to be deleted:
        drop_indices = list(np.where(
            (((self.chr_data[interval.chr].start + 2000 * self.chr_data[interval.chr].Strand) > interval.end) &
             (self.chr_data[interval.chr].start < interval.end) & (self.chr_data[interval.chr].Strand == 1)) |
            (((self.chr_data[interval.chr].end + 2000 * self.chr_data[interval.chr].Strand) < interval.end) &
             (self.chr_data[interval.chr].end > interval.end) & (self.chr_data[interval.chr].Strand == -1)))[0])
        # List of genes indices to be duplicated:
        dup_indices = list(self.chr_data[interval.chr].index[(self.chr_data[interval.chr].start > interval.start) &
                                                             (self.chr_data[interval.chr].end < interval.end)])
        debug = len(dup_indices) - len(drop_indices)
        dup_data = self.chr_data[interval.chr].loc[dup_indices]  # Duplicated genes as df
        dup_data[["mids", "start", "end"]] += interval.len  # Change coordinates of duplicated genes
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
        self.chr_data[interval.chr].drop(['Strand', 'Gene stable ID'], axis=1, inplace=True)
        assert len(self.chr_data[interval.chr]) - debug == old_length

    # Modifies data according to the relationship between the "start" and "end" of genes, transcription direction
    # (tss file) and the "start" and "end" of the duplication interval.
    def inverse_region_RNA(self, interval, tss_file):
        st, en = self.get_interval(interval, return_ids=True)
        old_length = len(self.chr_data[interval.chr])
        tss_file = pd.read_csv(tss_file, encoding='utf-8', sep='\t')
        # Добавить колонку с направлением транскрипции путем слияния данных RNAseq и tss файла:
        self.chr_data[interval.chr] = pd.merge(self.chr_data[interval.chr],
                                               tss_file.loc[:, ['Strand', 'Gene stable ID']], how="left",
                                               left_on="gene", right_on="Gene stable ID")
        # Поиск генов, затронутых инверсией. Создание списка их индексов:
        drop_indices = list(np.where(
            (((self.chr_data[interval.chr].Strand == 1) & ((self.chr_data[interval.chr].start < interval.start) &
                                                           ((self.chr_data[interval.chr].start + 2000 * self.chr_data[
                                                               interval.chr].Strand) > interval.start)) |
              ((self.chr_data[interval.chr].start < interval.end) &
               ((self.chr_data[interval.chr].start + 2000 * self.chr_data[interval.chr].Strand) > interval.end)))) |
            (((self.chr_data[interval.chr].Strand == -1) & ((self.chr_data[interval.chr].end > interval.start) &
                                                            ((self.chr_data[interval.chr].end + 2000 * self.chr_data[
                                                                interval.chr].Strand) < interval.start)) |
              ((self.chr_data[interval.chr].end > interval.end) &
               ((self.chr_data[interval.chr].end + 2000 * self.chr_data[interval.chr].Strand) < interval.end)))))[0])

        debug = len(self.chr_data[interval.chr])
        # новые координаты для инвертируемых генов
        starts = [interval.start + (interval.end - x) if x > (interval.start + interval.len / 2)
                  else interval.end - (x - interval.start) for x in
                  self.chr_data[interval.chr].iloc[st:en + 1, self.chr_data[interval.chr].columns.get_loc("end")]]
        # print(starts, 'starts')
        ends = [interval.start + (interval.end - x) if x > (interval.start + interval.len / 2)
                else interval.end - (x - interval.start) for x in
                self.chr_data[interval.chr].iloc[st:en + 1, self.chr_data[interval.chr].columns.get_loc("start")]]
        # print(ends, 'ends')
        # exit()
        if len(starts) > 0:
            self.chr_data[interval.chr].iloc[st:en + 1, self.chr_data[interval.chr].columns.get_loc("start")] = starts
        if len(ends) > 0:
            self.chr_data[interval.chr].iloc[st:en + 1, self.chr_data[interval.chr].columns.get_loc("end")] = ends

        self.chr_data[interval.chr].sort_values(by="start", inplace=True)

        # Delete genes affected by inversion:
        if len(drop_indices) > 0:
            self.chr_data[interval.chr].drop(drop_indices, inplace=True)

        # Set new indices according new "start" and "end":
        self.chr_data[interval.chr].set_index(self.chr_data[interval.chr].apply(
            lambda x: pd.Interval(x.start, x.end, closed="both"), axis="columns"), inplace=True)
        self.chr_data[interval.chr].drop(['Strand', 'Gene stable ID'], axis=1, inplace=True)
        assert old_length - debug == 0
