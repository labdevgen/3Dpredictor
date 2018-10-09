# Written by Minja, 09-2018
# Class to deal with Hi-C loops
# Stores Hi-C loops as data

import logging
import pandas as pd
from shared import FileReader
import numpy as np
import os


class LoopReader(FileReader):
    def __init__(self,fname,name=None):
        super(LoopReader,self).__init__(fname)
        if name is None:
            logging.getLogger(__name__).warning("Using filename as a name for predictor")
            self.name = os.path.basename(self.fname)
        else:
            self.name = name

    def process_data(self,data):
        #check duplicats
        duplicated = data.duplicated(subset=["chr1","x1","x2","chr2","y1","y2"])
        if sum(duplicated) > 0:
            logging.getLogger(__name__).warning(
                "Duplicates by genomic positions found in file " + self.fname)
        data.drop_duplicates(
            inplace=True)
        del duplicated

        # check inter-chromosomal loops
        assert len(data.query("chr1 != chr2"))==0

        # Set unique id for each loop
        data["id"] = list(range(len(data)))

        # Convert to chr-dict
        chr_data = dict([(chr,data[data["chr1"]==chr]) \
                         for chr in pd.unique(data["chr1"])])

        #sort
        sorted_data_l = {}
        sorted_data_r = {}
        for chr,data in chr_data.items():
            sorted_data_l[chr] = data.sort_values(by=["chr1","x1","x2","y1","y2"])
            sorted_data_r[chr] = data.sort_values(by=["chr1","y1","y2","x1","x2"])
        del chr_data

        # check for nested intervals
        nested_intevals_count = 0 #TODO check how is it work
        print_example = True
        for data in sorted_data_l.values():
            data_shifted = data.shift()
            nested_x = [False] + (data["x1"][1:] - data_shifted["x2"][1:] > 0) & \
                (data["x1"][1:] - data_shifted["x2"][1:] < 0)
            nested_y = [False] + (data["y1"][1:] - data_shifted["y2"][1:] > 0) & \
                (data["y1"][1:] - data_shifted["y2"][1:] < 0)

            nested_intevals_count += sum(nested_x) + sum(nested_y)

            if print_example and sum(nested_x) > 0:
                logging.getLogger(__name__).debug("Nested intervals found. Examples: ")
                logging.getLogger(__name__).debug(data[1:][nested_x].head(1))
                print_example = False
            elif print_example and sum(nested_y) > 0:
                logging.getLogger(__name__).debug("Nested intervals found. Examples: ")
                logging.getLogger(__name__).debug(data[1:][nested_y].head(1))
                print_example = False

        if nested_intevals_count > 0:
            logging.getLogger(__name__).warning("Number of nested intervals: "+str(nested_intevals_count))


        return sorted_data_l,sorted_data_r


    def read_loops(self, start_chr_name_with_chr=True):
        logging.getLogger(__name__).info(msg="Reading Loops file "+self.fname)
        data = pd.read_csv(self.fname,sep="\t")

        if start_chr_name_with_chr:
            #check wheather chr names include "chr"
            already_correct = np.any(data["chr1"].apply(lambda x: x.startswith("chr")).values)
            if not already_correct:
                data["chr1"] = data["chr1"].apply(lambda x: "chr"+str(x))
                data["chr2"] = data["chr2"].apply(lambda x: "chr"+str(x))

        #save
        self.chr_data_l, self.chr_data_r = self.process_data(data)
        del data

    def empty_like(self):
        return pd.DataFrame(columns = list(self.chr_data_l.values())[0].columns)


    def getLoops(self,chr):
        if not chr in self.chr_data_l.keys():
            return {chr:self.empty_like()}
        return {chr:self.chr_data_l[chr]}

    def getLeftLoopAncors(self,chr):
        if not chr in self.chr_data_l.keys():
            return {chr:self.empty_like()}
        return {chr:self.chr_data_l[chr][["x1","x2","id"]].rename(columns={"x1":"start","x2":"end"})}

    def getRightLoopAncors(self,chr):
        if not chr in self.chr_data_r.keys():
            return {chr:self.empty_like()}
        return {chr:self.chr_data_r[chr][["y1","y2","id"]].rename(columns={"y1":"start","y2":"end"})}