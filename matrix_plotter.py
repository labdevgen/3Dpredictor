import matplotlib.pyplot as plt
import pandas as pd
import logging
import numpy as np

class MatrixPlotter(): #A class that plot fragments of heatmap
                        #Note - not optimized for large datasets such as whole chromosome
    def process(self,data):
        return data.loc[:,["chr","contact_st","contact_en","contact_count"]]

    def set_control(self,ctrl):
        self.control = self.process(ctrl)

    def read_control(self,ctrl_fname):
        self.control = self.process(pd.read_csv(ctrl_fname,delimiter="\t"))

    def set_data(self,data):
        self.data = self.process(data)

    def read_data(self,data_fname):
        self.data = self.process(pd.read_csv(data_fname,delimiter="\t"))

    def convert2binned(self, data, interval, binsize):
        data["contact_st_bin"] = data["contact_st"].apply(lambda x: (x - interval.start) // binsize)
        data["contact_en_bin"] = data["contact_en"].apply(lambda x: (x - interval.start) // binsize)
        return data

    def getMatrix4plot(self, interval, binsize = None):
        def appendSeries2matrix(x,matrix,triangle = "both"):
            if triangle == "both":
                matrix[x["contact_en_bin"], x["contact_st_bin"]] = x["contact_count"]
                matrix[x["contact_st_bin"], x["contact_en_bin"]] = x["contact_count"]
            elif triangle == "upper":
                matrix[x["contact_st_bin"], x["contact_en_bin"]] = x["contact_count"]
            elif triangle == "lower":
                matrix[x["contact_en_bin"], x["contact_st_bin"]] = x["contact_count"]
            else:
                raise

        try:
            self.data
        except:
            logging.error("Please provide data firts")
            return

        if binsize == None: #infer binsize from data
            dist = pd.unique(self.data["contact_en"]-self.data["contact_st"])
            dist = dist[np.nonzero(dist)]
            assert len(dist) > 0
            binsize = min(dist)
            logging.info("Using binsize "+str(binsize))

        Interval_size_bins = (interval.end - interval.start) // binsize + 1
        matrix = np.zeros(shape = (Interval_size_bins , Interval_size_bins ))
        data = self.data.query("@interval.start <= contact_st <= @interval.end &"
                               "@interval.start <= contact_en <= @interval.end")
        data = self.convert2binned(data, interval, binsize)

        try:
            len(self.control)
            with_control = True
        except:
            with_control = False

        if with_control:
            logging.debug("Running with control")
            control = self.control.query("@interval.start <= contact_st <= @interval.end &"
                                   "@interval.start <= contact_en <= @interval.end")
            control = self.convert2binned(control, interval, binsize)
            logging.debug(control.head())
            logging.debug(data.head())
            data.apply(appendSeries2matrix,matrix=matrix,triangle="upper",axis = "columns")
            control.apply(appendSeries2matrix,matrix=matrix,triangle="lower",axis = "columns")
        else:
            data.apply(appendSeries2matrix,matrix=matrix,triangle="both",axis = "columns")

        return matrix
