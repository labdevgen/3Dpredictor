import logging
import pandas as pd
import numpy as np
import os

def fileName2binsize(name):
    last = name[-1]
    multiply = 1
    if last.lower() == "k":
        multiply = 1000
    elif last.lower() == "m":
        multiply = 1000000
    if multiply != 1:
        return int(os.path.basename(name.split(".")[-1][:-1]))*multiply
    else:
        return int(os.path.basename(name.split(".")[-1]))


class E1Reader():
    def __init__(self,binsize=0):
        self.binsize = binsize
        self.data = {}

    def read_file(self,fname,chrName):
        if self.binsize == 0:
            logging.error("Cannot load file "+fname+" : bin size not set")
        E1 = pd.read_csv(fname,header=None,names=["E1"])
        E1.fillna(value=0,inplace=True)
        E1["mid"] = np.arange(len(E1))*self.binsize + (self.binsize // 2)
        E1["chr"] = [chrName]*len(E1)
        if chrName in self.data:
            logging.warning("Overwriting data for chr"+chrName)
        self.data[chrName] = E1

    def read_files(self,fnames,**kwargs):
        if "binsize" in kwargs:
            self.binsize = kwargs["binsize"]
        for file in fnames:
            if "binSizeFromName" in kwargs:
                binsize = kwargs["binSizeFromName"](file)
                if self.binsize != 0 and self.binsize != binsize:
                    logging.warning("Binsize in file "+file
                                        +" does not match binsize "+str(binsize))
                    raise
                self.binsize = binsize
            chrName = os.path.basename(file).split(".")[0]
            self.read_file(file,chrName)

    def get_E1inInterval(self,interval,make_consistent_bins=True):
        #How it works:
        #Find the beginning of the interval
        #Add a number of bins equal to interval_size // binsize
        start_id = (interval.start // self.binsize)
        data = self.data[interval.chr]
        if start_id >= len(data):
            logging.error("Begining of the interval over the chromsome end")
            logging.debug(str(interval))
            logging.debug("start id "+str("start_id"))
            logging.debug("binsize "+str(self.binsize))
            raise
        end_id = start_id + ((interval.end - interval.start) // self.binsize)
        if end_id == start_id:
            assert (interval.end - interval.start) < self.binsize
            logging.warning("Interval length < binsize")
            end_id += 1

        if end_id > len(data):
            logging.warning("End of the interval over the chromsome end")
            logging.debug(str(interval))
            logging.debug("start id "+str("start_id"))
            logging.debug("end id "+str("end_id"))
            logging.debug("binsize "+str(self.binsize))
            result = data.iloc[start_id:len(data),:]

            if make_consistent_bins:
                #add mock data
                mock_len = (end_id - start_id - len(result))
                if mock_len <= 0:
                    logging.debug("Mock len = "+str(mock_len))
                    logging.debug(" ".join(["Start",str(start_id),"End",str(end_id),"Len",
                                           str(len(data)),str(interval.len // self.binsize)]))
                assert mock_len > 0
                mock = pd.DataFrame(columns=result.columns)
                mock["chr"] = [interval.chr]*mock_len
                mock["mid"] = (np.arange(mock_len) + 1)*self.binsize + result["mid"].iat[-1]
                mock["E1"] = [0]*mock_len
                result = result.append(mock)

            return result
        else:
            return data[start_id:end_id]
