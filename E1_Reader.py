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
        self.warnings = {"Interval length < binsize":0,"End of the interval over the chromsome end":0}
        self.deleted_leftover = 0

    def get_binsize(self):
        return self.binsize

    def read_file(self,fname,chrName):
        if self.binsize == 0:
            logging.error("Cannot load file "+fname+" : bin size not set")
        E1 = pd.read_csv(fname,header=None,names=["E1"])
        E1.fillna(value=0,inplace=True)
        E1["mid"] = np.arange(len(E1))*self.binsize + (self.binsize // 2)
        E1["chr"] = [chrName]*len(E1)
        if chrName in self.data:
            logging.getLogger(__name__).warning("Overwriting data for chr"+chrName)
        self.data[chrName] = E1

    def read_files(self,fnames,**kwargs):
        if "binsize" in kwargs:
            self.binsize = kwargs["binsize"]
        for file in fnames:
            if "binSizeFromName" in kwargs:
                binsize = kwargs["binSizeFromName"](file)
                if self.binsize != 0 and self.binsize != binsize:
                    logging.getLogger(__name__).warning("Binsize in file "+file
                                        +" does not match binsize "+str(binsize))
                    raise
                self.binsize = binsize
            chrName = os.path.basename(file).split(".")[0]
            self.read_file(file,chrName)

    def interval2ids(self,interval,verbose):
        start_id = (interval.start // self.binsize)
        data = self.data[interval.chr]
        if start_id >= len(data):
            logging.error("Begining of the interval over the chromsome end")
            logging.getLogger(__name__).debug(str(interval))
            logging.getLogger(__name__).debug("start id "+str("start_id"))
            logging.getLogger(__name__).debug("binsize "+str(self.binsize))
            raise
        end_id = start_id + ((interval.end - interval.start) // self.binsize)
        if end_id == start_id:
            assert (interval.end - interval.start) < self.binsize
            logging.log(msg="Interval length < binsize",level = verbose)
            self.warnings["Interval length < binsize"] += 1
            end_id += 1
        return start_id,end_id


    def get_E1inInterval(self,interval,make_consistent_bins=True,verbose = logging.NOTSET):
        #How it works:
        #Find the beginning of the interval
        #Add a number of bins equal to interval_size // binsize
        start_id,end_id = self.interval2ids(interval,verbose)
        data = self.data[interval.chr]

        if end_id > len(data):
            logging.log(msg="End of the interval over the chromsome end",
                        level=verbose)
            logging.log(msg=str(interval),level=verbose)
            self.warnings["End of the interval over the chromsome end"] += 1
            #logging.getLogger(__name__).debug("start id "+str("start_id"))
            #logging.getLogger(__name__).debug("end id "+str("end_id"))
            #logging.getLogger(__name__).debug("binsize "+str(self.binsize))
            result = data.iloc[start_id:len(data),:]

            if make_consistent_bins:
                #add mock data
                mock_len = (end_id - start_id - len(result))
                if mock_len <= 0:
                    logging.getLogger(__name__).debug("Mock len = "+str(mock_len))
                    logging.getLogger(__name__).debug(" ".join(["Start",str(start_id),"End",str(end_id),"Len",
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

    def print_warnings(self,verbose=logging.INFO,reset=True):
        if sum(self.warnings.values()) > 0:
            for (key, val) in self.warnings.items():
                logging.log(msg=key+" warning occured "+str(val)+" times",
                        level=verbose)
            if reset:
                self.warnings = {"Interval length < binsize": 0, "End of the interval over the chromsome end": 0}


    def delete_region(self,interval,verbose = logging.NOTSET):
        start_id,end_id = self.interval2ids(interval,verbose)
        data = self.data[interval.chr]
        self.deleted_leftover += interval.len - (end_id-start_id)*self.binsize
        if self.deleted_leftover > self.binsize:
            logging.error("Multuple deletion are not yet implemented")
            raise Exception()
        self.data[interval.chr].drop(self.data[interval.chr].index[start_id:end_id],inplace=True)