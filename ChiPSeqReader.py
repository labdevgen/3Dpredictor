import logging,os
from shared import Position,FileReader
import pandas as pd
import numpy as np


class ChiPSeqReader(FileReader): #Class process files with ChipSeq peaks
    def __init__(self,fname,name=None):
        self.data = None
        if name == None:
            logging.warning("Using filename as a name for predictor")
            self.proteinName = os.path.basename(fname)
        else:
            self.proteinName = name

        super(ChiPSeqReader,self).__init__(fname)

    def read_file(self): # store CTCF peaks as sorted pandas dataframe
        logging.log(msg="Reading CTCF file "+self.fname, level=logging.INFO)

        # set random temporary labels
        Nfields = len(open(self.fname).readline().strip().split())
        names = list(map(str,list(range(Nfields))))
        data = pd.read_csv(self.fname,sep="\t",header=None,names=names)

        # subset and rename
        data = data.iloc[:,[0,1,2,6]]
        data.rename(columns={"0":"chr","1":"start","2":"end","6":"sigVal"},
                    inplace=True)

        #check duplicats
        duplicated = data.duplicated(subset = ["chr","start","end"])
        if sum(duplicated) > 0:
            logging.warning("Duplicates by genomic positions found in file "+self.fname) #TODO check why this happens
        data.drop_duplicates(inplace=True) #Keep peaks with same coordinate and different sigVal, if such peask exist
        del duplicated

        #get peak mids
        data["mids"] = (data["start"] + data["end"]) // 2

        #convert to chr-dict
        chr_data = dict([(chr,data[data["chr"]==chr]) \
                         for chr in pd.unique(data["chr"])])
        del data

        logging.warning(
            msg="Filling orientation with mock values!") #TODO fill with real data

        for data in chr_data.values():
            data.sort_values(by=["chr","start"],inplace=True)
            data["orientation"] = [1]*len(data)

        #save
        self.chr_data = chr_data

    def get_nearest_peaks(self,point,side,N,N_is_strict=True):
        def get_mock_data(N,midpos):
            data = pd.DataFrame(columns=self.chr_data[point.chr].columns)
            data["mids"] = [midpos]*N
            data["sigVal"] = [0]*N
            return data

        #Some checks removed to speed up proccess
        #assert point.start == point.end
        #assert point.chr in self.chr_data.keys()
        #try:
        #    self.chr_data
        #except:
        #    logging.error("Please read data first")
        #    return None

        search_id = np.searchsorted(self.chr_data[point.chr]["mids"].values,point.start,side)
        if side == "left":
            result = self.chr_data[point.chr].iloc[max(search_id-N,0):search_id,:]
            if len(result) == 0:
                return get_mock_data(N,midpos=0)
            if len(result) < N and N_is_strict:
                return result.append(get_mock_data(N - len(result),midpos=0))
            return result
        elif side == "right":
            if search_id == len(self.chr_data[point.chr]):
                return get_mock_data(N=N,midpos=(point.start + 1))
            result = self.chr_data[point.chr].iloc[search_id:min(search_id+N,
                                                                 len(self.chr_data[point.chr])),:]
            if len(result) < N and N_is_strict:
                return result.append(get_mock_data(N=(N - len(result)),
                                     midpos=self.chr_data[point.chr]["mids"].iat[-1]))
            return result
        else:
            raise Exception()

    def get_interval(self,interval,return_ids=False): #return all peaks intersecting interval as pd.dataframen
                                                    #if return_ids=True returns int positions of first and last peak
        try:
            self.chr_data
        except:
            logging.error("Please read data first")
            return None

        if not interval.chr in self.chr_data:
            logging.error("Chromosome ",interval.chr,"not found in keys")
            return pd.DataFrame()

        left = np.searchsorted(self.chr_data[interval.chr]["mids"].values,interval.start,"left")
        right = np.searchsorted(self.chr_data[interval.chr]["mids"].values,interval.end,"right")
        assert left <= right
        if not return_ids:
            return self.chr_data[interval.chr].iloc[min(len(self.chr_data[interval.chr])-1,left):max(0,right),:]
        else:
            return min(len(self.chr_data[interval.chr]) - 1, left),max(0, right)

    def get_binned_interval(self,interval,binsize,extend=True): #return all peaks intersecting interval as list
                                                            #list len = interval len / bindsize
                                                            #if extend = True last bin is extended over interval's end
        try:
            self.chr_data
        except:
            logging.error("Please read data first")
            return None

        if not interval.chr in self.chr_data:
            logging.error("Chromosome ",interval.chr,"not found in keys")
            return pd.DataFrame()

        left_array = np.arange(interval.start,interval.end-1,binsize)
        right_array = left_array + binsize

        if (right_array[-1] >= interval.end) and (not extend):
            right_array[-1] = right_array.end+1

        left_ids = np.searchsorted(self.chr_data[interval.chr]["mids"].values,left_array,"left")
        right_ids = np.searchsorted(self.chr_data[interval.chr]["mids"].values,right_array,"right")

        result_strength = [0]*len(left_ids)

        for ind,(left,right) in enumerate(zip(left_ids,right_ids)):
            left = min(len(self.chr_data[interval.chr])-1,left)
            right = max(0,right)
            assert left <= right
            if left != right:
                result_strength[ind] = self.chr_data[interval.chr].sigVal.iloc[left:right].sum()
        return result_strength

    def delete_region(self,interval):
        debug = len(self.get_interval(interval))
        data = self.chr_data[interval.chr]
        st,en = self.get_interval(interval,return_ids=True)
        logging.debug(self.chr_data[interval.chr].iloc[data.columns.get_loc("start"),en:].head())
        self.chr_data[interval.chr].iloc[data.columns.get_loc("start"),en:] -= interval.len
        self.chr_data[interval.chr].iloc[data.columns.get_loc("end"),en:] -= interval.len
        self.chr_data[interval.chr].iloc[data.columns.get_loc("mids"),en:] -= interval.len
        old_length = len(self.chr_data[interval.chr])
        self.chr_data[interval.chr].drop(data.index[st:en],inplace=True)
        assert len(self.chr_data[interval.chr]) + debug == old_length