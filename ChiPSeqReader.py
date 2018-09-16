import logging,os
from shared import Position,FileReader,intersect_intervals
import pandas as pd
import numpy as np



class ChiPSeqReader(FileReader): #Class process files with ChipSeq peaks
    def __init__(self,fname,name=None):
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
        duplicated = data.duplicated(subset=["chr", "start", "end"])
        if sum(duplicated) > 0:
            logging.warning(
                "Duplicates by genomic positions found in file " + self.fname)  # TODO check why this happens
        data.drop_duplicates(
            inplace=True)
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
            data["plus_orientation"] = [0]*len(data)
            data["minus_orientation"] = [0]*len(data)

        #save
        self.chr_data = chr_data

    def get_nearest_peaks(self,point,side,N,N_is_strict=True):
        def get_mock_data(N,midpos):
            data = pd.DataFrame(columns=self.chr_data[point.chr].columns)
            data["mids"] = [midpos]*N
            data["sigVal"] = [0]*N
            data["chr"] = [point.chr]*N
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
            if len(result) == 0: #TODO is it possible?
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

    def get_N_nearest_peaks_in_interval(self,point,N,N_is_strict=True ):
        def get_mock_data(N,midpos):
            data = pd.DataFrame(columns=self.chr_data[point.chr].columns)
            data["mids"] = [midpos]*N
            data["sigVal"] = [0]*N
            data["chr"] = [point.chr]*N
            return data


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

    def read_orient_file(self):  # store peaks with orientation from gimmeMotifs as sorted pandas dataframe
        logging.log(msg="Reading CTCF_orientation file " + self.fname, level=logging.INFO)

        # set random temporary labels
        Nfields = len(open(self.fname).readline().strip().split())
        names = list(map(str, list(range(Nfields))))
        data = pd.read_csv(self.fname, sep="\t", header=None, names=names)  #TODO check: CTCF fname == CTCF orient fname.split('-')[0]
        # print(data)

        # subset and rename
        data = data.iloc[:, [0, 1, 2, 4, 5]]
        data.rename(columns={"0": "chr", "1": "start", "2": "end", "4": "score", "5": "orientation"},
                    inplace=True)

        # check duplicats
        duplicated = data.duplicated(subset=["chr", "start", "end"])
        if sum(duplicated) > 0:
            logging.warning(
                "Duplicates by genomic positions found in file " + orient_reader.fname)  # TODO check why this happens
        data.drop_duplicates(
            inplace=True)
        del duplicated

        # convert to chr-dict
        chr_data = dict([(chr, data[data["chr"] == chr]) \
                         for chr in pd.unique(data["chr"])])
        del data

        for data in chr_data.values():
            data.sort_values(by=["chr", "start"], inplace=True)

        # save
        orient_reader.chr_data = chr_data


    def set_sites_orientation(self, orient_fname):
        try:
            self.chr_data
        except:
            logging.error("Please read data first")
            return None
        orient_chr_data = self.read_orient_file(orient_fname)
        result_intersection_data = intersect_intervals(self.chr_data, orient_chr_data)
        #print(result_intersection_data['chr1'])
        for chr in result_intersection_data.keys():
            #if chr != 'chr4':
             #   continue
            result_intersection_data[chr].sort_values(by=["orientation", "intersection", "score"], inplace=True)
            #print(result_intersection_data['chr1'])
            #duplicated = result_intersection_data[chr].duplicated(subset=["intersection", "orientation"])
            #print(duplicated)
            result_intersection_data[chr].drop_duplicates(subset=["intersection", "orientation"], keep="first", inplace=True)
            #print(result_intersection_data['chr1'])
            plus_orient_data = result_intersection_data[chr].query("orientation =='+'")
            plus_col_ind = self.chr_data[chr].columns.get_loc("plus_orientation")
            plus_row_list = list(plus_orient_data["intersection"])
            self.chr_data[chr].iloc[plus_row_list, plus_col_ind] = list(plus_orient_data["score"])

            minus_orient_data = result_intersection_data[chr].query("orientation =='-'")
            minus_col_ind = self.chr_data[chr].columns.get_loc("minus_orientation")
            minus_row_list = list(minus_orient_data["intersection"])
            self.chr_data[chr].iloc[minus_row_list, minus_col_ind] = list(minus_orient_data["score"])

    def delete_region(self,interval):
        debug = len(self.get_interval(interval))
        data = self.chr_data[interval.chr]
        st,en = self.get_interval(interval,return_ids=True)
        self.chr_data[interval.chr].iloc[en:,data.columns.get_loc("start")] -= interval.len
        self.chr_data[interval.chr].iloc[en:,data.columns.get_loc("end")] -= interval.len
        self.chr_data[interval.chr].iloc[en:,data.columns.get_loc("mids")] -= interval.len
        old_length = len(self.chr_data[interval.chr])
        self.chr_data[interval.chr].drop(data.index[st:en],inplace=True)
        assert len(self.chr_data[interval.chr]) + debug == old_length