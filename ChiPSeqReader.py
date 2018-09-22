import logging,os
from shared import Position,FileReader,intersect_intervals
import pandas as pd
import numpy as np


class ChiPSeqReader(FileReader): #Class process files with ChipSeq peaks
    def __init__(self,fname,name=None):
        self.data = None
        if name == None:
            logging.getLogger(__name__).warning("Using filename as a name for predictor")
            self.proteinName = os.path.basename(fname)
        else:
            self.proteinName = name
        self.orient_data_real = False
        self.only_orient_peaks = False
        super(ChiPSeqReader,self).__init__(fname)

    #check duplicates, set mids, and split by chromosomes and sort
    def process_data(self,data):
        #check duplicats
        duplicated = data.duplicated(subset=["chr", "start", "end"])
        if sum(duplicated) > 0:
            logging.getLogger(__name__).warning(
                "Duplicates by genomic positions found in file " + self.fname)  # TODO check why this happens
        data.drop_duplicates(
            inplace=True)
        del duplicated
        #get peak mids
        data["mids"] = (data["start"] + data["end"]) // 2

        #convert to chr-dict
        chr_data = dict([(chr,data[data["chr"]==chr]) \
                         for chr in pd.unique(data["chr"])])

        #sort
        sorted_data = {}
        for chr,data in chr_data.items():
            sorted_data[chr] = data.sort_values(by=["chr","start"])
        del chr_data

        return sorted_data

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

        chr_data = self.process_data(data)
        del data

        logging.getLogger(__name__).warning(
            msg="Filling orientation with mock values!") #TODO fill with real data

        #Fill orientation with mock values
        for data in chr_data.values():
            data["plus_orientation"] = [0]*len(data)
            data["minus_orientation"] = [0]*len(data)

        #save
        self.chr_data = chr_data
        self.orient_data_real = False
        self.only_orient_peaks = False

    def get_mock_data(self, N,interval, midpos):
        data = pd.DataFrame(columns=self.chr_data[interval.chr].columns)
        data["mids"] = [midpos] * N
        data["sigVal"] = [0] * N
        data["chr"] = [interval.chr] * N
        data["plus_orientation"] = [0] * N
        data["minus_orientation"] = [0] * N
        return data

    def get_nearest_peaks(self,point,side,N,N_is_strict=True):
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
                return self.get_mock_data(N=N - len(result),interval=point, midpos=0) #TODO why midpos is 0?
            if len(result) < N and N_is_strict:
                return result.append(self.get_mock_data(N=N - len(result),interval=point, midpos=0))
            return result
        elif side == "right":
            if search_id == len(self.chr_data[point.chr]):
                return self.get_mock_data(N=N,interval=point, midpos=(point.start + 1))
            result = self.chr_data[point.chr].iloc[search_id:min(search_id+N,
                                                                 len(self.chr_data[point.chr])),:]
            if len(result) < N and N_is_strict:
                return result.append(self.get_mock_data(N=(N - len(result)),interval=point,
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

    def get_N_peaks_near_interval_boundaries(self,interval,N,N_is_strict=True ):
        all_peaks_in_interval = self.get_interval(interval)
        if not self.orient_data_real:
               logging.error("please set_orientation first")
        if len(all_peaks_in_interval) >= N*2:
            result_right = all_peaks_in_interval.head(N)
            result_left = all_peaks_in_interval.tail(N)
            return [result_right, result_left]
        elif 0 < len(all_peaks_in_interval) < N*2:
            if len(all_peaks_in_interval) == 1 and N_is_strict:
                result_right = self.get_mock_data(N=N, interval=interval,midpos=interval.start)
                result_left = all_peaks_in_interval.iloc[(len(all_peaks_in_interval)//2):, :]
                result_left = result_left.append(self.get_mock_data(N=(N - len(result_left)), interval=interval,
                                                                    midpos=result_left["mids"].iat[0]))
                assert len(result_right) == len(result_right)
                return [result_right, result_left]
            result_right = all_peaks_in_interval.iloc[0:(len(all_peaks_in_interval)//2), :]
            result_left = all_peaks_in_interval.iloc[(len(all_peaks_in_interval)//2):, :] # if len(all_peaks_in_interval)%2==1, peak always append to result_left
            if N_is_strict:
                result_right = result_right.append(self.get_mock_data(N=(N - len(result_right)), interval= interval,
                                            midpos=result_right["mids"].iat[-1])) #mids is the coordinates of the last peak in result_right
                result_left = result_left.append(self.get_mock_data(N=(N - len(result_left)), interval= interval,
                                         midpos=result_left["mids"].iat[0]))
                result_left.sort_values(['mids', 'sigVal'], inplace= True)
                assert len(result_right) == len(result_right)
                return[result_right, result_left]
        elif len(all_peaks_in_interval) == 0:
            if N_is_strict:
                result_right = self.get_mock_data(N=N, interval=interval,midpos=interval.start)
                result_left = self.get_mock_data(N=N, interval=interval, midpos=interval.end)
                return [result_right, result_left]


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

    def read_orient_file(self, orient_fname):  # store peaks with orientation from gimmeMotifs as sorted pandas dataframe
        logging.log(msg="Reading CTCF_orientation file " + orient_fname, level=logging.INFO)

        # set random temporary labels
        Nfields = len(open(orient_fname).readline().strip().split())
        names = list(map(str, list(range(Nfields))))
        data = pd.read_csv(orient_fname, sep="\t", header=None, names=names)  #TODO check: CTCF fname == CTCF orient fname.split('-')[0]
        # print(data)

        # subset and rename
        data = data.iloc[:, [0, 1, 2, 4, 5]]
        data.rename(columns={"0": "chr", "1": "start", "2": "end", "4": "score", "5": "orientation"},
                    inplace=True)

        # check duplicats
        duplicated = data.duplicated(subset=["chr", "start", "end"])
        if sum(duplicated) > 0:
            logging.warning(
                "Duplicates by genomic positions found in file " + orient_fname)  # TODO check why this happens
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
        return chr_data


    def set_sites_orientation(self, orient_fname): #Set orientation of sites based on gimmeMotifsData
                                    #It will fill plus_orient_data and minus_orient_data cols
                                    #And set orient_data_real to True
        try:
            self.chr_data
        except:
            logging.error("Please read data first")
            return None
        orient_chr_data = self.read_orient_file(orient_fname)
        result_intersection_data = intersect_intervals(self.chr_data, orient_chr_data)
        for chr in result_intersection_data.keys():
            #if chr != 'chr4':
             #   continue
            result_intersection_data[chr].sort_values(by=["orientation", "intersection", "score"], inplace=True)
            #duplicated = result_intersection_data[chr].duplicated(subset=["intersection", "orientation"])
            #print(duplicated)
            result_intersection_data[chr].drop_duplicates(subset=["intersection", "orientation"], keep="first", inplace=True)
            #print(result_intersection_data['chr4'])
            plus_orient_data = result_intersection_data[chr].query("orientation =='+'")
            plus_col_ind = self.chr_data[chr].columns.get_loc("plus_orientation")
            plus_row_list = list(plus_orient_data["intersection"])
            self.chr_data[chr].iloc[plus_row_list, plus_col_ind] = list(plus_orient_data["score"])

            minus_orient_data = result_intersection_data[chr].query("orientation =='-'")
            minus_col_ind = self.chr_data[chr].columns.get_loc("minus_orientation")
            minus_row_list = list(minus_orient_data["intersection"])
            self.chr_data[chr].iloc[minus_row_list, minus_col_ind] = list(minus_orient_data["score"])
        self.orient_data_real = True

    def export2bed_files_with_orientation(self, out_folder): #Export data in bed-graph format
        if not self.orient_data_real:
            logging.error('please set orientation first')
        orient_plus_chr_data = pd.DataFrame()
        orient_minus_chr_data = pd.DataFrame()
        for chr in self.chr_data:
            orient_plus_chr_data = orient_plus_chr_data.append(self.chr_data[chr][['chr', 'start', 'end', 'plus_orientation']])
            orient_minus_chr_data = orient_minus_chr_data.append(self.chr_data[chr][['chr', 'start', 'end', 'minus_orientation']])
        orient_plus_chr_data.to_csv(out_folder + "orient_plus.bedGraph", sep='\t', header=False, index=False)
        orient_minus_chr_data.to_csv(out_folder + "orient_minus.bedGraph", sep='\t', header=False, index=False)
    #keep only sites with orientation and replace orientation score by zero for the lowest score for this motif
    def keep_only_with_orient_data(self):
        if not self.orient_data_real:
               logging.error("please set_orientation first")
        for chr in self.chr_data:
            if chr != 'chr1':
                continue
            print(self.chr_data[chr])
            print(self.chr_data[chr].query("plus_orientation!='0'&"
                                     "minus_orientation!='0'"))
            print(self.chr_data[chr].query("plus_orientation!='0'&"
                                     "minus_orientation!='0'").index)
            print(self.chr_data[chr].loc[27550])
            #print(self.chr_data[chr])
            line_names = self.chr_data[chr].query("plus_orientation!='0'&"
                                     "minus_orientation!='0'").index
            print(line_names[0])
            for line_name in line_names:
                if self.chr_data[chr].loc[line_name, 'plus_orientation'] > self.chr_data[chr].loc[line_name, 'minus_orientation']:
                    self.chr_data[chr].loc[line_name, 'minus_orientation'] = 0
                elif  self.chr_data[chr].loc[line_name, 'plus_orientation'] < self.chr_data[chr].loc[line_name, 'minus_orientation']:
                    self.chr_data[chr].loc[line_name, 'plus_orientation'] = 0
                else:
                    logging.error('')
            #print(line_indexes)
            # for row in self.chr_data[chr].query("plus_orientation!='0'&"
            #                          "minus_orientation!='0'"):
            #     print('j')
            # self.chr_data[chr].query("plus_orientation!='0'|"
            #                          "minus_orientation!='0'",inplace=True)
            # self.only_orient_peaks = True


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

    def toXMLDict(self):
        res = {"ProteinName":self.proteinName,"fname":self.fname,
               "orintation_set":self.orient_data_real,
                "only_orient_peaks_kept":self.only_orient_peaks}
        return res

