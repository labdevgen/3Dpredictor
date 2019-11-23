# Written by Minja and Polina, 08-09 2018
# A general module purpose is generate predictors
# Class DataGenerator does it by applying contact2file to a dataframe
# which contains contacts information in a format contact_st -- contact_en -- contact_count
# To generate predictors for each contact, we use instances of PredictorGenerators
# Each instance should be able to generate value(s) of specific predicor, e.g. ChiPSeq signal,
# RNAseq of 1st Eigenvecor
# Function contact2file aggregates these predictors and writes data to file

import logging
import numbers
import numpy as np
from shared import Interval
from collections import  OrderedDict

class PredictorGenerator(object):
    def __init__(self,**kwargs):
        for key,val in kwargs.items():
            self.__setattr__(key,val)
        self.vectorizable = False

    def get_header(self,contact): #Returns predictor header line, i.e. name(s) of predictors
        pass

    def get_predictors(self,contact): #predictor value(s)
        pass

    def symmetric_window_around_contact(self,contact):
        #Returns interval with a size self.window_size center around region spanning contact
        window_start = max(0, contact.contact_st - ((self.window_size - contact.dist) // 2))
        window_end = window_start + self.window_size
        contacts_relative_start = contact.contact_st - window_start
        assert contacts_relative_start >= 0
        contacts_relative_end = contact.contact_en - window_start
        assert contacts_relative_end >= contacts_relative_start
        return window_start,window_end,contacts_relative_start,contacts_relative_end

    def intevals_around_ancor(self,contact):
        half = self.window_size // 2
        assert contact.contact_en - contact.contact_st > half
        return (Interval(contact.chr,
                         contact.contact_st - half,
                         contact.contact_st + half),
                Interval(contact.chr,
                         contact.contact_st + half,
                         contact.contact_en - half),
                Interval(contact.chr,
                         contact.contact_en - half,
                         contact.contact_en + half ),
                )
    def print_warnings_occured_during_predGeneration(self):
        pass

    def toXMLDict(self,contact):
        res = OrderedDict()
        res["name"] = self.name
        res["predictors"] = ",".join(self.get_header(contact))

        # Here we want to run toXMLDict(contact) from each reader
        # but we don't know the reader variable name
        # To find reader we iter over all class fields
        # Untill we get a member with .toXMLDict function avaliable
        for k,v in self.__dict__.items():
            try:
                res[k] = self.__dict__[v].toXMLDict
            except:
                pass
            if (isinstance(v, numbers.Number)): #Also save all numeric parameters
                res[k] = v
        return res


class ChipSeqPredictorGenerator(PredictorGenerator):
    def __init__(self, chipSeq_reader, binsize, window_size, **kwargs):
        self.name =  chipSeq_reader.proteinName
        self.chipSeq_reader = chipSeq_reader
        self.binsize = binsize
        self.window_size = window_size
        super(ChipSeqPredictorGenerator, self).__init__(**kwargs)

    def get_header(self,contact):
        self.header = [self.name + "win_st", self.name + "win_en",
                       self.name + "win_rel_st", self.name + "win_rel_en"] \
                      + [self.name + "_bin" + str(i) for i in range(self.window_size // self.binsize)]
        return self.header

    def get_predictors(self,contact):
        window_start, window_end, contacts_relative_start, contacts_relative_end = \
                                        self.symmetric_window_around_contact(contact)
        interval = Interval(contact.chr, window_start, window_end)
        return [window_start,window_end,contacts_relative_start,contacts_relative_end] \
               +self.chipSeq_reader.get_binned_interval(interval, binsize=self.binsize)

    def delete_region(self,interval):
        self.chipSeq_reader.delete_region(interval)

class SmallChipSeqPredictorGenerator(ChipSeqPredictorGenerator):
    # Less predictors:
    # 1. Total CTCF in 3 regions: L_ancor + winsize/2 ... R_ancor - winsize/2
    #                      -winsize/2 Lancor +winsize/2
    #                      -winsize/2 Rancor +winsize/2
    # 2. Distances and sigVral for K nearest CTCF peaks from L_ancor and R_ancor
    def __init__(self, chipSeq_reader, window_size, N_closest, **kwargs):
        self.name = chipSeq_reader.proteinName + "_SmallChip"
        self.N_closest = N_closest
        #print(self.name)
        super(SmallChipSeqPredictorGenerator, self).__init__(chipSeq_reader, 0, window_size, **kwargs)

    def get_header(self,contact):
        if self.name == "CTCF_SmallChip":
            self.header = [self.name + "_L", self.name + "_R"]
        else:
            self.header = [self.name + "_L", self.name + "_W", self.name + "_R"]
            for side in ["L", "R"]:
                for metric in ["Sig", "Dist"]:
                    for i in range(self.N_closest):
                        self.header += [self.name + "_" + side + metric + "_" + str(i)]
        return self.header

    def get_predictors(self,contact):
        #print(self.name)
        assert contact.contact_st < contact.contact_en
        if contact.chr not in set(self.chipSeq_reader.chr_data.keys()):
            return [0]*len(self.header)
        else:
            intL,intM,intR = self.intevals_around_ancor(contact)
            sig_L = self.chipSeq_reader.get_interval(intL).sigVal.sum()
            sig_R = self.chipSeq_reader.get_interval(intR).sigVal.sum()
            sig_mid = self.chipSeq_reader.get_interval(intM).sigVal.sum()
            Left_top = self.chipSeq_reader.get_nearest_peaks(Interval(contact.chr,
                                                                      contact.contact_st - (self.window_size // 2),
                                                                      contact.contact_st - (self.window_size // 2)),
                                                             N=self.N_closest, side="left")

            Left_top = Left_top["sigVal"].values.tolist() + \
                       (contact.contact_st-Left_top["mids"]).values.tolist()

            Right_top = self.chipSeq_reader.get_nearest_peaks(Interval(contact.chr,
                                                                       contact.contact_en + (self.window_size // 2),
                                                                       contact.contact_en + (self.window_size // 2)),
                                                              N=self.N_closest, side="right")
            Right_top = Right_top["sigVal"].values.tolist() + \
                        (Right_top["mids"]-contact.contact_en).values.tolist()
            if self.name == "CTCF_SmallChip":
                return [sig_L,sig_R]
            else:
                return [sig_L,sig_mid,sig_R]+Left_top+Right_top


class E1PredictorGenerator(PredictorGenerator):
    def __init__(self,eig_reader,window_size,name="E1"):
        self.eig_reader = eig_reader
        self.window_size = window_size
        self.name = name
        self.vectorizable = False

    def get_header(self,contact):
        self.header = [self.name + "win_st",self.name + "win_en",self.name + "win_rel_st",self.name + "win_rel_en"]
        self.header += [self.name + "_bin" + str(i) \
                         for i in range(self.window_size // self.eig_reader.get_binsize())]
        return self.header

    def get_predictors(self,contact):
        assert contact.contact_st < contact.contact_en
        window_start, window_end, contacts_relative_start, contacts_relative_end = \
                                        self.symmetric_window_around_contact(contact)
        interval = Interval(contact.chr, window_start, window_end)
        return [window_start, window_end, contacts_relative_start, contacts_relative_end] + \
                self.eig_reader.get_E1inInterval(interval)["E1"].tolist()

    def print_warnings_occured_during_predGeneration(self):
        self.eig_reader.print_warnings()

    def delete_region(self,interval):
        self.eig_reader.delete_region(interval)

class SmallE1PredictorGenerator(E1PredictorGenerator):
    # Less predictors:
    # 1. E1 in +/- winsize region around contact ancors
    def get_header(self,*args):
        self.header = [self.name + "_L",self.name + "_R",
                       self.name + "_N_Changes", self.name + "_SD"]
        return self.header

    def get_predictors(self,contact):
        intL,intM,intR = self.intevals_around_ancor(contact)
        E1_M = self.eig_reader.get_E1inInterval(intM)["E1"].values
        if len(E1_M) == 0:
            E1_SD = 0
            E1_N_Changes = 0
        else:
            E1_SD = np.std(E1_M) #STD error of E1 values
            E1_N_Changes = len(E1_M) - 1 - sum(np.equal(np.sign(E1_M[1:]),np.sign(E1_M[:-1]))) #Number of changes of sign
        return list(map(np.average,[self.eig_reader.get_E1inInterval(intL)["E1"].tolist(),
               self.eig_reader.get_E1inInterval(intR)["E1"].tolist()]))+[E1_N_Changes,E1_SD]

class SitesOrientPredictorGenerator(PredictorGenerator):
    def __init__(self, chipSeq_reader, N_closest, **kwargs):
        self.chipSeq_reader = chipSeq_reader
        if self.chipSeq_reader.only_orient_peaks:
            self.name = chipSeq_reader.proteinName + '_OnlySitesOrient'
        else:
            self.name = chipSeq_reader.proteinName + '_SitesOrient'
        self.N_closest = N_closest
        if not self.chipSeq_reader.orient_data_real:
            logging.error('please set orientation first')
            raise Exception("Can't generate predictions")
        self.vectorizable = False

    def get_header(self,contact):
        self.header = []
        for side in "L", "W_L", "W_R", "R":
            for metric in ["plus_orientation", "minus_orientation","sigVal", "dist"]:
                for i in range(self.N_closest):
                    self.header += [self.name + "_" + side + "_" + metric + "_" + str(i)]
        self.header += [self.name + "_W_sumSigVal"]
        self.header += [self.name + "_W_N+orient", self.name + "_W_N-orient"]

        #print('=================================')
        #print(len(self.header))
        return self.header

    def get_predictors(self,contact):
        #print(self.name)
        assert contact.contact_st < contact.contact_en

        # Peaks outside of the window
        Left_peaks = self.chipSeq_reader.get_nearest_peaks(Interval(contact.chr, contact.contact_st, contact.contact_st ),
                                                         N=self.N_closest, side="left")

        Left_peaks = Left_peaks["plus_orientation"].values.tolist() + \
                   Left_peaks["minus_orientation"].values.tolist() + Left_peaks["sigVal"].values.tolist() + \
                     (contact.contact_st - Left_peaks["mids"]).values.tolist()

        Right_peaks = self.chipSeq_reader.get_nearest_peaks(Interval(contact.chr, contact.contact_en, contact.contact_en),
                                                            N=self.N_closest, side="right")

        Right_peaks = Right_peaks["plus_orientation"].values.tolist() + \
                     Right_peaks["minus_orientation"].values.tolist() + Right_peaks["sigVal"].values.tolist() + \
                     (Right_peaks["mids"] - contact.contact_en).values.tolist()

        #Next statmetn will return list of 2 dataframes
        #1st DF with first N peaks on the right side of left interval boundary
        #2nd DF with first N peaks on the left side of right interval boundary
        Window_peaks = self.chipSeq_reader.get_N_peaks_near_interval_boundaries(Interval(contact.chr, contact.contact_st, contact.contact_en),
                                                                           N=self.N_closest)

        Window_peaks_left = Window_peaks[0]["plus_orientation"].values.tolist() +Window_peaks[0]["minus_orientation"].values.tolist() + \
                             Window_peaks[0]["sigVal"].values.tolist() + (Window_peaks[0]["mids"] - contact.contact_st).values.tolist()

        #Get properties of peaks
        Window_peaks_right = Window_peaks[1]["plus_orientation"].values.tolist() + Window_peaks[1]["minus_orientation"].values.tolist() + \
            Window_peaks[1]["sigVal"].values.tolist() + (contact.contact_en - Window_peaks[1]["mids"]).values.tolist()

        #if there are no peaks in window, set sigVal and other params to 0 TODO add if/else for onlyOrient
        if len(self.chipSeq_reader.get_interval(Interval(contact.chr, contact.contact_st, contact.contact_en))) == 0:
            Window_sigVal = 0
            N_plus_orient_W = 0
            N_minus_orient_W = 0
        else:
            Window_sigVal = self.chipSeq_reader.get_interval(Interval(contact.chr, contact.contact_st, contact.contact_en)).sigVal.sum()
            plus_orient_data = self.chipSeq_reader.get_interval(Interval(contact.chr, contact.contact_st, contact.contact_en)).query("plus_orientation!='0'")
            N_plus_orient_W = len(plus_orient_data)
            minus_orient_data = self.chipSeq_reader.get_interval(Interval(contact.chr, contact.contact_st, contact.contact_en)).query("minus_orientation!='0'")
            N_minus_orient_W = len(minus_orient_data)

        predictors = Left_peaks + Window_peaks_left + Window_peaks_right + Right_peaks + [Window_sigVal] + [N_plus_orient_W, N_minus_orient_W]
        return predictors

class OrientBlocksPredictorGenerator(PredictorGenerator): #this PG
    def __init__(self, chipSeq_reader, window_size, **kwargs):
            self.name = chipSeq_reader.proteinName + '_OrientBlock'
            self.chipSeq_reader = chipSeq_reader
            self.window_size = window_size
            if not self.chipSeq_reader.orient_data_real:
                logging.error('please set orientation first')
            if not self.chipSeq_reader.only_orient_peaks:
                logging.error('please get data with orientations only first')
            self.vectorizable = False
    def get_header(self,contact):
        self.header = [self.name + "_W_NBlocks", self.name + "_HasDivergOrient"]
        return self.header
    def get_predictors(self,contact):
        assert contact.contact_st < contact.contact_en

        # Get peaks in window and count "blocks"
        # Blocks are CTCF sites with divergent orientation, i.e. --> <-- <-- is a block
        all_Window_peaks = self.chipSeq_reader.get_interval(Interval(contact.chr, contact.contact_st + self.window_size//2, \
                                                                     contact.contact_en - self.window_size//2))
        N_blocks_W = 0
        plus_ori_idx = all_Window_peaks.columns.get_loc('plus_orientation')
        minus_ori_idx = all_Window_peaks.columns.get_loc('minus_orientation')
        # print(len(all_Window_peaks))
        # print(all_Window_peaks[["minus_orientation", "plus_orientation"]])
        for i in range(len(all_Window_peaks) - 1):
            if all_Window_peaks.iloc[i,plus_ori_idx ]!=0 and all_Window_peaks.iloc[i+1,minus_ori_idx] !=0:
                N_blocks_W +=1
        # print(N_blocks_W)

        # Check wheather we have CTCF sites in divergent orientation in contact ancors
        intL, intM, intR = self.intevals_around_ancor(contact)
        L_peaks = self.chipSeq_reader.get_interval(intL)
        R_peaks = self.chipSeq_reader.get_interval(intR)
        has_convergent_peak = 0
        if len(L_peaks) > 0 and len(R_peaks) > 0:
            has_convergent_peak = L_peaks.plus_orientation.sum()*R_peaks.minus_orientation.sum()
        return [N_blocks_W,has_convergent_peak]

class ConvergentPairPredictorGenerator(PredictorGenerator):
    def __init__(self, chipSeq_reader,binsize, **kwargs):
            self.name = chipSeq_reader.proteinName + '_ConvergentPair'
            self.chipSeq_reader = chipSeq_reader
            self.binsize=binsize
            if not self.chipSeq_reader.orient_data_real:
                logging.error('please set orientation first')
            if not self.chipSeq_reader.only_orient_peaks:
                logging.error('please get data with orientations only first')
            self.vectorizable = False
    def get_header(self,contact):
        self.header = [self.name + "_HasConvergentPair"]
        return self.header
    def get_predictors(self,contact):
        assert contact.contact_st < contact.contact_en
        # get the nearest right and left peak to the start and to the end of contact
        Left_start_peaks = self.chipSeq_reader.get_nearest_peaks(
            Interval(contact.chr, contact.contact_st, contact.contact_st),N=1, side="left")
        Left_start_peak = Left_start_peaks["minus_orientation"].values.tolist()[0] \
            if abs(Left_start_peaks["mids"].values.tolist()[0]-contact.contact_st) <= self.binsize*2 else 0
        Right_start_peaks = self.chipSeq_reader.get_nearest_peaks(
            Interval(contact.chr, contact.contact_st, contact.contact_st), N=1, side="right")
        Right_start_peak = Right_start_peaks["minus_orientation"].values.tolist()[0] \
            if abs(Right_start_peaks["mids"].values.tolist()[0]-contact.contact_st) <= self.binsize*2 else 0
        Left_end_peaks = self.chipSeq_reader.get_nearest_peaks(
            Interval(contact.chr, contact.contact_en, contact.contact_en), side="left", N=1)
        Left_end_peak = Left_end_peaks["minus_orientation"].values.tolist()[0] \
            if abs(Left_end_peaks["mids"].values.tolist()[0] - contact.contact_en) <= self.binsize*2 else 0
        Right_end_peaks = self.chipSeq_reader.get_nearest_peaks(
            Interval(contact.chr, contact.contact_en, contact.contact_en), side="right", N=1)
        Right_end_peak = Right_end_peaks["minus_orientation"].values.tolist()[0] \
            if abs(Right_end_peaks["mids"].values.tolist()[0] - contact.contact_en) <= self.binsize*2 else 0
        #minus orientation is orientation of CTCF to the right, plus to the left
        # 1 if CTCF sites in the end and start of contact have convergent orientation else 0
        start_minus_orientation = [Left_start_peak,Right_start_peak]
        end_plus_orientation = [Left_end_peak,Right_end_peak]
        if len(np.nonzero(start_minus_orientation)[0]) != 0 and len(np.nonzero(end_plus_orientation)[0]) != 0:
            predictors = [1]
        else:
            predictors = [0]
        return predictors


class SitesOnlyOrientPredictorGenerator(PredictorGenerator):
    def __init__(self, chipSeq_reader, N_closest, **kwargs):
            self.name = chipSeq_reader.proteinName + '_OnlyOrient'
            self.chipSeq_reader = chipSeq_reader
            self.N_closest = N_closest
            if not self.chipSeq_reader.orient_data_real:
                logging.error('please set orientation first')
            if not self.chipSeq_reader.only_orient_peaks:
                logging.error('please get data with orientations only first')
    def get_header(self,contact):
        self.header = []
        for contact_point in"start", "end":
            for side in "L", "R":
                for metric in ["+_orient", "-_orient","sigVal", "dist"]:
                    for i in range(self.N_closest):
                        self.header += [self.name + "_" + contact_point + "_" + side + "_" + metric + "_" + str(i)]
        # print("header", self.header)
        return self.header
    def get_predictors(self,contact):
        assert contact.contact_st < contact.contact_en
        #Peaks outside the window
        Left_start_peaks = self.chipSeq_reader.get_nearest_peaks(
            Interval(contact.chr, contact.contact_st, contact.contact_st),
            N=self.N_closest, side="left")
        Left_start_peaks = Left_start_peaks["plus_orientation"].values.tolist() + \
                           Left_start_peaks["minus_orientation"].values.tolist() + Left_start_peaks["sigVal"].values.tolist() + \
                     (contact.contact_st - Left_start_peaks["mids"]).values.tolist()
        Right_end_peaks = self.chipSeq_reader.get_nearest_peaks(
            Interval(contact.chr, contact.contact_en, contact.contact_en), side="right", N=self.N_closest)

        Right_end_peaks = Right_end_peaks["plus_orientation"].values.tolist() + \
                          Right_end_peaks["minus_orientation"].values.tolist() + Right_end_peaks[
                              "sigVal"].values.tolist() + \
                          (Right_end_peaks["mids"] - contact.contact_en).values.tolist()
        # Next statmetn will return list of 2 dataframes
        # 1st DF with first N peaks on the right side of left interval boundary
        # 2nd DF with first N peaks on the left side of right interval boundary
        Window_peaks = self.chipSeq_reader.get_N_peaks_near_interval_boundaries(
            Interval(contact.chr, contact.contact_st, contact.contact_en),
            N=self.N_closest)

        Right_start_peaks = Window_peaks[0]["plus_orientation"].values.tolist() + Window_peaks[0][
            "minus_orientation"].values.tolist() + \
                            Window_peaks[0]["sigVal"].values.tolist() + (
                            Window_peaks[0]["mids"] - contact.contact_st).values.tolist()

        Left_end_peaks = Window_peaks[1]["plus_orientation"].values.tolist() + Window_peaks[1][
            "minus_orientation"].values.tolist() + \
                             Window_peaks[1]["sigVal"].values.tolist() + (
                             contact.contact_en - Window_peaks[1]["mids"]).values.tolist()
        predictors = Left_start_peaks + Right_start_peaks + Left_end_peaks + Right_end_peaks
        return predictors

class Distance_to_TSS_PG(PredictorGenerator):
    def __init__(self, TSS_reader, **kwargs):
        self.name = '_Distance_toTSS'
        self.TSS_reader = TSS_reader
        self.vectorizable = False
    def get_header(self,contact):
        self.header=[]
        for contact_point in"start", "end":
            for side in "L", "R":
                self.header +=[self.name + "_"+contact_point+"_"+side]
        return self.header
    def get_predictors(self,contact):
        assert contact.contact_st < contact.contact_en
        # get the nearest right and left peak to the start and to the end of contact
        Left_start_peaks = self.TSS_reader.get_nearest_peaks(
            Interval(contact.chr, contact.contact_st, contact.contact_st),N=1, side="left")
        Left_start_peak_dist = (contact.contact_st - Left_start_peaks["TSS"]).values.tolist()
        Right_start_peaks = self.TSS_reader.get_nearest_peaks(
            Interval(contact.chr, contact.contact_st, contact.contact_st), N=1, side="right")
        Right_start_peak_dist = (Right_start_peaks["TSS"] - contact.contact_st).values.tolist()
        Left_end_peaks = self.TSS_reader.get_nearest_peaks(
            Interval(contact.chr, contact.contact_en, contact.contact_en), side="left", N=1)
        Left_end_peak_dist = (contact.contact_en - Left_end_peaks["TSS"]).values.tolist()
        Right_end_peaks = self.TSS_reader.get_nearest_peaks(
            Interval(contact.chr, contact.contact_en, contact.contact_en), side="right", N=1)
        Right_end_peak_dist = (Right_end_peaks["TSS"] - contact.contact_en).values.tolist()
        return Left_start_peak_dist+Right_start_peak_dist+Left_end_peak_dist+Right_end_peak_dist