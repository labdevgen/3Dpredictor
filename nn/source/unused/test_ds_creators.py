import torch
from torch.utils.data import Dataset, DataLoader
import numpy as np
import random
import sys
import os
sys.path.append(os.path.dirname(os.getcwd()))
from shared import get_bin_size
import pandas as pd
import matplotlib.pyplot as plt
import logging
import datetime

# Returns: predicted, real
class DatasetMaker(Dataset):
    def make_test(self, **kwargs):
        size = kwargs["size"]
        r = np.ones(shape=size)
        p = np.ones(shape=size)
        st, en = kwargs["st"], kwargs["en"]
        for i in range(st, en):
            for j in range(st, en):
                r[i, j] += kwargs["TAD"]

        p = np.reshape(p, (1,) + p.shape) # Add one "color" dimension
        r = np.reshape(r, (1,) + r.shape)

        return torch.from_numpy(p).float(),torch.from_numpy(r).float()

    def __init__(self, numSamples, **kwargs):
        self.reals = []
        self.predicted = []
        for i in range(numSamples):
            p,r = self.make_test(**kwargs)
            self.reals.append(r)
            self.predicted.append(p)


    def __len__(self):
        return len(self.reals)

    def __getitem__(self, idx):
        return self.predicted[idx], self.reals[idx]

class DatasetMaker_normalNoize(DatasetMaker):
    def make_test(self, **kwargs):
        size = kwargs["size"]
        r = np.abs(np.random.normal(0,1,size=size))
        p = np.abs(np.random.normal(0,1,size=size))
        st, en = kwargs["st"], kwargs["en"]
        for i in range(st, en):
            for j in range(st, en):
                r[i, j] += kwargs["TAD"]
        p = np.reshape(p, (1,) + p.shape) # Add one "color" dimension
        r = np.reshape(r, (1,) + r.shape)

        return torch.from_numpy(p).float(),torch.from_numpy(r).float()

class DatasetMaker_normalNoize_withLoop(DatasetMaker):
    def make_test(self, **kwargs):
        size = kwargs["size"]
        r = np.abs(np.random.normal(0,1,size=size))
        p = np.abs(np.random.normal(0,1,size=size))
        st, en = kwargs["st"], kwargs["en"]
        for i in range(st, en):
            for j in range(st, en):
                r[i, j] += kwargs["TAD"]
        # Add loop
        for i in list(range(st, en)):
            for j in list(range(st, en)):
                r[i,j] += kwargs["TAD"]*2
                r[j,i] += kwargs["TAD"]*2
        p = np.reshape(p, (1,) + p.shape) # Add one "color" dimension
        r = np.reshape(r, (1,) + r.shape)
        return torch.from_numpy(p).float(),torch.from_numpy(r).float()

class DatasetMaker_moving_TAD(DatasetMaker):
    def __init__(self, **kwargs):
        random.seed()
        super().__init__(**kwargs)

    def make_test(self, **kwargs):
        size = kwargs["size"]
        r = np.ones(shape=size)
        p = np.ones(shape=size)
        st, en, maxlen = kwargs["st"], kwargs["en"], kwargs["maxlen"]
        maxlen = min(maxlen,size[0],size[1])
        TAD_len = random.randint(4,maxlen)
        st = random.randint(0,en-TAD_len-1)
        en = st + TAD_len
        for i in list(range(st, en)):
            for j in list(range(st, en)):
                r[i, j] += kwargs["TAD"]
                p[i, j] += kwargs["TAD"]

        # Add loop
        for i in range(st,st+1):
            for j in range(en, en+1):
                r[i,j] += kwargs["TAD"]*2
                r[j,i] += kwargs["TAD"] * 2

        p = np.reshape(p, (1,) + p.shape) # Add one "color" dimension
        r = np.reshape(r, (1,) + r.shape)

        return torch.from_numpy(p).float(),torch.from_numpy(r).float()

class DatasetMaker_moving_TAD_with_Noize(DatasetMaker_moving_TAD):
    def make_test(self, **kwargs):
        size = kwargs["size"]
        r = np.abs(np.random.normal(0,1,size=size))
        p = np.abs(np.random.normal(0,1,size=size))
        st, en, maxlen = kwargs["st"], kwargs["en"], kwargs["maxlen"]
        maxlen = min(maxlen,size[0],size[1])
        TAD_len = random.randint(4,maxlen)
        st = random.randint(0,en-TAD_len-1)
        en = st + TAD_len
        for i in range(st, en):
            for j in range(st, en):
                r[i, j] += kwargs["TAD"]
                p[i, j] += kwargs["TAD"]

        # Add loop
        for i in range(st-1,st+2):
            for j in range(en-1, en+2):
                r[i,j] += kwargs["TAD"]*2
                r[j,i] += kwargs["TAD"] * 2
        p = np.reshape(p, (1,) + p.shape) # Add one "color" dimension
        r = np.reshape(r, (1,) + r.shape)

        return torch.from_numpy(p).float(),torch.from_numpy(r).float()


class DatasetMaker_moving_TAD_with_flying_Loop(DatasetMaker):
    def __init__(self, **kwargs):
        random.seed()
        super().__init__(**kwargs)

    def make_test(self, **kwargs):
        size = kwargs["size"]
        r = np.ones(shape=size)
        p = np.ones(shape=size)
        st, en, maxlen = kwargs["st"], kwargs["en"], kwargs["maxlen"]
        maxlen = min(maxlen,size[0],size[1])
        TAD_len = random.randint(6,maxlen)
        st = random.randint(st,en-TAD_len-1)
        en = st + TAD_len
        for i in list(range(st, en)):
            for j in list(range(st, en)):
                r[i, j] += kwargs["TAD"]
                p[i, j] += kwargs["TAD"]

        # Add loop
        fly_distance = 2
        r[st + TAD_len // 3, en + 2] = kwargs["TAD"]*3
        #r[en - 2, en + 2] = kwargs["TAD"] * 3
        #for i in range(st+fly_distance,st+fly_distance+1):
        #    for j in range(en + fly_distance, en + fly_distance + 1):
        #        r[i,j] += kwargs["TAD"]*3
        #        r[j,i] += kwargs["TAD"]*3

        coeff = np.sum(p)
        #p = p / coeff
        #r = r / coeff

        p = np.reshape(p, (1,) + p.shape) # Add one "color" dimension
        r = np.reshape(r, (1,) + r.shape)

        return torch.from_numpy(p).float(),torch.from_numpy(r).float()

class DatasetMaker_moving_TAD_with_flying_Loop_and_Noize(DatasetMaker):
    def __init__(self, **kwargs):
        random.seed()
        super().__init__(**kwargs)

    def make_test(self, **kwargs):
        size = kwargs["size"]
        r = np.random.normal(0,1,size=size)
        p = np.random.normal(0,1,size=size)
        st, en, maxlen = kwargs["st"], kwargs["en"], kwargs["maxlen"]
        maxlen = min(maxlen,size[0],size[1])
        TAD_len = random.randint(3,maxlen)
        st = random.randint(st,en-TAD_len-1)
        en = st + TAD_len
        for i in list(range(st, en)):
            for j in list(range(st, en)):
                r[i, j] += kwargs["TAD"]
                p[i, j] += kwargs["TAD"]

        # Add loop
        fly_distance = 2
        for i in range(st+fly_distance,st+fly_distance+1):
            for j in range(en + fly_distance, en + fly_distance + 1):
                r[i,j] += kwargs["TAD"]*3
                r[j,i] += kwargs["TAD"]*3

        coeff = np.sum(p)
        #p = p / coeff
        #r = r / coeff

        p = np.reshape(p, (1,) + p.shape) # Add one "color" dimension
        r = np.reshape(r, (1,) + r.shape)

        return torch.from_numpy(p).float(),torch.from_numpy(r).float()

class DatasetFromRealAndPredicted():
    def split2pairs(self, matrix, length, step):
        # split large matrix to regions which would be feeded to nn
        # length - length (along diagonal) of each region
        # step - distance (along diagonal) between start of each region

        # TODO this is now not really memory efficient, since overlapping part of array are stored as independent elements
        # if this will be the point, one may move core of this function to the .__get__ func, to get views from matrix
        # on the fly as they needed

        assert matrix.shape[0] == matrix.shape[1] # matrix supposed to be squared
        self.reals = []
        self.predicted = []
        I,J = np.triu_indices(length)

        for i in range(0, len(matrix)-length, step):
            real = np.array(matrix[i:i+length,i:i+length])
            predicted = np.array(matrix[i:i+length,i:i+length])

            real[J,I] = real[I,J]
            predicted[I,J] = predicted[J,I]

            # Add a 'color' dimension
            real = np.reshape(real,((1,) + real.shape))
            predicted = np.reshape(predicted,((1,) + predicted.shape))
            self.reals.append(torch.from_numpy(real).float())
            self.predicted.append(torch.from_numpy(predicted).float())

    def __init__(self, filename):
        # def add_contact(series):
        #     i = int((series["contact_st"] - data_min) // binsize)
        #     j = int((series["contact_en"]- data_min) // binsize)
        #     predicted = series["contact_count"]
        #     real = series["0"]
        #     matrix[i,j] = predicted
        #     matrix[j,i] = real

        logging.basicConfig(format='%(asctime)s %(name)s: %(message)s',
                                                datefmt='%I:%M:%S',
                                                level=logging.DEBUG)

        assert np.iinfo(np.dtype("uint32")).max > 250000000
        data = pd.read_csv(filename,sep=" ",dtype={"contact_st" : np.uint32,
                                                   "contact_en" : np.uint32,
                                                   "contact_count": np.float32,
                                                   "0": np.float32,
                                                   "IsLoop": np.uint8})
        logging.getLogger(__name__).debug("Reading data")
        data_min, data_max = data["contact_st"].min(), data["contact_en"].max()
        binsize = int(get_bin_size(data))
        logging.getLogger(__name__).debug("Using bin size "+str(binsize))
        assert data_min % binsize == data_max % binsize == 0
        assert data_max > data_min
        chr_len = int((data_max - data_min) // binsize)
        logging.getLogger(__name__).debug("Data size: " + str(chr_len) + " bins")
        matrix = np.zeros(shape=(chr_len + 1, chr_len + 1))


        # Fill matrix
        logging.getLogger(__name__).debug("Filling matrix")

        data["contact_st"] = ((data["contact_st"] - data_min) // binsize).astype(int)
        data["contact_en"] = ((data["contact_en"] - data_min) // binsize).astype(int)
        i = data["contact_st"].values
        j = data["contact_en"].values
        matrix[(i,j)] = data["contact_count"].values
        matrix[(j,i)] = data["0"].values

        #data.apply(add_contact,axis = "columns")

        diag_sums = [np.trace(matrix,i) for i in range(len(matrix))]

        assert diag_sums[0] + diag_sums[1] == 0 # check that first 2 diagonals are empty
        assert np.trace(matrix, 1500000 // binsize - 1) != 0 # check that we have disctances up to 1.5 Mb

        #TODO apply transformations to data here
        logging.getLogger(__name__).debug("Normalizing data")

        mean = -100 # find mean along first non-zero diagonal
        for ind,val in enumerate(diag_sums):
            if val != 0:
                mean = val / (len(matrix) - ind)

        if mean - 1 <= 1: # For oe values mean is ~1
            logging.getLogger(__name__).info("Assuming contacts, going to convert to o/e values.")

            for i in range(len(matrix)):
                if diag_sums[i] == 0:
                    continue
                else:
                    for j in range(i,len(matrix)):
                        matrix[j-i, j] = matrix[j-i,j] / diag_sums[i]
                        matrix[j,j-i] = matrix[j,j-i] / diag_sums[i]

        self.split2pairs(matrix,length=1500000 // binsize, step=750000 // binsize)

    def __len__(self):
        return len(self.reals)

    def __getitem__(self, idx):
        return self.predicted[idx], self.reals[idx]
