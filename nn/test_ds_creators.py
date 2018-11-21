import torch
from torch.utils.data import Dataset, DataLoader
import numpy as np
import random

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
        TAD_len = random.randint(4,maxlen)
        st = random.randint(st,en-TAD_len-1)
        en = st + TAD_len
        for i in list(range(st, en)):
            for j in list(range(st, en)):
                r[i, j] += kwargs["TAD"]
                p[i, j] += kwargs["TAD"]

        # Add loop
        fly_distance = 0
        for i in range(st+fly_distance,st+fly_distance+1):
            for j in range(en + fly_distance, en + fly_distance + 1):
                r[i,j] += kwargs["TAD"]*3
                #r[j,i] += kwargs["TAD"]*3

        coeff = np.sum(p)
        #p = p / coeff
        #r = r / coeff

        p = np.reshape(p, (1,) + p.shape) # Add one "color" dimension
        r = np.reshape(r, (1,) + r.shape)

        return torch.from_numpy(p).float(),torch.from_numpy(r).float()
