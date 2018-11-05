import sys, os
import math
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import sys
sys.path.append(os.path.dirname(os.getcwd()))
from matrix_plotter import MatrixPlotter
from shared import Interval

# Here comes torch
import torch.nn as nn
import torch.nn.functional as F



# Generate test sample
def generate_sample(length = 1500, binsize = 25, modifier = lambda x: x + 0.001):
    length = int(length / binsize)

    # Square matrix shape length * length, conts = exp(distance)
    counts = np.ones(shape=(length,length))
    for i in range(length):
        idxs = np.arange(length - i)
        counts[idxs,idxs + i] = (1 / math.exp(i))*10000000

    # Add some noise
    real = np.abs(np.random.normal(counts, counts / 2)) + np.min(counts)

    # Calculate "predicted"
    predicted = modifier(counts)

    #convert to pd dataframe
    for counts in [real, predicted]:
        data = np.array(list(np.ndenumerate(counts)))
        data = pd.DataFrame(data, columns= ["coord","contact_count"])
        data[["contact_st","contact_en"]] = data["coord"].apply(pd.Series)
        data.drop(["coord"],axis="columns",inplace=True)
        data["contact_st"] = (data["contact_st"] + 1) * binsize
        data["contact_en"] = (data["contact_en"] + 1) * binsize
        data = data.query("contact_st < contact_en")
        data["chr"] = ["chrT"]*len(data)
        yield data

def sample2vector(data):
    #Define binsize
    binsize = np.sort((np.concatenate((data["contact_en"].values,data["contact_st"].values))))
    binsize = binsize[1:] - binsize[:-1]
    binsize = min(binsize)
    assert binsize > 0

    l = max(data["contact_en"]) - min(data["contact_st"])
    assert l % binsize == 0
    l = l / binsize

    i = torch.LongTensor([data["contact_st"] / binsize, data["contact_en"] / binsize])
    v = torch.FloatTensor(data["contact_count"])
    dense = torch.sparse.FloatTensor(i, v, torch.Size([l,l])).to_dense()
    print (dense)


def generate_batch(size = 50):


# Helper fn to draw maps
def draw_matrix(real,predicted):
    mp = MatrixPlotter()
    mp.set_data(real)
    mp.set_control(predicted)
    matrix = np.log(matrix)
    tick_pos, tick_labels = mp.get_bins_strart_labels(maxTicksNumber=15)
    plt.xticks(tick_pos, tick_labels, rotation=45)
    plt.imshow(matrix, cmap="OrRd")
    plt.show()


class SimplestNet(nn.Module):
    # Lets make a simple Net, 1
    def __init__(self):
        super(Net, self).__init__()
        self.conv1 = nn.Conv2d(3, 6, 5)
        self.pool = nn.MaxPool2d(2, 2)
        self.conv2 = nn.Conv2d(6, 16, 5)
        self.fc1 = nn.Linear(16 * 5 * 5, 120)
        self.fc2 = nn.Linear(120, 84)
        self.fc3 = nn.Linear(84, 10)

    def forward(self, x):
        x = self.pool(F.relu(self.conv1(x)))
        x = self.pool(F.relu(self.conv2(x)))
        x = x.view(-1, 16 * 5 * 5)
        x = F.relu(self.fc1(x))
        x = F.relu(self.fc2(x))
        x = self.fc3(x)
        return x


net = Net()

real,predicted = list(generate_sample())
draw_matrix(real,predicted)
#print (real.head())
#print (predicted.head())

