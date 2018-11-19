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
import torch
import torch.nn as nn
import torch.nn.functional as F
import torch.optim as optim
from torch.utils.data import Dataset, DataLoader



class sampleDataset():
    # Generate test sample
    def generate_sample(self,length=750, binsize=25, modifier=lambda x: x):
        length = int(length / binsize)

        # Square matrix shape length * length, conts = exp(distance)
        counts = np.ones(shape=(length, length))
        for i in range(length):
            idxs = np.arange(length - i)
            #counts[idxs, idxs + i] = (1 / math.exp(i)) * 10000000#np.abs(np.random.normal(0,1)) #(1 / math.exp(i)) * 10000000
            #counts[idxs, idxs + i] = 1.#np.abs(np.random.normal(0,0.5,size=(length - i,)))

            # Add some noise
        #real = counts + 0.0001 #+ np.abs(np.random.normal(counts, counts / 2)) + np.min(counts)
        real = np.array(counts)

        # Add 2 TADs to real
        for i in range(length//4):
            for j in range(length//4):
                real[i,j] *= 2


        for i in range(length//4, length//2):
            for j in range(length//4, length//2):
                real[i,j] *= 2


        # Calculate "predicted"
        #predicted = modifier(counts)
        predicted = np.array(counts)

        # Add one TAD to predicted
        for i in range(length//4,length//2):
            for j in range(length//4,length//2):
                predicted[i,j] *= 2


        # normalize values to be in 0..1 range
        m = max(np.max(real),np.max(predicted))
        print (m)

        predicted /= m
        real /= m
        #plt.imshow(real)
        #plt.show()
        assert np.all(predicted <= 1)
        assert np.all(real <= 1)
        # convert to pd dataframe
        for counts in [real, predicted]:
            data = np.array(list(np.ndenumerate(counts)))
            data = pd.DataFrame(data, columns=["coord", "contact_count"])
            data[["contact_st", "contact_en"]] = data["coord"].apply(pd.Series)
            data.drop(["coord"], axis="columns", inplace=True)
            data["contact_st"] = (data["contact_st"] + 1) * binsize
            data["contact_en"] = (data["contact_en"] + 1) * binsize
            data = pd.DataFrame(data.query("contact_st < contact_en"))
            data["chr"] = ["chrT"] * len(data)
            yield data

    def sample2vector(self,data):  # returns Torch tensor of the data
        # Define binsize
        binsize = np.sort(np.unique(np.concatenate((data["contact_en"].values, data["contact_st"].values))))
        binsize = binsize[1:] - binsize[:-1]
        binsize = min(binsize)
        assert binsize > 0

        # switch to 0-based coordinates
        minpos = min(data["contact_en"].min(), data["contact_st"].min())
        assert minpos >= 0
        data["contact_en"] -= minpos
        data["contact_st"] -= minpos

        # Define length of region
        l = max(data["contact_en"]) - min(data["contact_st"])
        assert l % binsize == 0
        l = int(l / binsize) + 1

        # Convert coords to bins and create torch Tensor
        i = torch.LongTensor([(data["contact_st"] / binsize).astype(int).values,
                              (data["contact_en"] / binsize).astype(int).values])
        v = torch.FloatTensor(data["contact_count"].tolist())
        dense = torch.sparse.FloatTensor(i, v, torch.Size([l, l])).to_dense()

        for i in range(l):
            for j in range(i,l):
                dense[j,i] = dense[i,j]
        #plt.imshow(dense.numpy())
        #plt.show()
        return dense

    def __init__(self,size):
        self.data = []
        for i in range(size):
            print (i)
            real, predicted = list(map(self.sample2vector, self.generate_sample()))
            self.data.append({"real":real.unsqueeze_(0),
                              "predicted":predicted.unsqueeze_(0)})

    def __len__(self):
        return len(self.data)

    def __getitem__(self, item):
        return self.data[item]

# Helper fn to draw maps
def draw_matrix(real,predicted):
    mp = MatrixPlotter()
    mp.set_data(real)
    mp.set_control(predicted)
    matrix = mp.getMatrix4plot(Interval(real["chr"].iloc[0],
                                        min(real["contact_st"].values),
                                        max(real["contact_en"].values)))
    matrix = np.log(matrix)
    tick_pos, tick_labels = mp.get_bins_strart_labels(maxTicksNumber=15)
    plt.xticks(tick_pos, tick_labels, rotation=45)
    plt.imshow(matrix, cmap="OrRd")
    plt.show()

class SimplestNet(nn.Module):
        # Lets make a simple Net, 1
        def __init__(self):
            super(SimplestNet, self).__init__()
            self.conv1 = nn.Conv2d(1, 8, 7, padding=3)  # in_channels, out_channels, kernel_size
            self.conv2 = nn.Conv2d(8, 1, 1)

        def forward(self, x):
            print("start forwarding")
            print (x.shape)
            x = self.conv1(x)
            x = F.relu(x)
            print (x.shape)
            x = self.conv2(x)
            x = F.relu(x)
            #print (x.shape)
            #x = self.conv3(x)
            #x = F.relu(x)
            print (x.shape)
            return x


if __name__ == '__main__':


    net = SimplestNet()
    print(net)
    criterion = nn.MSELoss()
    optimizer = optim.SGD(net.parameters(), lr=0.001)
    #print ("Generating data")
    testset = sampleDataset(size=4)
    #print ("Seting TrainLoader")
    trainloader = DataLoader(testset, batch_size=4,shuffle=False, num_workers=1)

    learning_rate = 0.1
    net.zero_grad()
    r,p = testset.generate_sample()
    r = testset.sample2vector(r)
    #r.requires_grad_(True)
    r = r.unsqueeze(0).unsqueeze(0)

    p = testset.sample2vector(p)
    #p.requires_grad_(True)
    p = p.unsqueeze(0).unsqueeze(0)

    print (p.shape)
    print (r.shape)

    #p = torch.randn(10,10,dtype=torch.float)
    #p = p.unsqueeze(0).unsqueeze(0)

    #r = torch.randn(10,10,dtype=torch.float)
    #r = r.unsqueeze(0).unsqueeze(0)

    print ("p is : ",p)
    print (p.shape)
    print ("r is : ",r)
    print (r.shape)
    for i in range(5000):
        out = net.forward(p)
        loss = criterion(out, r)
        print (i,"loss=",loss.item())
        loss.backward()
        all_grads = []
        if True:#torch.no_grad():
           for y in net.parameters():
                #y.data.sub_(y.grad.data * learning_rate)
                all_grads.append(y.grad.numpy().flatten())
           all_zero = np.all(np.concatenate(tuple(all_grads))==0)
           optimizer.step()
           net.zero_grad()
           print ("-----------All gradients are zero: ",all_zero)
           if all_zero:
                break
    plt.subplot(131)
    vmin, vmax = np.min(r.detach().numpy()), np.max(r.detach().numpy())
    plt.imshow(r.detach().numpy()[0,0,:,:], vmin=vmin, vmax=vmax)
    plt.subplot(132)
    plt.imshow(p.detach().numpy()[0,0,:,:], vmin=vmin, vmax=vmax)
    plt.subplot(133)
    plt.imshow(net.forward(p).detach().numpy()[0,0,:,:], vmin=vmin, vmax=vmax)
    print(p)
    print(r)
    print(net.forward(p))
    plt.show()
    print("-----------All gradients are zero: ", all_zero)
    sys.exit(0)



    print ("Training")
    for epoch in range(10):  # loop over the dataset multiple times
        running_loss = 0.0
        for i, data in enumerate(trainloader, 0):
            # get the inputs
            inputs, labels = data["real"],data["predicted"]

            # zero the parameter gradients
            optimizer.zero_grad()
            net.zero_grad()

            # forward + backward + optimize
            outputs = net(inputs)
            loss = criterion(outputs, labels)
            print ("----------out---------")
            print (outputs)
            print ("----------labels---------")
            print (labels)
            print ("---------params---------")
            print (loss.item())
            loss.backward()
            optimizer.step()

            # print statistics
            running_loss += loss.item()
#            if i % 2000 == 1999:    # print every 2000 mini-batches
            if True:    # print every 2000 mini-batches
                print('[%d, %5d] loss: %.3f' %
                      (epoch + 1, i + 1, running_loss / 2000))
                running_loss = 0.0

    print('Finished Training')



    #real,predicted = list(generate_sample())
    #draw_matrix(real,predicted)
    #sample2vector(real)
    #print (real.head())
    #print (predicted.head())

