import sys, os
import math
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import sys
sys.path.append(os.path.dirname(os.getcwd()))
from matrix_plotter import MatrixPlotter
from test_ds_creators import *
from test_nets import *
from shared import Interval
import datetime
from multiprocessing import  freeze_support

# Here comes torch
import torch
import torch.nn as nn
import torch.optim as optim
import torch.nn.functional as F
from torch.utils.data import Dataset, DataLoader


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

def draw_dataset(dataset,net,subset = 1, title=""):
    def _draw_dataset(dataset,net,title):
        N = len(dataset)
        mi = np.finfo(np.float64).max
        ma = 0
        for i in range(len(dataset)):
            p,r = dataset[i]
            mi = min(mi,np.min(p.detach().numpy()),np.min(r.detach().numpy()))
            ma = max(ma, np.max(p.detach().numpy()), np.max(r.detach().numpy()))
        total_min, total_max = mi, ma
        N += 1 # one more line for major title
        plt.subplot(N,3,2)
        plt.title(title)
        plt.axis('off')

        for i in range(N-1):
            p,r = dataset[i]
            p_reduced = p.reshape(p.shape[-2],p.shape[-1]).detach().numpy()
            r_reduced = r.reshape(r.shape[-2],r.shape[-1]).detach().numpy()
            plt.subplot(N, 3,(i+1)*3 + 1)
            plt.imshow(r_reduced,vmin=total_min,vmax=total_max)
            plt.axis('off')
            if i==0: plt.title("Real")
            plt.subplot(N, 3, (i+1)*3 + 2)
            plt.imshow(p_reduced,vmin=total_min,vmax=total_max)
            plt.axis('off')
            if i == 0: plt.title("Predicted")
            plt.subplot(N, 3, (i+1)*3 + 3)
            p = torch.reshape(p,tuple([1]+list(p.shape))) # make +one dimension to imitate batch
            nn_p = net(p)
            nn_p = nn_p.reshape(nn_p.shape[-2], nn_p.shape[-1]).detach().numpy()
            plt.imshow(nn_p,vmin=total_min,vmax=total_max)
            plt.axis('off')
            if i == 0: plt.title("Predicted + NN")
        plt.show()

    if subset > 1:
        subset = min(1,subset / len(dataset))
    elif subset <= 0:
        raise Exception("Wrong value for subset")

    N = int(len(dataset)*subset)
    subset_of_dataset = torch.utils.data.Subset(dataset,np.random.choice(len(dataset),size=N))
    print (subset, " Subsetting to ",len(subset_of_dataset))
    _draw_dataset(subset_of_dataset, net, title=title)

def plot_matrixes(arrs,showFig=True,**kwargs):
    if "titles" in kwargs:
        titles = kwargs["titles"]
    else:
        titles = list(map(str,range(arrs)))

    total_min = min(map(np.min,arrs))
    total_max = max(map(np.max,arrs))
    for ind,val in enumerate(arrs):
        plt.subplot(len(arrs),1,ind+1)
        plt.imshow(val,vmin=total_min,vmax=total_max)
        plt.title(titles[ind])
    if showFig:
        plt.show()
        plt.clf()

def train_and_show(validation_dataset, train_dataloader,net,num_epochs,title="", subset = 1):
    lr = 0.05
    #optimizer = optim.SGD(net.parameters(), lr=lr)
    optimizer = optim.Adagrad(net.parameters(),lr=lr)
    criterion = nn.MSELoss()

    for epoch in range(num_epochs):
        raw_losses = []
        losses = []
        for i_batch, (p,r) in enumerate(train_dataloader):
            optimizer.zero_grad()
            net_out = net(p)
            raw_losses = criterion(p, r)
            loss = criterion(net_out, r)
            loss.backward()
            optimizer.step()
            losses.append(loss.item())
        if (epoch % 100 == 0 and epoch >= 100):
                if epoch == 100:
                    print("Raw loss = ",np.average(raw_losses))
                print (datetime.datetime.now()," ", epoch, " loss = ", np.average(losses))
    print(datetime.datetime.now(), " ", epoch, " loss = ", np.average(losses))
    draw_dataset(validation_dataset,net,title=title + " num epochs = "+str(num_epochs) + " lr = " + str(lr), subset=subset)

def fixed_placed_triangle():
    size = (10,10)
    dataset = DatasetMaker(numSamples=4,size=size,st=2,en=5,TAD=2)
    dataloader = DataLoader(dataset,batch_size=2)
    net = SimplestNet(input_size=size,output_size=size)
    train_and_show(dataset,dataloader=dataloader,net=net,num_epochs=1000,
                   title="Fixed place triangle")

def fixed_placed_triangle_with_random_noize():
    size = (10,10)
    dataset = DatasetMaker_normalNoize(numSamples=4,size=size,st=2,en=5,TAD=3)
    dataloader = DataLoader(dataset,batch_size=2)
    net = SimplestNet(input_size=size,output_size=size)
    train_and_show(dataset,dataloader=dataloader,net=net,num_epochs=1000,
                   title="Fixed place triangle with noize")

def fixed_placed_triangle_with_random_noize_and_loop():
    size = (20,20)
    dataset = DatasetMaker_normalNoize_withLoop(numSamples=8,size=size,st=2,en=14,TAD=3)
    dataloader = DataLoader(dataset,batch_size=2)
    net = SimplestNet(input_size=size,output_size=size)
    train_and_show(dataset,dataloader=dataloader,net=net,num_epochs=1000,
                   title="Fixed place triangle with noize and loop",subset = 5)


def moving_triangle_with_random_noize_and_loop():
    size = (40,40)
    dataset = DatasetMaker_moving_TAD(numSamples=40,size=size,st=2,en=38,maxlen=15,TAD=3)
    dataloader = DataLoader(dataset,batch_size=2)
    net = SimplestNet(input_size=size,output_size=size)
    train_and_show(dataset,dataloader=dataloader,net=net,num_epochs=1000,
                   title="Moving triangle with noize and loop",subset = 5)

def moving_triangle_and_loop_convOnlyNet():
    size = (40,40)
    dataset = DatasetMaker_moving_TAD(numSamples=80,size=size,st=1,en=25,maxlen=12,TAD=3)
    validation_dataset = DatasetMaker_moving_TAD(numSamples=10, size=size, st=25, en=35, maxlen=15, TAD=3)
    dataloader = DataLoader(dataset,batch_size=4)
    net = SimpleConvOnlyNet(input_size=size,output_size=size)
    train_and_show(validation_dataset = validation_dataset,train_dataloader=dataloader,net=net,num_epochs=600,
                   title=str(net) + "\nMoving triangle with or w/o loop",subset = 5)

    p,r = validation_dataset.__getitem__(0)
    print (p)
    print (p.size())
    plt.clf()
    plt.subplot(5,3,1)
    plt.imshow(p[0].detach().numpy()) # predicted figure
    p = p.unsqueeze(0)
    print (p)
    print(p.size())
    plt.subplot(5,3,2)
    plt.imshow(r[0].detach().numpy()) # original figure
    pnn = net(p)
    plt.subplot(5,3,3)
    plt.imshow(pnn[0][0].detach().numpy()) # original figure

    c1 = net.conv1(p)
    print (c1.size())
    c2 = F.relu(c1)
    for i in range(c1.size()[1]):
        plt.subplot(5, 3, 4 + 3*i)
        plt.imshow(c1.detach().numpy()[0][i])
        plt.subplot(5, 3, 5 + 3*i)
        plt.imshow(c2.detach().numpy()[0][i])
        plt.subplot(5,3, 6 + 3*i,)
        plt.imshow(net.conv1.weight.data[i][0].numpy())
    plt.show()
    #print(net.conv1.weight.data[0][0])
    #plt.imshow(net.conv1.weight.data.numpy()[0,0])
    #plt.show()

def moving_triangle_and_flying_loop_convOnlyNet():
    size = (40,40)
    dataset = DatasetMaker_moving_TAD_with_flying_Loop(numSamples=80,size=size,st=1,en=20,maxlen=12,TAD=3)
    validation_dataset = DatasetMaker_moving_TAD_with_flying_Loop(numSamples=10, size=size, st=19, en=34, maxlen=13, TAD=3)
    dataloader = DataLoader(dataset,batch_size=4)
    net = SimpleConvOnlyNet(input_size=size,output_size=size)
    train_and_show(validation_dataset = validation_dataset,train_dataloader=dataloader,net=net,num_epochs=600,
                   title=str(net) + "\nMoving triangle with or w/o loop",subset = 5)

    p,r = validation_dataset.__getitem__(0)
    print (p)
    print (p.size())
    plt.clf()
    nrows = (len(net.conv1.weight.data)+2)
    ncols = 3

    plt.subplot(nrows, ncols,1)
    plt.imshow(p[0].detach().numpy()) # predicted figure
    p = p.unsqueeze(0)
    print (p)
    print(p.size())
    plt.subplot(nrows, ncols,2)
    plt.imshow(r[0].detach().numpy()) # original figure
    pnn = net(p)
    plt.subplot(nrows, ncols,3)
    plt.imshow(pnn[0][0].detach().numpy()) # original figure

    c1 = net.conv1(p)
    print (c1.size())
    c2 = F.relu(c1)
    for i in range(c1.size()[1]):
        plt.subplot(nrows, ncols, 4 + 3*i)
        plt.imshow(c1.detach().numpy()[0][i])
        plt.subplot(nrows, ncols, 5 + 3*i)
        plt.imshow(c2.detach().numpy()[0][i])
        plt.subplot(nrows, ncols, 6 + 3*i,)
        plt.imshow(net.conv1.weight.data[i][0].numpy())
    plt.show()
    #print(net.conv1.weight.data[0][0])
    #plt.imshow(net.conv1.weight.data.numpy()[0,0])
    #plt.show()


def test():
    size = (40,40)
    dataset = DatasetMaker_moving_TAD_with_flying_Loop(numSamples=80,size=size,st=1,en=20,maxlen=12,TAD=3)
    validation_dataset = DatasetMaker_moving_TAD_with_flying_Loop(numSamples=10, size=size, st=19, en=34, maxlen=13, TAD=3)
    dataloader = DataLoader(dataset,batch_size=4)
    net = SimpleConvOnlyNet(input_size=size,output_size=size)

    p,r = validation_dataset.__getitem__(0)
    print (p)
    print (p.size())
    plt.clf()
    nrows = (len(net.conv1.weight.data)+2)
    ncols = 3

    plt.subplot(nrows, ncols,1)
    plt.imshow(p[0].detach().numpy()) # predicted figure
    p = p.unsqueeze(0)
    print (p)
    print(p.size())
    plt.subplot(nrows, ncols,2)
    plt.imshow(r[0].detach().numpy()) # original figure

    weight = np.zeros(shape=(5,5))
    weight[4,4] = 1
    #weight[1, :] = 1
    #bias = np.zeros(shape=(5,5))
    print (net.conv1.weight.data.size())
    net.conv1.weight.data[0][0] = torch.torch.from_numpy(weight).float()
    #net.conv1.bias.data[0][0] = torch.torch.from_numpy(bias).float()


    pnn = net.forward(p)
    plt.subplot(nrows, ncols,3)
    plt.imshow(pnn[0][0].detach().numpy()) # original figure

    c1 = net.conv1(p)
    print (c1.size())
    c2 = F.relu(c1)


    for i in range(c1.size()[1]):
        plt.subplot(nrows, ncols, 4 + 3*i)
        plt.imshow(c1.detach().numpy()[0][i])
        plt.subplot(nrows, ncols, 5 + 3*i)
        plt.imshow(c2.detach().numpy()[0][i])
        plt.subplot(nrows, ncols, 6 + 3*i,)
        plt.imshow(net.conv1.weight.data[i][0].numpy())
    plt.show()
    print(net.conv1.weight.data[0][0])
    #plt.imshow(net.conv1.weight.data.numpy()[0,0])
    #plt.show()


def test2():
    a = np.zeros(shape=(10,10))
    i,j = np.triu_indices(3)
    st = 5
    i += st
    j += st
    a[i,j] = 3
    a[j,i] = 3
    a[st-3,st-2] = 5
    a[st-2,st-3] = 5
    a = torch.from_numpy(a).float()
    a = a.unsqueeze(0).unsqueeze(0)
    print(a)

    net = test_conv(nfilters = 1, inshape=10)

    row = 1
    column = 4
    plt.subplot(row,column,1)
    plt.imshow(a[0][0])

    filter_size = 5
    weights = np.zeros(shape=(filter_size,filter_size))
    #weights[:,0] = 1
    #weights[:,1] = -1
    #weights[-1,:] = 1
    #weights[-2,1:] = -1
    #weights[:,0] = 1
    #weights[2, 3] = 1
    #weights[2, 1] = -1
    #weights[2, 2] = -1
    #weights[1, 3] = -1
    #weights[0,0] = 1
    #weights[-2,:] = -1
    weights[:,-1] = +1
   # weights[:,-2] = -1
    weights[-1,:] = 1

    weights2 = np.zeros(shape=tuple(net.linear.weight.data.size()))
    weights2[np.diag_indices(len(weights2))] = 1
    net.linear.weight.data = torch.from_numpy(weights2).float()
    net.linear.bias.data = torch.from_numpy(np.zeros(shape=tuple(net.linear.bias.size()))).float()


    plt.subplot(row,column,2)
    plt.imshow(weights)
    plt.title("weights")
    net.conv1.weight.data[0] = torch.from_numpy(weights).float().unsqueeze(0)
    net.conv1.bias.data[0] = torch.from_numpy(np.array([-14])).float().unsqueeze(0)

    b = net.forward(a).detach()
    plt.subplot(row,column,3)
    plt.title("result")
    plt.imshow(b[0][0])

    import torch.nn.functional as F
    c = F.relu(b)
    plt.subplot(row,column,4)
    plt.title("result + relu")
    plt.imshow(c[0][0])


    plt.gca().grid(which='minor', color='w', linestyle='-', linewidth=2)

    plt.show()

if __name__ == "__main__":
    #freeze_support()
    #fixed_placed_triangle()
    #fixed_placed_triangle_with_random_noize()
    #fixed_placed_triangle_with_random_noize_and_loop()
    #moving_triangle_with_random_noize_and_loop()
    #moving_triangle_with_random_noize_and_loop_convNet()
    #moving_triangle_and_flying_loop_convOnlyNet()
    test()