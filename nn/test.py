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

def draw_dataset(dataset,net,subset = 1, title=""):
    def _draw_dataset(dataset,net,title):
        #plt.clf()
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
            plt.imshow(r_reduced,vmin=total_min,vmax=total_max, cmap="OrRd")
            plt.axis('off')
            if i==0: plt.title("Real")
            plt.subplot(N, 3, (i+1)*3 + 2)
            plt.imshow(p_reduced,vmin=total_min,vmax=total_max, cmap="OrRd")
            plt.axis('off')
            if i == 0: plt.title("Predicted")
            plt.subplot(N, 3, (i+1)*3 + 3)
            p = torch.reshape(p,tuple([1]+list(p.shape))) # make +one dimension to imitate batch
            nn_p = net(p)
            nn_p = nn_p.reshape(nn_p.shape[-2], nn_p.shape[-1]).detach().numpy()
            plt.imshow(nn_p,vmin=total_min,vmax=total_max, cmap="OrRd")
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

def train_and_show(validation_dataset, train_dataloader,net,num_epochs,title="", subset = 1, lr = 0.05,
                   interactive=True):
    #optimizer = optim.SGD(net.parameters(), lr=lr)
    optimizer = optim.Adagrad(net.parameters(),lr=lr)
    criterion = nn.MSELoss()

    print (datetime.datetime.now(), " Start iteration")
    calc_raw_losses = True
    start_time = datetime.datetime.now()

    epoch = 0
    while epoch < num_epochs:
        raw_losses = []
        losses = []
        for i_batch, (p,r) in enumerate(train_dataloader):
            optimizer.zero_grad()
            net_out = net(p)
            if calc_raw_losses:
                raw_losses.append(criterion(p, r).item())
            loss = criterion(net_out, r)
            loss.backward()
            optimizer.step()
            losses.append(loss.item())
        time_delta = datetime.datetime.now() - start_time
        if calc_raw_losses:
            calc_raw_losses = False
            print("Raw loss = ", np.average(raw_losses))
        if (epoch % 100 == 0 and epoch >= 100) or time_delta.seconds > 60:
                print(datetime.datetime.now()," ", epoch, " loss = ", np.average(losses))
                start_time = datetime.datetime.now()
        epoch += 1
        if epoch == num_epochs:
            if interactive:
                print(datetime.datetime.now(), " ", epoch, " loss = ", np.average(losses))
                draw_dataset(validation_dataset, net, title=title + " num epochs = " + str(num_epochs) + " lr = " + str(lr),
                             subset=subset)
                add_epoch = input("How many inputs more to compute? (enter 0 to stop iterating) \n")
                try:
                    num_epochs += int(add_epoch)
                except:
                    add_epoch = input("Please eneter integer value\n")
                    num_epochs += int(add_epoch)

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
    train_and_show(validation_dataset = validation_dataset,train_dataloader=dataloader,net=net,num_epochs=6,
                   title=str(net) + "\nMoving triangle with or w/o loop",subset = 5, lr = 0.05)
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
    train_and_show(validation_dataset = validation_dataset,train_dataloader=dataloader,net=net,num_epochs=800,
                   title=str(net) + "\nMoving triangle with or w/o loop",subset = 5, lr=0.01)

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

def moving_triangle_and_flying_loop_and_Noize_convOnlyNet():
    size = (40,40)
    dataset = DatasetMaker_moving_TAD_with_flying_Loop_and_Noize(numSamples=80,size=size,st=1,en=20,maxlen=12,TAD=3)
    validation_dataset = DatasetMaker_moving_TAD_with_flying_Loop_and_Noize(numSamples=10, size=size, st=19, en=34, maxlen=13, TAD=5)
    dataloader = DataLoader(dataset,batch_size=4)
    net = SimpleConvOnlyNet(input_size=size,output_size=size)
    train_and_show(validation_dataset = validation_dataset,train_dataloader=dataloader,net=net,num_epochs=1250,
                   title=str(net) + "\nMoving triangle with or w/o loop",subset = 5, lr=0.01)
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

    row = 2
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

    weights2 = np.zeros(shape=(filter_size,filter_size))
    weights2[3,3] = 1


    plt.subplot(row,column,2)
    plt.imshow(weights)
    plt.title("weights")
    net.conv1.weight.data[0] = torch.from_numpy(weights).float().unsqueeze(0)
    net.conv1.bias.data[0] = torch.from_numpy(np.array([-14])).float().unsqueeze(0)

    net.conv2.weight.data[0] = torch.from_numpy(weights2).float().unsqueeze(0)
    net.conv2.bias.data[0] = torch.from_numpy(np.array([0])).float().unsqueeze(0)

    b = net.conv1(a).detach()
    plt.subplot(row,column,3)
    plt.title("conv1")
    plt.imshow(b[0][0])

    import torch.nn.functional as F
    c = F.relu(b)
    plt.subplot(row,column,4)
    plt.title("conv1 + relu")
    plt.imshow(c[0][0])

    d = net.conv2(c).detach()
    plt.subplot(row,column,5)
    plt.title("conv2")
    plt.imshow(d[0][0])

    plt.subplot(row,column,6)
    plt.title("weights")
    plt.imshow(weights2)


    #plt.gca().grid(which='minor', color='w', linestyle='-', linewidth=2)

    plt.show()

def test_real_data():
    dataset = DatasetFromRealAndPredicted("input/model5823236996.validation..equal.h=2.Interval_chr9_0_141100000validatingOrient.contacts.gz.6.1500000.50001.863649.25000.txt.scc")
    size = dataset[0][0][0].size() # [first_sample_in_dataset][predicted][first_color]
    validation_dataset = dataset
    dataloader = DataLoader(dataset,batch_size=12)
    #net = ConvNet_1(input_size=size,output_size=size)
    net = ConvNet_2(input_size=size,output_size=size)

    train_and_show(validation_dataset = validation_dataset,train_dataloader=dataloader,net=net,num_epochs=2,
                   title=str(net) + "\nReal data",subset = 5, lr=0.005)


if __name__ == "__main__":
    #freeze_support()
    #fixed_placed_triangle()
    #fixed_placed_triangle_with_random_noize()
    #fixed_placed_triangle_with_random_noize_and_loop()
    #moving_triangle_with_random_noize_and_loop()
    #moving_triangle_with_random_noize_and_loop_convNet()
    #moving_triangle_and_flying_loop_and_Noize_convOnlyNet()
    #test2()
    test_real_data()