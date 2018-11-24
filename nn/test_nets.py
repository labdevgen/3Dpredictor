import numpy as np
import torch
import torch.nn as nn
import torch.nn.functional as F

class SimplestNet(nn.Module):
    # Lets make a simple Net, 1
    def __init__(self,input_size,output_size):
        super(SimplestNet, self).__init__()
        mult = 1
        for i in input_size:
            mult *= i
        self.l_inp = nn.Linear(mult,100)

        mult = 1
        for i in output_size:
            mult *= i
        self.l_hid = nn.Linear(100,mult)

        self.out_size = output_size


    def forward(self, x):
        #print ("-----")
        #print(x.shape)
        batch_size = x.shape[0]
        x = torch.reshape(x,(batch_size,1,-1)) # reshape to make linear seq
        #print (x.shape)
        x = F.relu(self.l_inp(x))
        #print(x.shape)
        x = self.l_hid(x)
        #print(x.shape)
        x = torch.reshape(x,tuple([batch_size]+list(self.out_size))) # reshape to requered size
        #print(x.shape)
        return x

# 2D image --> 1 conv layer --> reshape ---> fully connected layer size 100 --> linear output layer --> reshape
class SimpleConvNet(nn.Module):
    def __init__(self, input_size, output_size):
        super(SimpleConvNet, self).__init__()
        filter_size = 4
        n_filters = 1
        self.conv1 = nn.Conv2d(1, n_filters, filter_size) # 3x3 square to find TAD loop
        self.triu_ids = np.triu_indices(input_size[0]-filter_size + 1,input_size[1]-filter_size+1)

        self.linear_size = (input_size[0]-filter_size + 1)*(input_size[1]-filter_size+1)*n_filters
        self.linear = nn.Linear(self.linear_size,self.linear_size)
        self.out = nn.Linear(self.linear_size,output_size[0]*output_size[1])
        self.out_shape = list(output_size)

    def forward(self, x):
        #print ("------")
        batch_size = x.shape[0]
        #print (x.shape)
        x = self.conv1(x)
        x = F.relu(x)
        #print (x.shape)
        #print (self.linear_size)
        x = x.reshape((batch_size,1,self.linear_size))
        #print (x.shape)
        x = self.linear(x)
        x = F.relu(x)
        x = self.out(x)
        x = x.reshape(tuple([batch_size,1]+self.out_shape))
        return x

class SimpleConvOnlyNet(nn.Module):
    def __init__(self, input_size, output_size):
        super(SimpleConvOnlyNet, self).__init__()
        filter_size = 3
        n_filters = 4
        self.conv1 = nn.Conv2d(1, n_filters, filter_size) # 3x3 square to find TAD loop
        self.conv2 = nn.Conv2d(n_filters,1,1)
        self.linear_size = (input_size[0]-filter_size + 1)*(input_size[1]-filter_size+1)
        assert (output_size[0]*output_size[1] + output_size[0]) % 2 == 0
        self.out = nn.Linear(self.linear_size, (output_size[0]*output_size[1] + output_size[0]) // 2)
        self.out_shape = list(output_size)
        self.triu_i, self.triu_j = np.triu_indices(output_size[0])
    def forward(self, x):
        #print ("------")
        batch_size = x.shape[0]
        #print (x.shape)
        x = self.conv1(x)
        x = F.relu(x)
        #print (x.shape)
        #print (self.linear_size)
        #x = x.reshape((batch_size,1,self.linear_size))
        #print (x.shape)
        x = self.conv2(x)
        #x = F.relu(x)
        #print(x.shape)
        x = x.reshape((batch_size,-1))
        #print(x.shape)
        #print (self.linear_size)
        x = self.out(x)
        #print (x.shape)
        y = torch.empty(([batch_size,1] + self.out_shape),dtype=x.dtype)
        #print (y.shape)
        #print (len(self.triu_i),len(self.triu_j))
        #print (y[0].shape)
        for i in range(batch_size):
            y[i][0][self.triu_i,self.triu_j] = x[i]
            y[i][0][self.triu_j,self.triu_i] = x[i]
        #x = x.reshape(tuple([batch_size,1]+self.out_shape))
        return y


class test_conv(nn.Module):
    def __init__(self, nfilters, inshape):
        super(test_conv, self).__init__()
        filter_size = 5
        n_filters = 1
        self.conv1 = nn.Conv2d(1, n_filters, filter_size, padding=2)
        self.linear = nn.Linear(inshape**2,inshape**2)
    def forward(self, x):
        batch = x.size()[0]
        original_shape = x.size()
        x = self.conv1(x)
        x = F.relu(x)
        x = x.reshape((batch,-1))
        x = self.linear(x)
        x = x.reshape(original_shape)
        return x

class ConvNet_1(nn.Module):
    def __init__(self, input_size, output_size):
        super(ConvNet_1, self).__init__()
        filter_size = 5
        n_filters = 10
        assert filter_size % 2 == 1
        assert filter_size > 2
        padding = (filter_size - 1 ) // 2
        self.conv1 = nn.Conv2d(1, n_filters, filter_size, padding=padding) # 3x3 square to find TAD loop
        self.conv2 = nn.Conv2d(n_filters, n_filters, filter_size, padding=padding)
        self.conv3 = nn.Conv2d(n_filters, n_filters // 2, filter_size, padding=padding)
        self.conv4 = nn.Conv2d(n_filters // 2, 1, filter_size, padding=padding)
    def forward(self, x):
        x1 = self.conv1(x)
        x1 = F.relu(x1)
        x1 = self.conv2(x1)
        x1 = F.relu(x1)
        x1 = self.conv3(x1)
        x1 = F.relu(x1)
        return self.conv4(x1)

class ConvNet_2(nn.Module):
    def __init__(self, input_size, output_size):
        super(ConvNet_2, self).__init__()
        self.linear = nn.Linear(input_size[0]*input_size[1],input_size[0]*input_size[1]*2)
        self.out = nn.Linear(input_size[0]*input_size[1]*2, output_size[0]*output_size[1])
        self.out_size = output_size
    def forward(self, x):
        x = x.reshape(tuple(x.size()[0:2])+(-1,))
        x = self.linear(x)
        x = F.relu(x)
        return self.out(x).reshape(tuple(x.size()[0:2])+tuple(self.out_size))
