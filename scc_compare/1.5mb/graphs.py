import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib as mpl
import seaborn
th = np.linspace(0, 2*np.pi, 128)

# def demo(sty):
#     mpl.style.use(sty)
#     fig, ax = plt.subplots(figsize=(3, 3))
#
#     ax.set_title('style: {!r}'.format(sty), color='0.001')
#
#     ax.plot(th, np.cos(th), 'c1', label='c1')
#     ax.plot(th, np.sin(th), 'c2', label='c2')
#     ax.legend()
#
# #demo('default')
# demo('seaborn-dark')
# print(mpl.style.available)
# with h5py.File(/home/evgeniy/PycharmProjects/pyth01 + /TP.h5, r+) as h5_file:
#         was_second_intersection = False
#         we_are_above_01 = False
#         action_potential_duration = 0
#
#         fy.append((7*3)*2*3.14/(120*0.5*0.13))
#         fx.append(x)
#         x = x + 6
#
#         fy.append(0)
#         fx.append(x)
#         x = x + 6
#
#         fy.append((8.6*3)*2*3.14/(75*0.5*0.13))
#         fx.append(x)
#         x = x + 6
# d =
# with open('K562_IMR90_1.5mb.txt') as csvfile:
#     spamreader = csv.reader(csvfile, delimiter=' ', quotechar='|')
#     print(spamreader)
#     for row in spamreader:
#         print(np.array(spamre))
df1 = pd.read_csv('K562_IMR90_1.5mb.txt',sep = ' ')
df2 = pd.read_csv('Mcyte_IMR90_1.5mb.txt',sep = ' ' )
df3 = pd.read_csv('Mcyte_K562_1.5mb.txt',sep =  ' ')
df4 = pd.read_csv('Mcyte_NHEK_1.5mb.txt',sep = ' ' )
df5 = pd.read_csv('Mcyte_rep2_1.5mb.txt',sep = ' ' )
df6 = pd.read_csv('NHEK_IMR90_1.5mb.txt',sep =  ' ')
df7 = pd.read_csv('K562_NHEK_1.5mb.txt',sep =  ' ')
df8 = pd.read_csv('rep1_rep2_1.5mb.txt',sep = ' ' )
df9 = pd.read_csv('rep2_K562_1.5mb.txt',sep =  ' ')
df10 = pd.read_csv('rep2_NHEK_1.5mb.txt',sep = ' ' )
df11 = pd.read_csv('rep2_IMR90_1.5mb.txt',sep = ' ' )
df1 = np.array(df1)
df2 = np.array(df2)
df3 = np.array(df3)
df4 = np.array(df4)
df5 = np.array(df5)
df6 = np.array(df6)
df7 = np.array(df7)
df8 = np.array(df8)
df9 = np.array(df9)
df10 = np.array(df10)
df11 = np.array(df11)
x0 = np.array(range(23))
x0 = x0 + 1
x0 = x0/33./23*15
x0 =  x0 - (max(x0) - min(x0))/2 + 1
df = []
#print(list(df1[:,2]) + (df2[:,2]))
df = df + list(df1[:,2]) + list(df2[:,2]) + list(df3[:,2]) + list(df4[:,2]) + list(df5[:,2]) + list(df6 [:,2]) + list(df7[:,2]) + list(df8[:,2]) + list(df9[:,2]) + list(df10[:,2]) + list(df11[:,2])
#print(x)
#print(df1[:,2])
fig=plt.figure()
# fx = list(x)
# fy = list(df1[:,2])
# fig = plt.figure()

graph1 = fig.add_subplot(111)
graph2 = fig.add_subplot(111)
graph3 = fig.add_subplot(111)
graph4 = fig.add_subplot(111)
graph5 = fig.add_subplot(111)
graph6 = fig.add_subplot(111)
graph7 = fig.add_subplot(111)
graph8 = fig.add_subplot(111)
graph9 = fig.add_subplot(111)
graph10 = fig.add_subplot(111)
graph11 = fig.add_subplot(111)
# #graph.plot(fx, fy)
# graph.plot(fx, fy,'rs')
# #print(df2[:,2])
# #print(fx)
# #print(fy)
# #graph1.plot(fx, fy)
# graph1.plot(fx, fy,'bs')



for j in range(6):
    fx = x0[j] + range(23) + 1
    fy = []
    for i in range(23):
        i = i + j*23
        print(i)
        fy = fy + list([df[i]])
    graph1.plot(fx, fy, 's')

fx = x0[j] + range(23) + 1
fy = []
for i in range(23):
    i = i + 7*23
    print(i)
    fy = fy + list([df[i]])
graph1.plot(fx, fy, 'ro')

#ax0=fig.add_subplot(111,frame_on=False)
#ax0.plot(fx,fy, 'bs')

# fx = x0[1] + range(11) + 1
# fy = []
# for i in range(11):
#     i = i*23
#     print([df[i]])
#     fy = fy + list([df[i]])
# graph2.plot(fx, fy, 'bs')
# #ax1=fig.add_subplot(111,frame_on=False)
# #ax1.plot(fx,fy, 'bs')
#
# fx = x0[2] + range(11) + 1
# fy = []
# for i in range(11):
#     i = i*23
#     fy = fy + list([df[i]])
# graph3.plot(fx, fy, 'bs')
# #ax2=fig.add_subplot(111,frame_on=False)
# #ax2.plot(fx,fy, 'bs')
#
# fx = x0[3] + range(11) + 1
# fy = []
# for i in range(11):
#     i = i*23
#     fy = fy + list([df[i]])
# graph4.plot(fx, fy, 'bs')
# #ax3=fig.add_subplot(111,frame_on=False)
# #ax3.plot(fx,fy, 'bs')
#
# fx = x0[4] + range(11) + 1
# fy = []
# for i in range(11):
#     i = i*23
#     fy = fy + list([df[i]])
# graph5.plot(fx, fy, 'bs')
# #ax4=fig.add_subplot(111,frame_on=False)
# #ax4.plot(fx,fy, 'bs')
#
# fx = x0[5] + range(11) + 1
# fy = []
# for i in range(11):
#     i = i*23
#     fy = fy + list([df[i]])
# graph6.plot(fx, fy, 'bs')
# #ax5=fig.add_subplot(111,frame_on=False)
# #ax5.plot(fx,fy, 'bs')
#
# fx = x0[6] + range(11) + 1
# fy = []
# for i in range(11):
#     i = i*23
#     fy = fy + list([df[i]])
# graph7.plot(fx, fy, 'bs')
# #ax6=fig.add_subplot(111,frame_on=False)
# #ax6.plot(fx,fy, 'bs')
#
# fx = x0[7] + range(11) + 1
# fy = []
# for i in range(11):
#     i = i*23
#     fy = fy + list([df[i]])
# graph8.plot(fx, fy, 'bs')
# #ax7=fig.add_subplot(111,frame_on=False)
# #ax7.plot(fx,fy, 'bs')
#
# fx = x0[8] + range(11) + 1
# fy = []
# for i in range(11):
#     i = i*23
#     fy = fy + list([df[i]])
# graph9.plot(fx, fy, 'bs')
# #ax8=fig.add_subplot(111,frame_on=False)
# #ax8.plot(fx,fy, 'bs')
#
# fx = x0[9] + range(11) + 1
# fy = []
# for i in range(11):
#     i = i*23
#     fy = fy + list([df[i]])
# graph10.plot(fx, fy, 'bs')
# #ax9=fig.add_subplot(111,frame_on=False)
# #ax9.plot(fx,fy, 'bs')
#
# fx = x0[10] + range(11) + 1
# fy = []
# for i in range(11):
#     i = i*23
#     fy = fy + list([df[i]])
# graph11.plot(fx, fy, 'bs')
#ax10=fig.add_subplot(111,frame_on=False)
#ax10.plot(fx,fy, 'bs')





# fig = plt.figure()
# graph = fig.add_subplot(111)
#
# fy = [fy[2] for fy in y_epi1]
# fy = np.array(fy)
# graph.plot(range(fy.size), fy)
#
# fy = [fy[2] for fy in y_epi2]
# fy = np.array(fy)
# graph.plot(range(fy.size), fy)
#
# fx = fx + list(x + 1)
# fy = fy + list(df3[:,2])
# graph.plot(fx, fy,'rs')
#
# fx = fx + list(x + 1)
# fy = fy + list(df3[:,2])
# graph.plot(fx, fy,'cs')
#
# fx = fx + list(x + 1)
# fy = fy + list(df3[:,2])
# graph.plot(fx, fy,'ys')
#
# fx = fx + list(x + 1)
# fy = fy + list(df3[:,2])
# graph.plot(fx, fy,'ks')
#
# fx = fx + list(x + 1)
# fy = fy + list(df3[:,2])
#graph.plot(fx, fy,'ws')
#fig = plt.figure()
#graph = fig.add_subplot(111)


# # graph.plot(fx, fy,'C1')
# ax.plot(fx, fy,'bs')
# ax.plot(fx, fy,'C1')
# fy[0] = (17.6*3)*2*3.14/(240*0.5*0.13)
# fy[1] = (7.4*3)*2*3.14/(500*0.5*0.13)
# fy[2] = (7.4*3)*2*3.14/(1000*0.5*0.13)
# graph.plot(fx, fy)
# graph.plot(fx, fy,'bs')
#
# fy[0] = (6.6*3)*2*3.14/(98*0.5*0.13)
# fy[1] = (7.31*3)*2*3.14/(79*0.5*0.13)
# fy[2] = (7.81*3)*2*3.14/(75*0.5*0.13)
# graph.plot(fx, fy,'g')
# graph.plot(fx, fy,'ro')
#
# fy[0] = 0
# fy[1] = (7.26*3)*2*3.14/(80*0.5*0.13)
# fy[2] = (7.37*3)*2*3.14/(75*0.5*0.13)
# graph.plot(fx, fy,'g')
# graph.plot(fx, fy,'bs')
#
#
#
#
# graph.set_xlabel(u'Apical thickness h(mm)')
# graph.set_ylabel(u'Linear ratation speed(mm/s)')
# # plt.hist(l, 25*10/14)
plt.show()