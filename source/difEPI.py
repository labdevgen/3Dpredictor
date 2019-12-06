import numpy as np
import scipy.stats as sc
import matplotlib.pyplot as plt
import func
from DatadifEPI import *
reload(func)
c = []

sp = 'mouse5K'

c = func.OEdumping(EPP,c1,c2,c1m,c2m,c1r1,c1r2)


plt.cla()
lab=['hep_vs_npc','hep_vs_model1','npc_vs_model2','rep1_vs2'] # name of 
# lab=['Mcyte_vs_K562','Mcyte_vs_model1','K562_vs_model2','Mcyte_vs_model2','K562_vs_model1','model1_vs_model2']
print lab
# cc = c[:,1]-c[:,0],c[:,2]-c[:,0],c[:,3]-c[:,1],c[:,3]-c[:,0],c[:,2]-c[:,1],c[:,3]-c[:,2],c[:,5]-c[:,4]
# acc = np.abs(c[:,1]-c[:,0]),np.abs(c[:,2]-c[:,0]),np.abs(c[:,3]-c[:,1]),np.abs(c[:,3]-c[:,0]),np.abs(c[:,2]-c[:,1]),np.abs(c[:,3]-c[:,2]),np.abs(c[:,5]-c[:,4])
# print np.nanmedian(acc,axis=1)
#print np.nanmedian(cc,axis=1)

# plt.violinplot(cc,showmedians=True,showextrema=False)#,showfliers=False)
# plt.xticks(range(1,8),lab,rotation=15,size=8)
# plt.ylim([-6,6])
# plt.savefig('%s_viol_all' % sp,dpi=400)
# plt.clf()
# plt.violinplot(acc,showmedians=True,showextrema=False)#,showfliers=False)
# plt.xticks(range(1,8),lab,rotation=15,size=8)
# plt.ylim([0,5])
# plt.savefig('%s_viol_all_abs' % sp,dpi=400)


lim = 1
ch = c[((np.abs(c[:,1]-c[:,0]) > lim)&(np.abs(c[:,5]-c[:,4]) < lim)),:]
acc = np.abs(ch[:,1]-ch[:,0]),np.abs(ch[:,2]-ch[:,0]),np.abs(ch[:,3]-ch[:,1]),np.abs(ch[:,5]-ch[:,4])
st = sc.mannwhitneyu(acc[0],acc[1]),sc.mannwhitneyu(acc[0],acc[2]),sc.mannwhitneyu(acc[1],acc[2])
print np.nanmedian(acc,axis=1)
print st
# print np.nanmedian(cc,axis=1)
# plt.clf()
# plt.violinplot(cc,showmedians=True,showextrema=False)#,showfliers=False)
# plt.xticks(range(1,8),lab,rotation=15,size=8)
# plt.ylim([-6,6])
# plt.savefig('%s_viol_dif_rep_1' % sp,dpi=400)
# plt.clf()
# plt.violinplot(acc,showmedians=True,showextrema=False)#,showfliers=False)
# plt.xticks(range(1,8),lab,rotation=15,size=8)
# plt.ylim([0,5])
# plt.savefig('%s_viol_dif_rep_1_abs' % sp,dpi=400)