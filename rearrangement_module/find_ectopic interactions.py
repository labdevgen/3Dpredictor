import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import pickle

def find_ectopic_interactions(diag_array):
    ectopic_array = np.zeros_like(Mut_array)
    ectopic_array = ectopic_array + np.NaN
    percs=[]
    #################
    # values = diff_array[np.triu_indices(len(diff_array))].flatten()
    # values_normal = values[~np.isnan(values)]
    # values_normal = values_normal[np.logical_and(values_normal<np.percentile(values_normal,96),
    #                                              values_normal>np.percentile(values_normal, 4))]
    # values_normal = values_normal[np.nonzero(values_normal)]
    # mean = np.mean(values_normal)
    # std = np.std(values_normal)
    # print(mean, std)
    #
    # ectopic_ids = np.abs(diff_array) - mean > 2*std
    # #print (ectopic_ids)
    # ectopic_array[ectopic_ids] = diff_array[ectopic_ids]
    # return ectopic_array

    ###################

    for k,k_array in enumerate(diag_array):
        # print(k_array)
        nonzero = np.nonzero(k_array)
        # print("nonzero", nonzero)
        perc_upper = np.percentile(k_array[nonzero], 96)
        percs.append(perc_upper)
        # print(perc_upper)
        perc_bottom = np.percentile(k_array[nonzero], 4)
        # print(perc_bottom)
        # print(k_array)
        diag_96perc = [k for k in k_array if k < perc_upper  and k > perc_bottom and k != 0]
        if len(diag_96perc) < 10:
            continue

        # print(len(k_array), len(diag_96perc))
        # print(diag_96perc)
        diag_mean = np.mean(diag_96perc)
        diag_std = np.std(diag_96perc)
        print("diag_mean ", diag_mean, " diag_std", diag_std)
        count = 0
        for n, contact in enumerate(k_array):
            #if abs(contact-diag_mean) > 2*diag_std and contact != 0:
            if contact != 0:
                # print((contact-diag_mean)/diag_std)
                ectopic_array[k+n, n]=(contact-diag_mean)/diag_std
                count += 1
        print (count,len(ectopic_array)-k,count / (len(ectopic_array)-k))
    # print(np.min(percs), np.max(percs))
    return ectopic_array
def get_all_digonals(array_a):
    diagonals=[]
    for k in range(len(array_a)):
        diag=np.diag(array_a, k=k)
        diagonals.append(diag)
    return diagonals

# dataWT = pd.read_csv("input/Dactily/chr1.5MB.MusWT.contacts", sep="\t")
# dataMut = pd.read_csv("input/Dactily/chr1.5MB.MusDelB.contacts", sep="\t")
dataWT = pd.read_csv("input/Dactily/chr1.5MB.MusWT.contacts", sep="\t")
dataMut = pd.read_csv("input/Dactily/chr1.5MB.MusDelB.contacts", sep="\t")
int_st = 76400000
int_en = 78075000
pivotWT = dataWT.pivot(index="st", columns="end", values="count")
pivotWT.fillna(value=0, inplace=True)
pivotMut=dataMut.pivot(index="st", columns="end", values="count")

sortPivotMut = pivotMut[sorted(pivotMut.keys())]
sortPivotMut.fillna(value=0, inplace=True)
sortPivotMut.sort_index(inplace=True)
cols_to_drop = []

for n,key in enumerate(sortPivotMut.keys()):
    assert key == pivotWT.keys()[n]
    if int_st <= key <=int_en:
        sortPivotMut[key] = [0]*len(sortPivotMut)
        sortPivotMut.loc[key,:] = [0]*len(sortPivotMut)
        pivotWT[key]=[0]*len(sortPivotMut)
        pivotWT.loc[key, :] = [0] * len(pivotWT)
        cols_to_drop.append(key)
    if n!=0:
        if sortPivotMut.keys()[n]-sortPivotMut.keys()[n-1]!=25000:
            print(sortPivotMut.keys()[n],sortPivotMut.keys()[n-1])
            print(pivotWT.keys()[n], pivotWT.keys()[n - 1])
            # print(max(pivotWT.keys()))
            # print(min(pivotWT.keys()))
        elif pivotWT.keys()[n] - pivotWT.keys()[n-1]!=25000:
            print("wt")
            print(pivotWT.keys([n]), pivotWT.keys()[n-1])
        # assert sortPivotMut.keys()[n]-sortPivotMut.keys()[n-1]==25000
for n,key in enumerate(sortPivotMut.index):
    assert key == pivotWT.index[n]
# print(cols_to_drop)

# To take into account the genomic distance bias, we normalized
# the difference matrix by dividing each sub-diagonal by the average wt reads count at its
# corresponding pairwise genomic distance.
Mut_del=sortPivotMut.drop(columns=cols_to_drop, index=cols_to_drop)
WT_del=pivotWT.drop(columns=cols_to_drop, index=cols_to_drop)
WT_sum = sum(WT_del.sum())
Mut_sum=sum(Mut_del.sum())
# print(WT_sum, Mut_sum)
WT_array = pivotWT.as_matrix()
WT_diags = get_all_digonals(WT_array)


# print(len(WT_diags),WT_diags[6])
Mut_array = sortPivotMut.as_matrix()
read_coef = WT_sum/Mut_sum
Mut_array=Mut_array*read_coef
diff_array=Mut_array-WT_array

# print(np.diag(diff_array))

for i in range(len(diff_array)):
    X = np.array(range(0,len(diff_array)-i))
    Y = np.array(range(i,len(diff_array)))
    coeff = np.average(WT_array[X,Y])
    diff_array[X,Y] = diff_array[X,Y] / coeff
    diff_array[Y,X] = diff_array[X,Y]
print("diff_array")
print(diff_array)
#plt.imshow(diff_array)
#plt.colorbar()
#plt.draw()
# plt.show()

diff_diags = get_all_digonals(diff_array)
# diff_diags_norm = []
# for k,value in enumerate(diff_diags):
#     # if k<5:
#         # print(k)
#         # print(np.mean(WT_diags[k]))
#         # print(value)
#         # print(value/np.mean(WT_diags[k]))
#     x = value/np.mean(WT_diags[k])
#     diff_diags_norm.append(x)
# assert len(diff_diags_norm)==len(diff_diags)
np.nan_to_num(diff_diags, copy=False)
ectopic_array = find_ectopic_interactions(diff_diags)
print(np.min(ectopic_array))
print(np.max(ectopic_array))

# plt.imshow(ectopic_array, cmap="OrRd")
# plt.show()
#
# print(len(WT_array))
WT_array=np.log(WT_array)
plt.title("WT_array")
plt.imshow(WT_array, cmap="OrRd")
# plt.imshow(ectopic_array)
plt.show()
# #
# Mut_array=np.log(Mut_array)
# plt.imshow(Mut_array, cmap="bwr")
# plt.imshow(ectopic_array)
# plt.show()
plt.title("diff_array")
plt.imshow(diff_array, cmap="OrRd")
plt.show()
plt.clf()
# diff_array = np.log(diff_array)
# colored_data = plt.imshow(diff_array, cmap="OrRd")
with open("out/analysis/rearrangement/ect_array_real.pickle", 'wb') as f:
    pickle.dump(ectopic_array, f)
ectopic_array[ectopic_array<=2] = np.NaN
plt.title("ectopic_interactions_real")
plt.imshow(ectopic_array, cmap="OrRd")
# plt.colorbar(colored_data)
# plt.imshow(ectopic_array)

plt.show()




# plt.imshow(WT_array[0:199], cmap="OrRd")
# plt.scatter([1,1],[2,1])

# print(diagonal_array[1])
# print(len(diagonal_array), len(diff_array))


# print(pivotWT)
# print(WT_array)
# print("---------------------------")
# print(sortPivotMut)
# print(Mut_array)
# print(len(pivotWT), len(pivotMut))
# print(len(pivotWT.keys()), len(pivotMut.keys()))
# print(pivotMut.loc[70950000, 70950000])
# print(pivotWT.loc[70950000, 70950000])
# print(sorted(pivotMut.keys()))
# print(pivotMut[sorted(pivotMut.keys())])
# print(sortPivotMut.index[0])
# # print(dataWT.keys())
