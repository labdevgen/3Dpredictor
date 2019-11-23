import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import pickle
from numpy import inf


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

dataWT = pd.read_csv("MusDelB/5KB/chr1.5KB.DelB.contacts.txt", sep="\t")
dataMut = pd.read_csv("MusDelB/5KB/MusDelB_predicted.txt", sep="\t")
int_st = 76390000
int_en = 78060000
pivotWT = dataWT.pivot(index="st", columns="end", values="count")
pivotWT.fillna(value=0, inplace=True)
pivotMut=dataMut.pivot(index="st", columns="end", values="count")
pivotMut.fillna(value=0, inplace=True)
print(len(set(pivotWT.keys())), len( set(pivotMut.keys())))

#add keys to index and index to keys
for key in set(pivotMut.index):
    if key not in set(pivotMut.keys()):
        pivotMut[key]=[0]*len(pivotMut)
for key in set(pivotMut.keys()):
    if key not in set(pivotMut.index):
        pivotMut.loc[key,:]=[0]*len(pivotMut.keys())
for key in set(pivotWT.index):
    if key not in set(pivotWT.keys()):
        pivotWT[key]=[0]*len(pivotWT)
for key in set(pivotWT.keys()):
    if key not in set(pivotWT.index):
        pivotWT.loc[key,:]=[0]*len(pivotWT.keys())

Mut_keys = set(pivotMut.keys())
WT_keys = set(pivotWT.keys())

extra_keys_in_WT = []
extra_keys_in_Mut = []
for key in pivotWT.keys():
    if key not in Mut_keys:
        extra_keys_in_WT.append(key)
        pivotMut[key]=[0]*len(pivotMut)
        pivotMut.loc[key, :] = [0] * len(pivotMut.keys())
for key in pivotMut.keys():
    if key not in WT_keys:
        extra_keys_in_Mut.append(key)
        print(key)
        pivotWT[key] = [0] * len(pivotWT)
        pivotWT.loc[key, :] = [0] * len(pivotWT.keys())
        print(pivotWT[key])
# for key in extra_keys:
#     pivotMut.loc[key, :] = [0]*len(pivotMut.keys())
print("extra_keys_in_WT", sorted(extra_keys_in_WT))
print("extra_keys_in_Mut", sorted(extra_keys_in_Mut))
assert len(pivotWT)==len(pivotMut)
assert set(pivotMut.keys())==set(pivotWT.keys())
# print(len(set(pivotWT.keys())), len( set(pivotMut.keys())))

sortPivotMut = pivotMut[sorted(pivotMut.keys())]
sortPivotMut.fillna(value=0, inplace=True)
sortPivotMut.sort_index(inplace=True)
pivotWT = pivotWT[sorted(pivotWT.keys())]
pivotWT.sort_index(inplace=True)
pivotWT.fillna(value=0, inplace=True)
#delete columns with deletion from matrices
cols_to_drop = []
# for n,key in enumerate(sortPivotMut.keys()):
#     assert key ==pivotWT.keys()[n]
#     if int_st <= key <=int_en:
#         sortPivotMut[key] = [0] * len(sortPivotMut)
#         sortPivotMut.loc[key, :] = [0] * len(sortPivotMut)
#         pivotWT[key] = [0] * len(sortPivotMut)
#         pivotWT.loc[key, :] = [0] * len(pivotWT)
#         cols_to_drop.append(key)
#     if n!=0:
#         if sortPivotMut.keys()[n]-sortPivotMut.keys()[n-1]!=25000:
#             print(sortPivotMut.keys()[n],sortPivotMut.keys()[n-1])
#             print(pivotWT.keys()[n], pivotWT.keys()[n - 1])
#
#             # print(max(pivotWT.keys()))
#             # print(min(pivotWT.keys()))
#         elif pivotWT.keys()[n] - pivotWT.keys()[n-1]!=25000:
#             print("wt")
#             print(pivotWT.keys([n]), pivotWT.keys()[n-1])
#         assert sortPivotMut.keys()[n]-sortPivotMut.keys()[n-1]==25000
for n,key in enumerate(sortPivotMut.keys()):
    assert key ==pivotWT.keys()[n]
for n,key in enumerate(sortPivotMut.index):
    assert key == pivotWT.index[n]
# print(cols_to_drop)

# drop columns for picture
start_capture =60000000
end_capture = 61200000
drop_columns = []
drop_indeces = []
for n,key in enumerate(sortPivotMut.keys()):
    assert key ==pivotWT.keys()[n]
    if key < start_capture or key > end_capture:
        drop_columns.append(key)
for n,key in enumerate(sortPivotMut.index):
    assert key == pivotWT.index[n]
    if key < start_capture or key > end_capture:
        drop_indeces.append(key)
print("drop",
      len(drop_columns),drop_columns)
# we multiply the matrix corresponding to the mutation (experimental and simulated)
#  by a factor that equalizes the reads count equivalent of the regions that are not involved in the mutation
Mut_del=sortPivotMut.drop(columns=cols_to_drop, index=cols_to_drop)
WT_del=pivotWT.drop(columns=cols_to_drop, index=cols_to_drop)
WT_sum = sum(WT_del.sum())
Mut_sum=sum(Mut_del.sum())
# print(WT_sum, Mut_sum)
print("test", len(pivotWT.keys()), len(sortPivotMut.keys()), len(pivotWT), pivotWT.keys(), min(list(pivotWT.keys())),max(list(pivotWT.keys()))  )
pivotWT.drop(columns=drop_columns, index=drop_indeces, inplace=True)
WT_array = pivotWT.as_matrix()
print(pivotWT)
print(len(WT_array), WT_array)
WT_diags = get_all_digonals(WT_array)
sortPivotMut.drop(columns=drop_columns, index=drop_indeces, inplace=True)
Mut_array = sortPivotMut.as_matrix()
read_coef = WT_sum/Mut_sum
# Mut_array=Mut_array*read_coef
diff_array = Mut_array/WT_array
#To take into account the genomic distance bias, we normalized the difference matrix by dividing each sub-diagonal by the average
# wt reads count at its corresponding pairwise genomic distance.
# for i in range(len(diff_array)):
#     X = np.array(range(0,len(diff_array)-i))
#     Y = np.array(range(i,len(diff_array)))
#     coeff = np.average(WT_array[X,Y])
#     diff_array[X,Y] = diff_array[X,Y] / coeff
#     diff_array[Y,X] = diff_array[X,Y]
# print(diff_array)
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
# np.nan_to_num(diff_diags, copy=False)
# ectopic_array = find_ectopic_interactions(diff_diags)
# print("ectopic")
# print(ectopic_array)
# print(np.min(ectopic_array))
# print(np.max(ectopic_array))

# WT_array=np.log(WT_array)
# print(WT_array)
# plt.title("WT_array")
# plt.imshow(WT_array, cmap="OrRd")
# plt.colorbar()
# plt.show()
# # #
# Mut_array=np.log(Mut_array)
# print(Mut_array)
# plt.imshow(Mut_array, cmap="OrRd")
# plt.colorbar()
# plt.title("Mut_array")
# plt.show()

# diff_array=np.log(diff_array)
#draw line for gene coordinates
#unc5b chr10:60762610-60820641
#Slc29a3 chr10:60720861-60752794
start_unc5b = 60762610
unc5b_bin_left = start_unc5b//25000*25000
unc5b_bin_right = unc5b_bin_left+25000
unc5b_proportion = (start_unc5b-unc5b_bin_left)/25000
unc5b_index = int(np.where(pivotWT.keys()==unc5b_bin_left)[0])
unc5b_line_coord = unc5b_index + unc5b_proportion

end_Slc29a3 = 60720861
Slc29a3_bin_left = 60752794//25000*25000
Slc29a3_bin_right = Slc29a3_bin_left+25000
Slc29a3_proportion = (end_Slc29a3-Slc29a3_bin_left)/25000
Slc29a3_index = int(np.where(pivotWT.keys()==Slc29a3_bin_left)[0])
Slc29a3_line_coord = Slc29a3_index + Slc29a3_proportion
print("slc", Slc29a3_index, pivotWT.keys()[Slc29a3_index], Slc29a3_line_coord)
# print("unc5b", unc5b_index, pivotWT.keys()[54], unc5b_line_coord)
# print(unc5b_bin_left, unc5b_bin_right, unc5b_proportion)
# print("start", start_unc5b)
rows,cols = np.where(diff_array==1)
plt.title("diff_array")
print(diff_array)
print(rows, cols)
print(diff_array[rows, cols])
diff_array[rows, cols]='nan'
plt.imshow(diff_array, cmap="OrRd")
# plt.vlines([unc5b_line_coord, Slc29a3_line_coord], 0, len(diff_array),colors='dodgerblue')
# plt.vlines(Slc29a3   _line_coord, 0, len(diff_array),colors='k')
plt.colorbar()
plt.show()
# plt.clf()
# colored_data = plt.imshow(diff_array, cmap="OrRd")
# with open("out/analysis/rearrangement/ect_array_real.pickle", 'wb') as f:
#     pickle.dump(ectopic_array, f)
# ectopic_array[ectopic_array<=2] = np.NaN
# plt.title("ectopic_interactions")
# plt.imshow(ectopic_array, cmap="OrRd")
# plt.colorbar(colored_data)
# plt.show()
# plt.imshow(ectopic_array)

# plt.imshow(WT_array[0:199], cmap="OrRd")
# plt.scatter([1,1],[2,1])

