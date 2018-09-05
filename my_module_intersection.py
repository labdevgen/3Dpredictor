import numpy as np
import pandas as pd
def intersect_with_loops(intervals_file,loops_dict,start_field = "start",end_field = "end"):
    intervals_df = pd.read_table(intervals_file)
    assert start_field in intervals_df \
           and end_field in intervals_df \
           and "chr" in intervals_df \
           and "x1" in next(iter(loops_dict.values())) \
           and "x2" in next(iter(loops_dict.values())) \
           and "index" in next(iter(loops_dict.values()))

    # print(promoters_file)
    chr_intervals = dict()  # this dict looks like chr --> dataframe(chr|prom_start|prom_end)
    for chr in pd.unique(intervals_df['chr']):
        chr_intervals[chr] = pd.DataFrame(intervals_df[intervals_df['chr'] == chr])
        assert np.all(
            chr_intervals[chr][start_field][1:] - chr_intervals[chr][start_field][:-1]) >= 0  # check sorting
        #assert np.all(chr_intervals[chr][start_field][1:] - chr_intervals[chr][end_field][:-1]) > 5000
        # print(chr_promoters["chr10"])
    assert sum([len(chr_intervals[p]) for p in chr_intervals.keys()]) == len(intervals_df['chr'])
    result = {}
    statistics = [] #list store some temporary statistics info
    for chr in chr_intervals:
        print(chr)
        if not chr in loops_dict: #Take care about chrms where there are no loops
            result[chr] = pd.DataFrame([])
            print ("WARNING: not loops on chr ",chr)
            continue
        st_end_i = np.searchsorted(loops_dict[chr]['x2'], chr_intervals[chr][start_field])
        end_st_i = np.searchsorted(loops_dict[chr]['x1'], chr_intervals[chr][end_field])

        assert np.all(end_st_i - st_end_i) <= 2 #check that end_st_i always larger than st_end_i
        assert len(st_end_i) == len(end_st_i) == len(chr_intervals[chr][end_field])

        chr_result = [] #resultind df for current chr will be here
        chr_result_loops = [] #loops on this chr will be here
        #now construct list of rows, each row is promoter
        #number of rows for each promoter equal to number of loops this promoter overlap
        for ind,val in enumerate(st_end_i):
            chr_result += [chr_intervals[chr].iloc[ind]]*(end_st_i[ind] - st_end_i[ind])
            chr_result_loops += list([loops_dict[chr].iloc[i]["index"] for i in range(st_end_i[ind],end_st_i[ind])])
            assert end_st_i[ind] - st_end_i[ind] >= 0
        assert len(chr_result) == len(chr_result_loops)
        statistics += list(end_st_i - st_end_i)
        chr_result = pd.DataFrame(chr_result)
        chr_result["loop"] = chr_result_loops
        result[chr] = chr_result
    print ("INFO: number of loops related to promoter destribution:")
    print (np.bincount(statistics))
    print()
    return result
