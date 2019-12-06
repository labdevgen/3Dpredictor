import time
import os, errno
import subprocess
import math
import numpy as np
import pandas as pd
import pickle
import matplotlib.pyplot as plt
from matrix_plotter import MatrixPlotter
from termcolor import colored
from shared import  get_bin_size


#   This func is for visualising Predictor results.
#
#   matplot_obj             Matrix_Plotter Object
#   fname                   name of current visualisation
#   out_folder              folder where all the visualisation files are
#
#   ----------------------
#
#   Love, Emil

def MatPlot2HiC(matplot_obj, fname, out_folder):
    def Pandas2ChrSizes(chrsizes_filename,
                        pandas_df):  # This func takes all the chromosomes from pandas object, find out their sizes and write into file
        chromosomes = pandas_df.ix[:, 0].unique()
        chrsizes_table = pd.DataFrame(columns=chromosomes)

        for i in range(len(chromosomes)):
            buf = pandas_df.loc[pandas_df['chr'] == chromosomes[i]][['contact_st', 'contact_en']]
            max1 = buf.max().max()
            chrsizes_table.at[0, chromosomes[i]] = max1

            print('Completed: {}%'.format(i * 100 // len(chromosomes)), end='\r')

        chr_list = list(chrsizes_table)

        chrsizes_file = open(chrsizes_filename, 'w')

        for j in range(len(chr_list)):
            chrsizes_file.write(chr_list[j] + '\t' + str(chrsizes_table.iloc[0][chr_list[j]]) + '\n')

        chrsizes_file.close()

    def Pandas2Pre(pre_filename, pandas_df):  # This func makes pre-HiC file from the pandas object, control or data
        pre_file = open(pre_filename, 'w')
        data_rows = pandas_df.shape[0]

        pandas_df.columns = ["chr1", "start", "end", "count"]
        pandas_df['str1'] = 0
        assert len(pandas_df.loc[(pandas_df['count'] < 0.000001) & (pandas_df['count'] != 0)]) < (len(pandas_df['count']) / 10)
        pandas_df['exp'] = pandas_df['count'] * ( 1000000 )
        pandas_df['exp'] = round(pandas_df['exp']).astype(int)

        pandas_df.to_csv(pre_file, sep=" ",
                         columns=['str1', 'chr1', 'start', 'start', 'str1', 'chr1', 'end', 'end', 'exp'], header=False,
                         index=False)

        pre_file.close()

    # make dirs
    try:
        os.makedirs(out_folder + '/' + fname)
        os.makedirs(out_folder + '/' + fname + '/pre')
        os.makedirs(out_folder + '/' + fname + '/hic')
    except OSError as e:
        if e.errno != errno.EEXIST:
            raise

    # make filenames
    chromsizes_filename = out_folder + '/' + fname + '/pre/chrom.sizes'
    pre_data_filename = out_folder + '/' + fname + '/pre/pre_data.txt'
    hic_data_filename = out_folder + '/' + fname + '/hic/data.hic'
    pre_control_filename = out_folder + '/' + fname + '/pre/pre_control.txt'
    hic_control_filename = out_folder + '/' + fname + '/hic/control.hic'

    # make chrom.sizes, pre-Hic for data and control
    print('Make chromosome sizes file...')
    time1 = time.time()
    Pandas2ChrSizes(chromsizes_filename, matplot_obj.data)
    time2 = time.time()
    print('Time: ' + str(round(time2 - time1, 3)) + ' sec\n')
    print(colored("[SUCCESS]", 'green') + ' Chromosome sizes file created.\n')

    print('Make data pre-HiC file...')
    time1 = time.time()
    Pandas2Pre(pre_data_filename, matplot_obj.data)
    time2 = time.time()
    print('Time: ' + str(round(time2 - time1, 3)) + ' sec\n')
    print(colored("[SUCCESS]", 'green') + ' DATA pre-HiC file created.\n')

    print('Make control pre-HiC file...')
    time1 = time.time()
    Pandas2Pre(pre_control_filename, matplot_obj.control)
    time2 = time.time()
    matplot_obj.columns = ["chr1", "start", "end", "count"]
    binsize = str(get_bin_size(matplot_obj.control, fields=["start", "start"]))
    print('Time: ' + str(round(time2 - time1, 3)) + ' sec\n')
    print(colored("[SUCCESS]", 'green') + ' CONTROL pre-HiC file created.\n')

    #call juicer
    subprocess.call(
        ['java', '-jar', './juicer_tools.jar', 'pre', pre_data_filename, hic_data_filename, chromsizes_filename, '-n',
         '-r', binsize])
    print(colored("[SUCCESS]", 'green') + ' DATA HiC file created.\n')

    subprocess.call(
        ['java', '-jar', './juicer_tools.jar', 'pre', pre_control_filename, hic_control_filename, chromsizes_filename,
         '-n', '-r', binsize])
    print(colored("[SUCCESS]", 'green') + ' CONTROL HiC file created.\n')
