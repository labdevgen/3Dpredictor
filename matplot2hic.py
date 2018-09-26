import os, errno
import subprocess
import math
import numpy as np
import pandas as pd
import pickle
import matplotlib.pyplot as plt
from matrix_plotter import MatrixPlotter
from termcolor import colored

# This func is for visualising Predictor results.
#
#   matplot_obj             Matrix_Plotter Object
#   fname                   name of current visualisation
#   out_folder              folder where all the visualisation files are
#   juicer_tools_folder     folder where Juicer is 

def MatPlot2HiC(matplot_obj, fname, out_folder, juicer_tools_folder):
    def Pandas2ChrSizes(chrsizes_filename, pandas_df): # This func takes all the chromosomes from pandas object, find out their sizes and write into file
        data_rows = pandas_df.shape[0]
        chrsizes_table = pd.DataFrame()
        
        for i in range(data_rows):
            chr_buf = None
            contact_st_buf = None
            contact_en_buf = None
            
            chr_buf = pandas_df.iloc[i]['chr']
            contact_st_buf = pandas_df.iloc[i]['contact_st']
            contact_en_buf = pandas_df.iloc[i]['contact_en']
            
            if chr_buf not in chrsizes_table:
                chrsizes_table[chr_buf] = [0]
            
            if chrsizes_table.iloc[0][chr_buf] < contact_st_buf:
                chrsizes_table.at[0, chr_buf] = contact_st_buf
        
            if chrsizes_table.iloc[0][chr_buf] < contact_en_buf:
                chrsizes_table.at[0, chr_buf] = contact_en_buf
        
            print('Completed: {}%'.format(i * 100 // data_rows), end='\r')
        
        chr_list = list(chrsizes_table)
        
        chrsizes_file = open(chrsizes_filename, 'w')
        
        for j in range(len(chr_list)):
            chrsizes_file.write(chr_list[j] + '\t' + str(chrsizes_table.iloc[0][chr_list[j]]) + '\n')
        
        chrsizes_file.close()
            
    def Pandas2Pre(pre_filename, pandas_df): # This func makes pre-HiC file from the pandas object, control or data
        data_rows = pandas_df.shape[0]
        pre_file = open(pre_filename, 'w')
        
        for i in range(data_rows):
            chr_buf = None
            contact_st_buf = None
            contact_en_buf = None
            contact_count_exp_buf = None
    
            chr_buf = pandas_df.iloc[i]['chr']
            contact_st_buf = pandas_df.iloc[i]['contact_st']
            contact_en_buf = pandas_df.iloc[i]['contact_en']
            contact_count_exp_buf = round(math.exp(pandas_df.iloc[i]['contact_count']))
        
            for k in range(contact_count_exp_buf):
                pre_file.write('0 ' + chr_buf + ' ' + str(contact_st_buf) + ' ' + str(contact_st_buf) + ' 0 ' + chr_buf + ' ' + str(contact_en_buf) + ' ' + str(contact_en_buf) + '\n')
        
            print('Completed: {}%'.format(i * 100 // data_rows), end='\r')
    
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
    juicer_tools_filename = juicer_tools_folder + '/juicer_tools.jar'
    
    # make chrom.sizes, pre-Hic for data and control
    print('Make chromosome sizes file...')
    Pandas2ChrSizes(chromsizes_filename, matplot_obj.data)
    print(colored("[SUCCESS]", 'green') + ' Chromosome sizes file created.\n')
    
    print('Make data pre-HiC file...')
    Pandas2Pre(pre_data_filename, matplot_obj.data)
    print(colored("[SUCCESS]", 'green') + ' DATA pre-HiC file created.\n')
    
    print('Make control pre-HiC file...')
    Pandas2Pre(pre_control_filename, matplot_obj.control)
    print(colored("[SUCCESS]", 'green') + ' CONTROL pre-HiC file created.\n')

    # call the Juicer
    subprocess.call(['java', '-jar', juicer_tools_filename, 'pre', pre_data_filename, hic_data_filename, chromsizes_filename])
    print(colored("[SUCCESS]", 'green') + ' DATA HiC file created.\n')

    subprocess.call(['java', '-jar', juicer_tools_filename, 'pre', pre_control_filename, hic_control_filename, chromsizes_filename])    
    print(colored("[SUCCESS]", 'green') + ' CONTROL HiC file created.\n')
