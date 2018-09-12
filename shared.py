import os,sys
import logging
import inspect
from hashlib import sha224
import numpy as np
import pandas as pd

class Interval:
    def __init__(self,chr,start,end=None,strand=0):
        self.chr = chr
        self.start = start
        if end == None:
            end = start
        self.end = end
        self.len = end - start
        assert self.len >= 0
    def __eq__(self, other):
        return self.chr == other.chr and \
                self.start == other.start and \
                self.end == other.end

    def __gt__(self, other):
        return self.len > other.len

    def __repr__(self):
        return self.__class__.__name__+" "+" ".join(map(str,[self.chr,self.start,self.end]))

    def toFileName(self):
        return self.__class__.__name__ + "_" + "_".join(map(str, [self.chr, self.start, self.end]))

class Position(Interval):
    def __gt__(self, other):
        #if self.chr != other.chr:
        #    raise Exception("Cann't compare positions on different chromosomes")
        if self.chr != other.chr:
            return self.chr > other.chr
        return self.start > other.start

class FileReader(object):
    def __init__(self,fname):
        self.fname = fname
        if not os.path.isfile(fname):
            logging.error("File "+fname+" does not exists")
            sys.exit(11)

def expand_lnk_path(lnk):
    try:
        os.listdir(lnk)
        return lnk + "/"
    except:
        print ("FILE")
        if os.name == "nt":
            import win32com.client
            shell = win32com.client.Dispatch("WScript.Shell")
            shortcut = shell.CreateShortCut(lnk)
            return  shortcut.Targetpath + "/"
        else:
            raise Exception

class Parameters (object):
    def __repr__(self):
        members = inspect.getmembers(self, lambda a: not (inspect.isroutine(a)))
        members = [a[1] for a in members if not (a[0].startswith('__') and a[0].endswith('__')) and \
                   (isinstance(a[1],int) or isinstance(a[1],str))]
        return ".".join(map(str,members))

def str2hash(s,maxlen=100): # Used to hash long file names into shorter ones.
                            # Long fnames (> ~150 characters) are not supported by Windows
                            # This converts long part of the filename "s" which is >100 chars to short hash <= 9 chars
    if len(s) < maxlen:
        return s
    else:
        h = str(int(sha224(s.encode()).hexdigest(), 16) % (10 ** 10))
        logging.warning("Hashing string \n"+s+"\n to string "+h)
        return h

def intersect_intervals(chr_int_data1, chr_int_data2): #input: chr_int_datas are 2 dictionaries of pd.dataframes where 1,2,3 columns == chr, start, end,
                                                       #key of dict: chr output:
                                                       #return chr_int_data2 with additional column 'intersection'. It is column with row indices
                                                       #of chr_int_data1 which intersect intervals in chr_data_2
    if len(chr_int_data1.keys()) != len(chr_int_data2):
        logging.warning("Data have different number of chromosomes", chr_int_data1, '=', len(chr_int_data1.keys()), 'chrs', \
                      chr_int_data2, '=', len(chr_int_data2), 'chrs')
    for chr in chr_int_data2.keys():
        if chr != 'chr1':
            continue
        print(chr)
        if not chr in chr_int_data1:
            logging.warning("Not intervals on chr", chr)
            continue
        st_end_i = np.searchsorted(chr_int_data1[chr]['end'], chr_int_data2[chr]['start'])
        end_st_i = np.searchsorted(chr_int_data1[chr]['start'], chr_int_data2[chr]['end'])
        assert np.all(end_st_i - st_end_i) <= 2  # check that end_st_i always larger than st_end_i
        assert len(st_end_i) == len(end_st_i) == len(chr_int_data2[chr]['end'])
        intersection_result = []
        for ind,val in enumerate(st_end_i):
            if end_st_i[ind] == st_end_i[ind]:
                logging.warning(chr_int_data2[chr].iloc[ind] + 'do not intersect other data')
                intersection_result.append(None)
            elif end_st_i[ind] > st_end_i[ind]:
                indices = []
                [indices.append(i) for i in range(st_end_i[ind], end_st_i[ind])]
                intersection_result.append(indices)
            else:
                logging.error('st_end_i larger then end_st_i')
       # print(len(intersection_result), len(chr_int_data2[chr]['end']))
        assert len(intersection_result) == len(chr_int_data2[chr]['end'])
        chr_int_data2[chr]['intersection'] = intersection_result
        #break
    return chr_int_data2

