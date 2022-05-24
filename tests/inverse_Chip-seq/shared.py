import os,sys,numbers
import logging
import inspect
import requests
from hashlib import sha224
import numpy as np
import pandas as pd
import dicttoxml
from xml.dom.minidom import parseString
from collections import OrderedDict
from functools import partial
import pickle

# changing this values will result in rewriting all dumps of FileReader objects generated previously
# the only exception is for dumps of objects which overload standard toXMLDict function
# or do not use toXMLDict do generate full_name

version = 1

def get_version():
    return version

class Interval:
    def __init__(self,chr,start,end=None,strand=0):
        self.chr = chr
        assert int(start) == start
        self.start = int(start)
        if end == None:
            end = int(start)
        assert int(end) == end
        self.end = int(end)
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

class Genome(object):
    def __init__(self):
        pass

class FileReader(object):
    def check_local_path_exist(self,path):
        return os.path.isfile(path)
    def check_url_path_exist(self,path):
        try:
            r = requests.head(path)
            if not r.status_code in [200,302,304,307,308]:
                return False
            else:
                return True
        except:
            return False

    def __init__(self,fname):
        self.fname = fname
        if not self.check_local_path_exist(fname) and not self.check_url_path_exist(fname):
            logging.error("File "+fname+" does not exists")
            sys.exit(11)

    def toXMLDict(self, exludedMembers = ()):
        # how it works:
        # first, we get all class members, excluding only routine function, such as __eq__, __str__ and etc
        members = inspect.getmembers(self, lambda a: not (inspect.isroutine(a)))

        # next, we filter these items exclude those provided in excludeMembers
        members = OrderedDict((a[0],str(a[1])) for a in sorted(members) if not (a[0].startswith('__') and a[0].endswith('__')) and \
             not (a[0] in exludedMembers))

        # we have global parameter version for the whole code
        # this term is reverved and should not be used by classes
        assert not "version" in members.keys()

        members["version"] = version

        # this gives a dict key--val of all members of class
        return members

    def get_full_name(self):
        # full name is id, which should represent file reader with all its paramters
        # if any of parameters changed, is should change full name
        try:
            return self.full_name
        except:
            raise Exception("Dump is only possible if full name is defined")

    def get_dump_path(self):
        # construct dump path based on fname and full_name
        # full name is id, which should represent file reader with all its parameters

        if os.path.isfile(self.fname):
            self.dump_path = self.fname + "." + str2hash(self.get_full_name())
            return self.dump_path
        elif os.path.isdir(self.fname):
            self.dump_path = os.path.join(self.fname, str2hash(self.get_full_name()))
            return self.dump_path
        else:
            logging.getLogger(__name__).error("Cann't create dump: reader fileName is neither file nor dir: " + self.fname)
            raise Exception()

    def dump(self, descriptiveXML = None):
        dump_path = self.get_dump_path()
        if os.path.exists(dump_path):
            logging.warning("Overwriting file "+dump_path)
        logging.getLogger(__name__).debug("Saving dump to the file "+dump_path)
        pickle.dump(self,open(dump_path,"wb"))
        if descriptiveXML != None:
            XML = descriptiveXML["XML"]
            header = descriptiveXML["header"]
            write_XML(XML,header,dump_path+".xml")

    def load(self):
        logging.getLogger(__name__).debug("Loading dump from the file"+self.get_dump_path())
        return pickle.load(open(self.get_dump_path(),"rb"))

    def read_file(self,
                  renamer = {"0":"chr","1":"start","2":"end","6":"sigVal"}): # example for chip-seq file
        logging.getLogger(__name__).info(msg="Reading ChipSeq file "+self.fname)

        # set random temporary labels
        if self.fname.endswith(".gz"):  # check gzipped files
            import gzip
            temp_file = gzip.open(self.fname)
        else:
            temp_file = open(self.fname)
        Nfields = len(temp_file.readline().strip().split())
        temp_file.close()

        names = list(map(str, list(range(Nfields))))
        data = pd.read_csv(self.fname, sep="\t", header=None, names=names, comment='#')

        # subset and rename
        data_fields = list(map(int,renamer.keys()))
        data = data.iloc[:,data_fields]
        data.rename(columns=renamer,
                        inplace=True)
        #save
        self.data = data
        del data

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

    def toXMLDict(self):
        members = inspect.getmembers(self, lambda a: not (inspect.isroutine(a)))
        members = OrderedDict(a for a in members if not (a[0].startswith('__') and a[0].endswith('__')) and \
                   (isinstance(a[1], numbers.Number) or isinstance(a[1],str)))
        return members


def str2hash(s,maxlen=100): # Used to hash long file names into shorter ones.
                            # Long fnames (> ~150 characters) are not supported by Windows
                            # This converts long part of the filename "s" which is >100 chars to short hash <= 9 chars
    if len(s) < maxlen:
        return s
    else:
        h = str(int(sha224(s.encode()).hexdigest(), 16) % (10 ** 10))
        #logging.getLogger(__name__).debug("Hashing string \n"+s+"\n to string "+h)
        return h

def intersect_intervals(chr_int_data1, chr_int_data2, suppreseChrNumberCheck=False): #input: chr_int_datas are 2 dictionaries of pd.dataframes where 1,2,3 columns == chr, start, end,
                                                       #key of dict: chr
                                                       #output:
                                                       #returns chr_int_data2 with additional column 'intersection'. It is column with row indices
                                                       #of chr_int_data1 which intersect intervals in chr_data_2
    if (len(chr_int_data1.keys()) != len(chr_int_data2.keys())) and (not suppreseChrNumberCheck):
        logging.getLogger(__name__).warning("Data have different number of chromosomes: ")
        logging.getLogger(__name__).warning(str(chr_int_data1.keys()))
        logging.getLogger(__name__).warning(str(chr_int_data2.keys()))
    result = {}
    for chr in chr_int_data2.keys():
        result[chr] = pd.DataFrame([])
        if not chr in chr_int_data1:
            logging.getLogger(__name__).warning("Not intervals on chr", chr)
            result[chr] = pd.DataFrame(columns=chr_int_data2[chr].columns.values)
            continue
        st_end_i = np.searchsorted(chr_int_data1[chr]['end'], chr_int_data2[chr]['start'])
        end_st_i = np.searchsorted(chr_int_data1[chr]['start'], chr_int_data2[chr]['end'])

        # TODO what does this assert mean?
        # In fact, it checks that end_st == st_end occurs not more than 2 times...
        assert np.all(end_st_i - st_end_i) <= 2  # check that end_st_i always larger than st_end_i
        assert len(st_end_i) == len(end_st_i) == len(chr_int_data2[chr]['end'])
        intersection_result = []
        ids_column = []
        chr_intervals_result = []
        for ind,val in enumerate(st_end_i):
            if end_st_i[ind] == st_end_i[ind]:
                pass #it's common situtation, so no need for warning
                #logging.geotLogger(__name__).warning("do not intersect other data " + str(chr_int_data2[chr].iloc[ind]) + '      ' +  str(ind))
            elif end_st_i[ind] > st_end_i[ind]:
                chr_intervals_result += [chr_int_data2[chr].iloc[ind]] * (end_st_i[ind] - st_end_i[ind])
                [ids_column.append(ind) for index in range(st_end_i[ind], end_st_i[ind])]
                [intersection_result.append(index) for index in range(st_end_i[ind], end_st_i[ind])]
            else:
                logging.getLogger(__name__).error('st_end_i larger then end_st_i')
                logging.getLogger(__name__).error(str(st_end_i)+" "+str(end_st_i))
                #As it's an error, I assume raising exeption.
                raise Exception("Exception from intersect_intervals function")
        #print(len(intersection_result), len(chr_intervals_result))
        assert len(intersection_result) == len(chr_intervals_result)
        if len(chr_intervals_result)==0:
            chr_intervals_result = pd.DataFrame(columns=chr_int_data2[chr].columns.values)
        else:
            chr_intervals_result = pd.DataFrame(chr_intervals_result)
        chr_intervals_result["intersection"] = intersection_result
        chr_intervals_result["ids_column"] = ids_column
        result[chr] = chr_intervals_result
    return result

def intersect_with_interval(chr_int_data1,
                               interval,
                               return_ids=False):

    chr = interval.chr
    if not chr in chr_int_data1:
        logging.getLogger(__name__).warning("No intervals on chr", chr)
        return pd.DataFrame({}) #Return empty DataFrame
    intersection = chr_int_data1[chr].index.overlaps(pd.Interval(left=interval.start,
                                                                 right=interval.end,
                                                                 closed="both"))
    if return_ids:
        if len(np.where(intersection)[0]) == 0:
            end = chr_int_data1[chr]["start"].searchsorted(interval.start)
            start = -1
            return start, end-1
        else:
            start = np.min(np.where(intersection)[0])
            end = np.max(np.where(intersection)[0])
            return start,end
    else:
        return chr_int_data1[chr][intersection]

# File descriptions are saved in XML form
# Description should be dict-like
def write_XML(XML_report, header,fname="files_description.xml"):
    dicttoxml.LOG.setLevel(logging.WARNING) #Get rid of tonnes of INFO/DEBUG logs from dicttoxml

    #Add top level object representing file name
    XML_report = {header:XML_report}
    #get xml line without data-types
    xml_line = dicttoxml.dicttoxml(XML_report,attr_type=False)
    #convert to pretty string
    to_write = parseString(xml_line).toprettyxml()

    #write to file
    # if not os.path.exists(fname):
    f = open(fname,"w")
    # else:
    #     f=open(fname, "a")
    f.write(to_write)
    f.close()

def oe2obs(contacts, data, expected_folder, cell_type,coeff_fname, **kwargs): # dists is array, element[i] --> distance between ancors of contact i
    # read expected file
    # First number in this file is for diagonal elements, i.e. where distance = 0
    input_data = data
    expected_file = expected_folder + input_data.iloc[0, input_data.columns.get_loc("chr")] + "." + cell_type + ".expected.txt"
    expected = np.loadtxt(expected_file)
    expected = np.nan_to_num(expected)
    dists = np.array(input_data["contact_dist"])
    #get binsize
    dist = pd.unique(input_data["contact_en"] - input_data["contact_st"])
    sorted_starts = np.sort(input_data["contact_st"].values[:min(1000, len(input_data))])
    dist2 = np.unique(np.subtract(sorted_starts[1:], sorted_starts[:-1]))
    assert (dist2 >= 0).all()
    dist = np.unique(np.concatenate((dist, dist2)))
    dist = dist[np.nonzero(dist)]
    assert len(dist) > 0
    binsize = min(dist)

    expected_dist = dict([(ind*binsize,val) for ind,val in enumerate(expected)]) # dictionary, distance --> expected
    coeff_data = pd.read_csv(coeff_fname, delimiter="\t") #coeff for normalization, we do it because for contacts prediction we divide all contacts on coeff
    coeff = float(coeff_data["coeff"])
    result = []
    for ind,val in enumerate(contacts):
        result.append(val*expected_dist[dists[ind]]/coeff)
    assert len(result) == len(contacts)
    if kwargs['data_type']=='predicted':
        return result
    elif kwargs['data_type']=='validation':
        input_data['contact_count']=result
        return input_data

def decorate_oe2obs(func,expected_folder, cell_type, coeff_fname):
    result = partial(func,expected_folder=expected_folder, cell_type=cell_type, coeff_fname=coeff_fname)
    result.__name__ = str(cell_type) + func.__name__
    return result

def return_coordinates_after_deletion(contacts, data, intervals, **kwargs):
    if kwargs['data_type']=='predicted':
        return contacts
    elif kwargs['data_type'] == 'validation':
        if intervals.__repr__()=='no_deletion':
            return data
        else:
            for interval in intervals:
                data["contact_st"] = data["contact_st"].apply(lambda x: x if x < interval.start else x + interval.len)
                data["contact_en"] = data["contact_en"].apply(lambda x: x if x < interval.start else x + interval.len)
                return data


# TODO check this function!
def return_coordinates_after_duplication(contacts, data, intervals, **kwargs):
    data["contact_count"] = contacts
    for interval in intervals:
        dup_ids = data.query("@interval.start+interval.len < contact_st < @interval.end+interval.len | "
                             + "@interval.start+interval.len < contact_en < @interval.end+interval.len").index
        dup_second_copy_ids = data.query("@interval.start+interval.len < contact_st < @interval.end+interval.len | "
                                         + "@interval.start+interval.len < contact_en+interval.len < @interval.end+interval.len").index
        assert len(dup_ids) == len(dup_second_copy_ids)
        data.loc[dup_ids, "contact_count"] += data.loc[dup_second_copy_ids, "contact_count"]
        data.drop(dup_second_copy_ids, inplace=True)
        data["contact_st"] = data["contact_st"].apply(lambda x: x if x < interval.start else x - interval.len)
        data["contact_en"] = data["contact_en"].apply(lambda x: x if x < interval.start else x - interval.len)
    if kwargs['data_type'] == 'predicted':
        return data["contact_count"]
    elif kwargs['data_type'] == 'validation':
        return data

def decorate_return_coordinates_after_rearrangements(func, intervals):
    result = partial(func, intervals=intervals)
    result.__name__ = intervals.__repr__()
    return result

def get_bin_size(data, fields = ["contact_en","contact_st"]):
    f1 = fields[0]
    f2 = fields[1]
    dist = pd.unique(data[f1] - data[f2])
    sorted_starts = np.sort(data[f1].values[:min(1000, len(data))])
    dist2 = np.unique(np.subtract(sorted_starts[1:], sorted_starts[:-1]))
    assert (dist2 >= 0).all()
    dist = np.unique(np.concatenate((dist, dist2)))
    dist = dist[np.nonzero(dist)]
    assert len(dist) > 0
    binsize = min(dist)

    assert int(binsize) == binsize
    return int(binsize)

def sparse2dense(data,fields=["st","end","oe"], default_value=0, debug_mode = False):
    # Sparse to dense
    # IMPORTANT: 1.) data should be already binned
    # 2.) min(data["st"]) will be subtracted from all coordinates
    # 3.) all missing data will be filled with "default_value"
    # assume 'data' is pd.DataFrame
    # describing matrix in sparce format
    # where columns with names
    # fields[0] and fields[1] are row and column ids of matrix item
    # and column fields[2] is matrix value
    # sparse2dense will return numpy array corresponding to dense matrix
    f1 = fields[0]
    f2 = fields[1]
    f3 = fields[2]
    L = data[f2].max() - data[f1].min() + 1
    array = np.zeros(shape=(L, L)) + default_value
    X = (data[f1] - data[f1].min()).values
    Y = (data[f2] - data[f1].min()).values
    array[X, Y] = data[f3]
    array[Y, X] = data[f3]

    if debug_mode:
    # slow; use ones for debuging only
        global_min = data[f1].min()
        def check_func(series):
            x = int(series[f1] - global_min)
            y = int(series[f2] - global_min)
            return array[x,y] == array[y,x] == series[f3]

        data["check"] = data.apply(check_func,axis = "columns")
        assert np.all(data["check"].values)
        print ("passed check!")
    return array

def makedirs(dir_path):
    try:
        os.makedirs(dir_path)
    except OSError as exc:  # Python >2.5
        if os.path.isdir(dir_path):
            pass
        else:
            logging.error("Error "+str(exc.errno)+": cannot create directory "+dir_path)
            logging.error(str(exc))
            raise
