# Code by Minja Fishman
# Nov 2019, ICG SB RAS

import inspect
import logging
import os
import sys
import gzip
from collections import OrderedDict

import numpy as np

# Add main directory to import path
import sys, os

# Add main source directory to import path
root_dir = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
source_dir = os.path.join(root_dir,"source")
sys.path.append(source_dir)

from shared import FileReader, get_version

def default_chr_name_ranamer(chr):
    return chr

def rm_chr_from_chrName(chr):
    if chr.startswith("chr"):
        return chr[3:]
    else:
        return chr


class fastaReader(FileReader): #Reading, processing and storing the data from
                                # fasta/multifasta genome files
    def add_file(self,fname):
        self.files.append(fname)

    def read_file(self,fname):
        def add_chrm(chrm, seq): # add seq to data and chrmSizes dictionaries
            assert len(seq) > 0
            assert not chrm in self.chrmSizes.keys()
            assert not chrm in self.data.keys()
            self.chrmSizes[chrm] = len(seq)
            self.data[chrm] = np.array(seq, dtype=np.uint8)

        logging.getLogger(__name__).info("Processing file "+os.path.abspath(fname))
        if fname.endswith(".gz"):
            handle = gzip.open(fname)
        else:
            handle = open(fname)

        chrm = None
        seq = []
        exclude = False

        for line in handle:
            if line[0] == ">": # new chrm starts here
                if chrm != None and not exclude:
                    add_chrm(chrm,seq)
                elif exclude:
                    logging.getLogger(__name__).info("Skipping chrm "+chrm)

                chrm = line.strip().split()[0][1:]
                chrm = self.chrm_names_renamer(chrm)
                assert len(chrm) > 0
                logging.getLogger(__name__).info(str("Found chrm "+chrm))

                if chrm in self.data.keys():
                    raise Exception("chrm "+chrm+" found twice")

                seq = []
                if len(self.useOnlyChromosomes) != 0:
                    exclude = not (chrm in self.useOnlyChromosomes)
                elif len(self.excludeChr) !=0:
                    exclude = (chrm in self.excludeChr)
                if len(self.chrmSizes) == len(self.useOnlyChromosomes) and \
                        len(self.useOnlyChromosomes)>0: # we have read all chromosomes required
                    seq = []
                    break
            else:
               if exclude: continue
               if chrm == None:
                   raise Exception("Not correct fasta format")
               for temp in line.strip():
                    seq.append(self.converter[temp])
        assert chrm != None

        if seq != []:
            add_chrm(chrm,seq)
        else:
            assert exclude # if we found empty seq it should mean last chrm was in exclude list


        handle.close()


    def add_folder(self,fname):
        files = [f for f in os.listdir(fname) if f.endswith(".fa") or \
                                                 f.endswith(".fasta") or \
                                                 f.endswith(".fa.gz") or \
                                                 f.endswith(".fasta.gz")]
        for f in files:
            self.add_file(os.path.abspath(f))

    def add_fileitem(self,f):
        f = os.path.abspath(f)
        if os.path.isdir(f):
            self.add_folder(f)
        elif os.path.isfile(f):
            self.add_file(f)
        else:
            raise Exception("Item " + f + " is neither file nor folder")

    def __init__(self,files,
                 converter =  None,
                 excludeChromosomes = None,
                 useOnlyChromosomes = None,
                 chrm_names_renamer = default_chr_name_ranamer,
                 name=None):
        #fpath could be:
        #1.) path to single fasta/multifasta file
        #2.) path to folder with fasta/multifasta file(s) - in this case all files will be used
        #3.) iterable with files/folders
        # converter is used to converts letters to integers
        # One can provide either excludeChromsomes or useOnlyChromosomes list. Names are self-explanatory =(

        # set default inputs
        if excludeChromosomes is None:
            excludeChromosomes = []
        if useOnlyChromosomes is None:
            useOnlyChromosomes = []
        if converter is None:
            converter = OrderedDict({"A":0,"a":0,
                               "T":1,"t":1,
                               "G":2,"g":2,
                               "C":3,"c":3,
                               "N":4,"n":4})

        # initialize variables
        self.chrmSizes = OrderedDict({})
        self.files = []
        self.data = {}
        self.converter = converter
        self.chrm_names_renamer = chrm_names_renamer
        self.chrm_names_renamer_function = chrm_names_renamer.__name__

        if len(excludeChromosomes)*len(useOnlyChromosomes) != 0:
            logging.getLogger(__name__).error(
                "You can only provide one of these list: useOnlyChromosomes/excludeChromosomes")
            raise Exception("Wrong input")
        self.excludeChr = excludeChromosomes
        self.useOnlyChromosomes = useOnlyChromosomes

        if type(files)==str: # for single file/dir simply load it
            self.add_fileitem(files)
            self.fname = files
        else: # for lists iter over list and load each dir/tree
            try:
                for f in files:
                    self.add_fileitem(f)
            except TypeError:
                self.add_fileitem(files)

        assert len(self.files) != 0

        if name == None: # create some meaningful name
            if len(self.files) == 1:
                self.name = os.path.basename(self.files[0])
            else:
                self.name = "".join(sorted(self.files))
        else:
            self.name = name

        self.full_name = "".join(sorted([os.path.basename(f) for f in self.files])+sorted(self.excludeChr)+sorted(self.useOnlyChromosomes)+\
                                 [str(sorted(self.converter))]+[self.chrm_names_renamer.__name__]+[str(get_version())])

    def read_data(self):
        if os.path.exists(self.get_dump_path()):
            return self.load()
        for i in self.files:
            self.read_file(i)

        XML = self.toXMLDict(exludedMembers=("data"))
        header = self.name
        self.dump(descriptiveXML={"XML":XML,"header":header})
        return self

    def get_interval(self, interval):
        return self.data[interval.chr][interval.start:interval.end]

    def get_chr_sizes(self):
        return self.chrmSizes

    def __repr__(self):
        XMLrepresentation  = self.toXMLDict(exludedMembers=("data"))
        return "\n".join([key+"\t"+str(val) for key,val in XMLrepresentation.items()])