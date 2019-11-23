import logging
import os
import sys
import gzip
import numpy as np

# Add main directory to import path
import_path = os.path.dirname(os.path.dirname(os.getcwd()))
logging.getLogger(__name__).info("Appending import path: "+import_path)
sys.path.append(import_path)

from shared import str2hash

class fastaReader(object): #Reading, processing and storing the data from
                                # fasta/multifasta genome files
    def add_file(self,fname):
        logging.getLogger(__name__).info("Processing file "+os.path.abspath(fname))
        self.files.append(os.path.abspath(fname))
        if fname.endswith(".gz"):
            handle = gzip.open(fname)
        else:
            handle = open(fname)

        chrm = None
        seq = []
        for line in handle:
            if line[0] == ">":
                chrm = line.strip().split()[0][1:]
                if chrm in self.data.keys():
                    raise Exception("chrm "+chrm+" found twise")
                logging.getLogger(__name__).info(str("Processing chrm "+chrm))

                if chrm != None:
                    self.chrmSizes[chrm] = len(seq)
                    self.data[chrm] = np.array(seq,dtype=np.uint8)
                seq = []
            else:
               if chrm == None:
                   raise Exception("Not correct fasta format")
               for temp in line.strip():
                    seq.append(self.converter[temp])
        assert chrm != None
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
                 converter =  {"A":0,"a":0,
                               "T":1,"t":1,
                               "G":2,"g":2,
                                "C":3,"c":3,
                               "N":4,"n":4},
                 name=None):
        #fpath could be:
        #1.) path to single fasta/multifasta file
        #2.) path to folder with fasta/multifasta file(s) - in this case all files will be used
        #3.) iterable with files/folders
        # converter is used to converts letters to integers

        # initialize variables
        self.chrmSizes = {}
        self.files = []
        self.data = {}
        self.converter = converter

        if type(files)==str: # for single file/dir simply load it
            self.add_fileitem(files)
        else: # for lists iter over list and load each dir/tree
            try:
                for f in files:
                    self.add_fileitem(f)
            except TypeError:
                self.add_fileitem(files)

        assert len(self.files) != 0

        if name == None: # create some meaningful name
            self.name = str2hash("".join(sorted(self.files)))
        else:
            self.name = name

    def get_interval(self, interval):
        return self.data[interval.chr][interval.start:interval.end]

    def __repr__(self):
        s="genome "+self.name+"\n"
        for (chr,size) in self.chrmSizes.items():
            s+=chr+"\t"+str(size)+"\n"
        s+="Data loader from files: "
        s+="\n".join(self.files)
        return s