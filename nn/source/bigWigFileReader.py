import logging
import os
import sys
import pyBigWig

# Add main directory to import path
import_path = os.path.dirname(os.path.dirname(os.getcwd()))
logging.getLogger(__name__).info("Appending import path: "+import_path)
sys.path.append(import_path)

from shared import FileReader

class bigWigReader(FileReader): #Reading, processing and storing the data from bigWig genome tracks
        def __init__(self, fname, name=None):
            self.data = None # this is dict of numpy arrays with values
                             # only present if inMemory = True
            self.dataPointer = None # points to opened bigWig file
            self.inMemory = None
            if name == None:
                logging.getLogger(__name__).warning("Using filename as a name for predictor")
                self.proteinName = os.path.basename(fname)
            else:
                self.proteinName = name
            assert self.proteinName != ""
            super(bigWigReader, self).__init__(fname)

        def readData(self, inMemory = False):
            #inMemory = True - story decompressed bigWig in memory as dict of np arrays
            #inMemory = False - dianmically load it from resource

            self.inMemory = inMemory
            logging.debug(self.fname)
            self.dataPointer = pyBigWig.open(self.fname)
            if inMemory == True:
                self.chrms = self.dataPointer.chroms().keys()
                # returns smthg like dict_proxy({'1': 195471971L, '10': 130694993L})
                for chr in self.chrms.keys():
                    logging.getLogger(__name__).info("Loading chrm " + chr + " in memory")
                    self.chrms[chr] = self.dataPointer.values(chr,0,chrms[chr])

        def get_interval(self, interval):
            if self.inMemory:
                return self.chrms[interval.chr][interval.start:interval.end]
            else:
                return self.dataPointer.values(interval.chr,interval.start,interval.end)