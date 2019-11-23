import logging
import os
import sys
import pyBigWig
import numpy as np

# Add main directory to import path
import_path = os.path.dirname(os.path.dirname(os.getcwd()))
logging.getLogger(__name__).info("Appending import path: "+import_path)
sys.path.append(import_path)

from shared import FileReader

class bigWigReader(FileReader): #Reading, processing and storing the data from bigWig genome tracks
        def __init__(self, fname, genome, name=None):
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
            self.genome = genome
            if self.genome == None:
                logging.getLogger(__name__).warning("We HIGLY RECOMEND to provide genome \
                                                    object to avoid genome version conflicts")

            super(bigWigReader, self).__init__(fname)

        def readData(self, inMemory = False, chunkSize = 50000000):
            #inMemory = True - story decompressed bigWig in memory as dict of np arrays
            #inMemory = False - dianmically load it from resource

            self.inMemory = inMemory
            logging.debug(self.fname)
            self.dataPointer = pyBigWig.open(self.fname)
            chrmSizes = self.dataPointer.chroms()

            # check that all genomic chrms are in chrmSize
            # and their sizes are similar
            for chr in self.genome.chrmSizes:
               assert chrmSizes[chr] == self.genome.chrmSizes[chr]

            if inMemory == True:
                self.data = {}
                # returns smthg like dict_proxy({'1': 195471971L, '10': 130694993L})
                for chr in chrmSizes.keys(): # skip un chrms
                    if not chr in self.genome.chrmSizes:
                        logging.warning("Skipping chrm "+chr)
                        continue
                    logging.getLogger(__name__).info("Loading chrm " + chr + " into memory")

                    #allocate array
                    #TODO shell we assert for overflow? How to check it?
                    self.data[chr] = np.zeros(shape=(self.genome.chrmSizes[chr]),
                                              dtype = np.float32)

                    #split data per chunks
                    assert  self.genome.chrmSizes[chr] > 2
                    assert chunkSize > 0
                    if self.genome.chrmSizes[chr] % chunkSize == 0:
                        chunkStatrs = np.arange(0,self.genome.chrmSizes[chr]-1,chunkSize)
                    else:
                        chunkStatrs = np.arange(0,self.genome.chrmSizes[chr],chunkSize)
                    chunkStops = chunkStatrs + chunkSize
                    chunkStops[-1] = self.genome.chrmSizes[chr]
                    assert np.sum(chunkStops - chunkStatrs) == self.genome.chrmSizes[chr]

                    # fill array by chunks
                    for st,en in zip(chunkStatrs,chunkStops):
                        logging.debug("Filling chunk "+str(st))
                        self.data[chr][st:en] = self.dataPointer.values(chr,st,en)

        def get_interval(self, interval):
            if self.inMemory:
                return self.data[interval.chr][interval.start:interval.end]
            else:
                return self.dataPointer.values(interval.chr,interval.start,interval.end)