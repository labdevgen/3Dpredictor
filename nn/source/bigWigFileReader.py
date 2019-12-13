# Code by Minja Fishman
# Nov 2019, ICG SB RAS

import logging
import os
import sys
import pyBigWig
import numpy as np

# Add main source directory to import path
root_dir = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
source_dir = os.path.join(root_dir,"source")
sys.path.append(source_dir)

from shared import FileReader

class bigWigReader(FileReader): #Reading, processing and storing the data from bigWig genome tracks
        def __init__(self, fname, genome, name=None, inMemory = False):
            #inMemory = True - story decompressed bigWig in memory as dict of np arrays
            #inMemory = False - dianmically load it from resource

            self.data = None # this is dict of numpy arrays with values
                             # only present if inMemory = True
            self.dataPointer = None # points to opened bigWig file
            self.inMemory = inMemory
            if name == None:
                logging.getLogger(__name__).warning("Using filename as a name for predictor")
                self.name = os.path.basename(fname)
            else:
                self.name = name

            assert self.name != ""
            self.genome = genome
            self.full_name = str(self.toXMLDict(exludedMembers=("data","dataPointer")))
            super(bigWigReader, self).__init__(fname)

        def dump(self):
            # we don't want to dump large genome object
            # and pickl cann't dump datapointers

            descriptiveXML = {"XML": self.toXMLDict(exludedMembers=("data","dataPointer")),
                              "header": self.name}
            temp1 = self.genome
            try:
                temp2 = self.dataPointer
                del self.dataPointer
            except:
                temp2 = None

            super(bigWigReader, self).dump(descriptiveXML=descriptiveXML)
            self.genome = temp1
            self.dataPointer = temp2


        def load(self, genome):
            result = super(bigWigReader, self).load()
            result.genome = genome # we didn't dump large genome object
            result.dataPointer = pyBigWig.open(result.fname) # pickl cann't dump datapointers, so we didn't dump it
            return result

        def readData(self, chunkSize = 50000000, noDump = False, debugMode = False):
            # if noDump, will not try load data from dump file
            #  debug_mode = False will skip some assertions, useful only for intentionally incorrect input data

            # first try to load data from dump
            if os.path.isfile(self.get_dump_path()) and not noDump and self.inMemory:
                return self.load(self.genome)

            logging.debug(self.fname)
            self.dataPointer = pyBigWig.open(self.fname)
            chrmSizes = self.dataPointer.chroms()

            # check that all genomic chrms are in chrmSize
            # and their sizes are similar
            for chr in self.genome.chrmSizes:
                if not debugMode:
                    assert chrmSizes[chr] == self.genome.chrmSizes[chr]

            if self.inMemory == True:
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
                        # TODO this data contains NaN, and I don't transform it
                        # so be aware of NaNs!
                        #assert np.all(np.isfinite(self.data[chr][st:en]))
                self.dump()
            return self

        def get_interval(self, interval):
            if self.inMemory:
                return self.data[interval.chr][interval.start:interval.end]
            else:
                return self.dataPointer.values(interval.chr,interval.start,interval.end)