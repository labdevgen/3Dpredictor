# Code by Minja Fishman
# Nov 2019, ICG SB RAS

# The purpose of this module is to split genome onto intervals for training and validation
# For the very simple type of split, one could split by chromosomes
# Latter one can suggest more complex splitting strategies.

import logging

class Splitter():
    def __init__(self):
        pass

    def split(self,
              minsize = 1000,
              maxsize = 1000000,
              ):
        # minimum size of the interval
        # maximal size of the interval