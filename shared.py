import  os,sys
import logging

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