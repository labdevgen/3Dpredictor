# Written by Minja, 09-2018
# Class to deal with Hi-C loops
# Stores Hi-C loops as data


from shared import FileReader

#S
class LoopReader(FileReader):
    def read_loops(self):
        pass

    def getLoops(self,chr):
        return {chr:self.loops[chr]}