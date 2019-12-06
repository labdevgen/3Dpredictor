from fastaFileReader import fastaReader
from shared import Interval
import sys
import os
sys.path.append(os.path.dirname(os.path.dirname(os.getcwd())))
from PredictorGenerators import PredictorGenerator

class SequencePredictorGenerator(PredictorGenerator):
    def __init__(self, fastaReader=fastaReader, binsize=0, dist_from_anchor=0, **kwargs):
        self.name = "SequencePredictor"
        self.fastaReader = fastaReader
        self.dist_from_anchor = dist_from_anchor
        self.binsize = binsize
        self.dist_of_interval = self.binsize + 2*self.dist_from_anchor
        super(SequencePredictorGenerator, self).__init__(**kwargs)

    def get_header(self,contact):
        self.header = ["left" + str(i) for i in range(1, self.dist_of_interval+1)] + \
        ["right" + str(i) for i in range(1, self.dist_of_interval+1)]
        example_predictor = self.get_predictors(contact)
        assert len(example_predictor) == len(self.header)
        return self.header

    def get_predictors(self,contact):
        left_interval = Interval(contact.chr, contact.contact_st - self.dist_from_anchor, \
                                 contact.contact_st + self.binsize + self.dist_from_anchor)
        right_interval = Interval(contact.chr, contact.contact_en - self.dist_from_anchor, \
                                 contact.contact_en + self.binsize + self.dist_from_anchor)
        left_sequence = self.fastaReader.get_interval(left_interval)
        right_sequence = self.fastaReader.get_interval(right_interval)
        assert len(left_sequence) == self.dist_of_interval
        assert len(right_sequence) == self.dist_of_interval
        return left_sequence.tolist() + right_sequence.tolist()