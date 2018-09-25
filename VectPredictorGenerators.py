import pandas as pd
import numpy as np
from PredictorGenerators import PredictorGenerator
from shared import intersect_intervals

class VectPredictorGenerator(PredictorGenerator):
    def __init__(self,**kwargs):
        super(VectPredictorGenerator, self).__init__(**kwargs)
        self.vectorizable = True

class loopsPredictorGenerator(VectPredictorGenerator):
    def __init__(self,loopsReader,window_size):
        super(VectPredictorGenerator, self).__init__(loopsReader = loopsReader,
                                                    window_size = window_size)

    def get_header(self,contact):
        return ["IsLoop"]

    def get_predictors(self,contacts):
        result = pd.DataFrame({"IsLoop":[-1]*len(contacts)})
        left = {}
        right = {}
        for chr in np.unique(contacts["chr"].values):
            idxs = contacts["chr"] == chr
            left[chr] = pd.DataFrame({"start":contacts[idxs]["contact_st"] - self.window_size,
                                      "end": contacts[idxs]["contact_st"] + self.window_size})
            right[chr] = pd.DataFrame({"start":contacts[idxs]["contact_en"] - self.window_size,
                                      "end": contacts[idxs]["contact_en"] + self.window_size})

            left_loops = self.loopsReader.getLeftLoopAncors(chr)
            right_loops = self.loopsReader.getRightLoopAncors(chr)

            intersections_L = intersect_intervals(left_loops, left)["intersection"]
            intersections_R = intersect_intervals(right_loops, right)["intersection"]
            intersectios = pd.concat([intersections_L,intersections_R],axis=1)
            intersectios.columns = ["L","R"]
            predictor = intersectios.apply(lambda x: len(np.intersect1d(x.L,x.R)==0))
            contacts[idxs] = predictor.values
        assert len(result.query("IsLoop == -1"))==0
        return result