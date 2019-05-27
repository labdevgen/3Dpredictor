import logging
import pandas as pd
import numpy as np
from PredictorGenerators_new_edition import PredictorGenerator
from shared import intersect_intervals

class VectPredictorGenerator(PredictorGenerator):
    def __init__(self,**kwargs):
        super(VectPredictorGenerator, self).__init__(**kwargs)
        self.vectorizable = True

class loopsPredictorGenerator(VectPredictorGenerator):
    def __init__(self,loopsReader,window_size):
        self.name = loopsReader.name
        super(VectPredictorGenerator, self).__init__(loopsReader = loopsReader,
                                                    window_size = window_size)
        self.vectorizable = True

    def get_header(self,contact):
        return ["IsLoop"]

    def get_predictors(self,contacts):
        # Basic idea:
        # First, set predictor for all contacts to 0
        # Next, intersect left and right anchor with loops
        # Then, intersect these intersections to get contacts with both anchors belonging to same loop
        # Finally set predictors for these contacts to 1

        result = pd.DataFrame({"IsLoop": [0]*len(contacts)})
        left = {}
        right = {}
        for chr in np.unique(contacts["chr"].values):

            # Important think about index
            # Contacts are not assumned to belong to same chr or be sorted by chr
            # But other funcs, i.e. intersect_intervals operate on chr-based dicts of intervals
            # To solve this, we have idxs which is boolean index of contacts belonging to single chrm
            # And we will get intersections["ids_column"] which is idxs of those elements of idxs,
            # which have intersections. I.e. if we have chr1 in 3rd, 5th and 6th contact, and 6th
            # has intersection with loop, intersections["ids_column"] will be not 6, but 2.
            # To remap from intersections["ids_column"] to initial indexing of contacts we
            # use np.flatnonzero(idxs)[idxs2] statment (see blow

            idxs = contacts["chr"] == chr
            left[chr] = pd.DataFrame({"start":contacts[idxs]["contact_st"] - self.window_size,
                                      "end": contacts[idxs]["contact_st"] + self.window_size})
            right[chr] = pd.DataFrame({"start":contacts[idxs]["contact_en"] - self.window_size,
                                      "end": contacts[idxs]["contact_en"] + self.window_size})

            left_loops = self.loopsReader.getLeftLoopAncors(chr)
            if len(left_loops[chr]) == 0: # No loops on this chr
                logging.getLogger(__name__).warning("No loops on chr "+chr)
                continue
            right_loops = self.loopsReader.getRightLoopAncors(chr)
            if len(right_loops[chr]) == 0:
                continue

            #print (left)
            #print (left_loops)
            intersections_L = intersect_intervals(left_loops, left)[chr]
            intersections_L["Loop_id"] = intersections_L.intersection.apply(lambda x: left_loops[chr].id.iloc[x])

            intersections_R = intersect_intervals(right_loops, right)[chr]
            intersections_R["Loop_id"] = intersections_R.intersection.apply(lambda x: right_loops[chr].id.iloc[x])

            # id_of_element_in_left -- id_of_intersecting_element_in_right_loops

            intersections = intersections_L.merge(intersections_R, on=["Loop_id","ids_column"], how="inner")
            idxs2 = intersections["ids_column"].values
            global_idxs = np.flatnonzero(idxs)[idxs2]
            result.iloc[global_idxs,0] = 1
        return result