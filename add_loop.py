from LoopReader import LoopReader
from VectPredictorGenerators import loopsPredictorGenerator

def add_loop(validation_data,loop_file):
    window_size = 25000
    loopsReader = LoopReader(loop_file)
    loopsReader.read_loops()
    loopspg = loopsPredictorGenerator(loopsReader, window_size)
    contacts = validation_data[["chr", "contact_st", "contact_en", "contact_count"]].copy()
    isLoop_df = loopspg.get_predictors(contacts)
    validation_data["IsLoop"] = isLoop_df["IsLoop"].values
