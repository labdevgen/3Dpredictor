import logging
from Predictor import Predictor

logging.basicConfig(format='%(asctime)s %(message)s', datefmt='%I:%M:%S', level=logging.DEBUG)

contact_type = ["oe","contacts"]
suffix = ".gz.1000000.50001.500000.25000.txt"
training_file = "out/2018-09-23-trainingOrient.RandOnChr1."
validation_files = [
    "out/Interval_chr2_36000000_41000000validatingOrient.",
    "out/Interval_chr2_47900000_53900000validatingOrient.",
    "out/Interval_chr2_85000000_92500000validatingOrient."
            ]


#Some examples of predictors filtering:
#predictor.filter_predictors(".*CTCF.*|.*RNA.*", keep=False) #-- discard CTCF AND RNA
#predictor.filter_predictors(".*CTCF.*", keep=False) #-- discard CTCF
#predictor.filter_predictors(".*contact_dist.*|.*CTCF_W.*", keep=True) #-- keep only distance and CTCF in window

for contact_type,apply_log in zip(["contacts","oe"],[True,False]):
    for (filter,keep),shortcut in zip(zip([".*","E1"],[True,False]),
                                           ["all","no E1"]):
        predictor = Predictor()
        predictor.read_data_predictors(training_file + contact_type + suffix)
        predictor.filter_predictors(filter, keep)
        trained_predictor = predictor.train(shortcut=shortcut, apply_log = apply_log)
        for validation_file in validation_files:
            trained_predictor.validate(validation_file + contact_type + suffix)

