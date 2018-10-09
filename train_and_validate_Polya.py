import logging
from Predictor import Predictor
from Weight_funcs_modul import *
import numpy as np
from functools import partial

logging.basicConfig(format='%(asctime)s %(message)s', datefmt='%I:%M:%S', level=logging.DEBUG)


contact_type = ["oe"]#,"contacts"]
suffix = ".gz.3.1500000.50001.25000.25000.txt"
training_file = "out/GM12878/2018-10-02-training.RandOnchr2"
validation_files = [
     #"out/GM12878/Interval_chr10_10000000_60000000validatingOrient."
    #"out/Interval_chr1_100000000_110000000validatingOrient.",
    #"out/GM12878/Interval_chr10_65000000_70000000validatingOrient.",
    "out/GM12878/Interval_chr20_37000000_40000000validatingOrient.",
    #"out/Interval_chr2_47900000_53900000validatingOrient."
           ]

#training_file = "out/2018-09-23-trainingOrient.RandOnChr1."
#validation_files = [
#    "out/Interval_chr2_36000000_41000000validatingOrient.",
#    "out/Interval_chr2_47900000_53900000validatingOrient.",
#    "out/Interval_chr2_85000000_92500000validatingOrient."
           #]



#Some examples of predictors filtering:
#predictor.filter_predictors(".*CTCF.*|.*RNA.*", keep=False) #-- discard CTCF AND RNA
#predictor.filter_predictors(".*CTCF.*", keep=False) #-- discard CTCF
#predictor.filter_predictors(".*contact_dist.*|.*CTCF_W.*", keep=True) #-- keep only distance and CTCF in window

for contact_type,apply_log in zip(["oe"],[False]):
#for contact_type,apply_log in zip(["contacts"],[False]):
    for (filter,keep),shortcut in zip(zip([".*", "CTCF"]#"Loop","Loop|E1"]#,"Loop|E1|contact_dist","Loop|E1|contact_dist","Loop|E1|contact_dist|CTCF_L|CTCF_W|CTCF_R"] \
            ,[True, True]),#,False,True,True]),
                                           ["all", "CTCF_only"]):#, "no loop,no E1,no dist", \
                                            #"loop,E1,dist", "loop,E1,dist,LRW_CTCF,Nbl"]):
        if contact_type == "oe":
            weightFuncs = [#decorate_overweight_loops(overweight_loops,10), decorate_overweight_loops(overweight_loops,1000), \
            #                decorateContactWeither(contactWeitherFunction, coeff=5),decorateContactWeither(contactWeitherFunction, power=3), \
                            decorateContactWeither(contactWeitherFunction, power=3, abs=True)]#, abs_log]# decorateContactWeither(contactWeitherFunction, coeff=5, piecing=True), \
                            # decorateContactWeither(contactWeitherFunction, threshold=1.2, coeff=5, piecing=True, asymmetric=1), \
                            # decorateContactWeither(contactWeitherFunction, power=3,asymmetric=1)]
                            # [ones_like, array, abs_log, decorate_mult_abs_log(mult_abs_log,100), decorate_overweight_loops(overweight_loops,100)] \
                            # [decorateContactWeither(contactWeitherFunction, coeff=5),decorateContactWeither(contactWeitherFunction, power=3), \
                            # decorateContactWeither(contactWeitherFunction, power=3, abs=True), decorateContactWeither(contactWeitherFunction, coeff=5, piecing=True), \
                            # decorateContactWeither(contactWeitherFunction, threshold=1.2, coeff=5, piecing=True, asymmetric=1), \
                            # decorateContactWeither(contactWeitherFunction, power=3,asymmetric=1)]
                            # #decorate_mult_abs_log(mult_abs_log,10000)]
        else:
            weightFuncs = [ones_like]
        for weightFunc in weightFuncs:
            predictor = Predictor()
            predictor.read_data_predictors(training_file + contact_type + suffix)
            predictor.filter_predictors(filter, keep)
            trained_predictor = predictor.train(shortcut=shortcut, apply_log = apply_log,
                                                weightsFunc = weightFunc, show_plot=False)
            trained_predictor.out_dir = "out/models/"
            trained_predictor.draw_Feature_importances(show_plot=False)
            for validation_file in validation_files:
                trained_predictor.validate(validation_file + contact_type + suffix,
                                           show_plot=False)
                # trained_predictor.validate(validation_file + contact_type + suffix,
                #                            show_plot = False, validators=[trained_predictor.r2score, trained_predictor.plot_juicebox])
                #my_plot_matrix = partial(trained_predictor.plot_matrix,juicebox=True)


