import logging
from Predictor import Predictor
from Weight_funcs_modul import *
import numpy as np
from functools import partial
from shared import decorate_oe2obs, oe2obs

logging.basicConfig(format='%(asctime)s %(message)s', datefmt='%I:%M:%S', level=logging.DEBUG)
import pandas as pd

d = pd.read_csv("input/Loops/GM12878/GM12878.25000.loops")
expected_folder = "input/expected/GM12878/comb/"
not_cor_predictors_file =  "out/GM12878/validating_chrms_2/not_cor_predictors0.8"
loops_folder = "input/loops/GM12878/"
cell_type = 'GM12878'
suffix_val = ".gz.12.1500000.50001.1370946.25000.txt"
suffix = ".gz.12.1500000.50001.250000.25000.txt"#".gz.12.1500000.50001.250000.25000.txt"
training_files = [#"out/GM12878/validating_chrms_2/train23.txt",
                  # "out/GM12878/validating_chrms_2/train25000.txt",
                  "out/GM12878/validating_chrms_2/train250000.txt",
                  #  "out/GM12878/validating_chrms_2/train500000.txt",
                  #  "out/GM12878/validating_chrms_2/Interval_chr2_0_243175000validatingOrient.contacts.gz.8.1500000.50001.1850471.25000.txt"
                  ]
# coeff_fname = ""
validation_files = [
    "out/GM12878/validating_chrms_2/Interval_chr15_20000000_102500000validatingOrient.contacts.gz.8.1500000.50001.621051.25000.txt"
              ]

#Some examples of predictors filtering:
#predictor.filter_predictors(".*CTCF.*|.*RNA.*", keep=False) #-- discard CTCF AND RNA
#predictor.filter_predictors(".*CTCF.*", keep=False) #-- discard CTCF
#predictor.filter_predictors(".*contact_dist.*|.*CTCF_W.*", keep=True) #-- keep only distance and CTCF in window

#choose predictors from file
# predictors_data=pd.read_csv(not_cor_predictors_file,sep="\t")
# predictors=str("|".join(list(predictors_data["predictors"])))
# validation_data = pd.read_csv(validation_files[0] + "contacts" +suffix_val)
# print(validation_data.keys())
for training_file in training_files:
    for contact_type,apply_log in zip(["contacts"],[True]):
    #for contact_type,apply_log in zip(["contacts"],[False]):
        for (filter,keep),shortcut in zip(zip([".*"]#".*"]#, "CTCF|contact_dist", "CTCF|contact_dist|RNA|CAGE"]#".*", "CTCF"]#"Loop","Loop|E1"]#,"Loop|E1|contact_dist","Loop|E1|contact_dist","Loop|E1|contact_dist|CTCF_L|CTCF_W|CTCF_R"] \
                ,[True]),#,True,True]),#,False,True,True]),
                                               [
                                                "all"]):#, "CTCF, dist, RNA, cage"]):#, "no loop,no E1,no dist", \
                                                #"loop,E1,dist", "loop,E1,dist,LRW_CTCF,Nbl"]):
            if contact_type == "oe":
                weightFuncs = [ones_like]#, decorate_overweight_loops(overweight_loops,1000, loop_file="input/Loops/GM12878/GM12878.25000.loops") ]#, array, abs_log,
                               # decorate_overweight_loops(overweight_loops,10, loop_file="input/Loops/GM12878/GM12878.25000.loops"),
                               # decorate_overweight_loops(overweight_loops,1000, loop_file="input/Loops/GM12878/GM12878.25000.loops"), \
                               # decorateContactWeither(contactWeitherFunction, threshold=1.2, coeff=5, piecing=True, asymmetric=1),\
                               #                        decorateContactWeither(contactWeitherFunction, power=3, abs=True)]
                                #decorateContactWeither(contactWeitherFunction, power=3, abs=True)]#, abs_log]# decorateContactWeither(contactWeitherFunction, coeff=5, piecing=True), \
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
                # if train_on_few_chrms:
                #     predictor1 = Predictor()
                #     predictor1.read_data_predictors(training_file + contact_type + suffix)
                #     predictor1.filter_predictors(filter, keep)
                #     predictor2 = Predictor()
                #     predictor2.read_data_predictors(training_file + contact_type + suffix)
                #     predictor2.filter_predictors(filter, keep)
                #     predictor =
                # else:

                predictor = Predictor()
                predictor.read_data_predictors(training_file)
                predictor.filter_predictors(filter, keep)
                trained_predictor = predictor.train(shortcut=shortcut, apply_log = apply_log,
                                                    weightsFunc = weightFunc, show_plot=False)
                trained_predictor.out_dir = "out/models/"
                trained_predictor.draw_Feature_importances(show_plot=False)
                for validation_file in validation_files:
                    # trained_predictor.validate(validation_file + contact_type + suffix, show_plot=False,
                    #                            transformation=decorate_oe2obs(oe2obs, expected_folder=expected_folder, cell_type=cell_type))
                                               #transformation=xxxx)
                    for h in range(2, 3):
                        if apply_log:
                            trained_predictor.validate(validation_file , show_plot=False,
                                                       validators=[trained_predictor.decorate_scc(trained_predictor.scc, h=h, loop_file="input/Loops/GM12878/GM12878.25000.loops"),
                                                                   trained_predictor.plot_juicebox])
                        else:
                            trained_predictor.validate(validation_file, show_plot = False,
                                                   transformation=decorate_oe2obs(oe2obs, expected_folder=expected_folder, cell_type=cell_type, coeff_fname=""),
                                                   validators=[trained_predictor.decorate_scc(trained_predictor.scc, h=h, loop_file="input/Loops/GM12878/GM12878.25000.loops"),
                                                               trained_predictor.plot_juicebox])
                                                   # trained_predictor.r2score,trained_predictor.plot_matrix, trained_predictor.scc,
                                                   #         trained_predictor.plot_juicebox])
                                               #transformation=decorate_oe2obs(oe2obs, expected_folder=expected_folder, cell_type=cell_type) )#trained_predictor.plot_juicebox,
                    #my_plot_matrix = partial(trained_predictor.plot_matrix,juicebox=True)


