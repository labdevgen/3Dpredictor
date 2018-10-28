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
loops_folder = "input/loops/GM12878/"
cell_type = 'GM12878'
contact_type = ["oe", "contacts"]
suffix = ".gz.12.1500000.50001.25000.25000.txt"
training_file = "out/GM12878/2018-10-11-training.RandOnchr1"
validation_files = [
    # "out/GM12878/Interval_chr1_0_249225000validatingOrient.oe.gz.12.1500000.50001.1728185.25000.txt",
    # "out/GM12878/Interval_chr2_0_243175000validatingOrient.oe.gz.12.1500000.50001.1850471.25000.txt",
    # "out/GM12878/Interval_chr3_50000_197925000validatingOrient.oe.gz.12.1500000.50001.1518948.25000.txt",
    # "out/GM12878/Interval_chr4_0_191025000validatingOrient.oe.gz.12.1500000.50001.1458075.25000.txt",
    #"out/GM12878/Interval_chr5_0_180900000validatingOrient.",
    # "out/GM12878/Interval_chr6_75000_171050000validatingOrient.oe.gz.12.1500000.50001.1298441.25000.txt",
    "out/GM12878/Interval_chr7_0_159125000validatingOrient.",
    # "out/GM12878/Interval_chr8_25000_146300000validatingOrient.oe.gz.12.1500000.50001.1097984.25000.txt",
    # "out/GM12878/Interval_chr9_0_141100000validatingOrient.oe.gz.12.1500000.50001.863649.25000.txt",
    # "out/GM12878/Interval_chr10_50000_135500000validatingOrient.oe.gz.12.1500000.50001.991450.25000",
    # "out/GM12878/Interval_chr11_100000_134925000validatingOrient.oe.gz.12.1500000.50001.1006983.25000.txt",
    # "out/GM12878/Interval_chr12_50000_133825000validatingOrient.oe.gz.12.1500000.50001.1003008.25000.txt",
    # "out/GM12878/Interval_chr13_19000000_115100000validatingOrient.oe.gz.12.1500000.50001.743030.25000.txt",
    # "out/GM12878/Interval_chr14_19000000_107275000validatingOrient.oe.gz.12.1500000.50001.683266.25000.txt",
    # "out/GM12878/Interval_chr15_20000000_102500000validatingOrient.oe.gz.12.1500000.50001.621051.25000.txt",
    # "out/GM12878/Interval_chr16_50000_90275000validatingOrient.oe.gz.12.1500000.50001.578746.25000.txt",
    # "out/GM12878/Interval_chr17_0_81175000validatingOrient.oe.gz.12.1500000.50001.577715.25000.txt",
    # "out/GM12878/Interval_chr18_0_78000000validatingOrient.oe.gz.12.1500000.50001.559471.25000.txt",
    # "out/GM12878/Interval_chr19_75000_59100000validatingOrient.oe.gz.12.1500000.50001.406785.25000.txt",
    # "out/GM12878/Interval_chr20_50000_62950000validatingOrient.oe.gz.12.1500000.50001.438352.25000.txt",
    # "out/GM12878/Interval_chr21_9400000_48100000validatingOrient.oe.gz.12.1500000.50001.251796.25000.txt",
    # "out/GM12878/Interval_chrX_50000_155250000validatingOrient.oe.gz.12.1500000.50001.1157673.25000.txt"
              ]

#Some examples of predictors filtering:
#predictor.filter_predictors(".*CTCF.*|.*RNA.*", keep=False) #-- discard CTCF AND RNA
#predictor.filter_predictors(".*CTCF.*", keep=False) #-- discard CTCF
#predictor.filter_predictors(".*contact_dist.*|.*CTCF_W.*", keep=True) #-- keep only distance and CTCF in window

for contact_type,apply_log in zip(["contacts"],[True]):
#for contact_type,apply_log in zip(["contacts"],[False]):
    for (filter,keep),shortcut in zip(zip(["CTCF|contact_dist|RNA|CAGE"]#".*"]#, "CTCF|contact_dist", "CTCF|contact_dist|RNA|CAGE"]#".*", "CTCF"]#"Loop","Loop|E1"]#,"Loop|E1|contact_dist","Loop|E1|contact_dist","Loop|E1|contact_dist|CTCF_L|CTCF_W|CTCF_R"] \
            ,[True]),#,True,True]),#,False,True,True]),
                                           ["CTCF,contact_dist,RNA,CAGE"]):#, "CTCF, dist, RNA, cage"]):#, "no loop,no E1,no dist", \
                                            #"loop,E1,dist", "loop,E1,dist,LRW_CTCF,Nbl"]):
        if contact_type == "oe":
            weightFuncs = [ones_like]#, array, abs_log,
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
            predictor = Predictor()
            predictor.read_data_predictors(training_file + contact_type + suffix)
            predictor.filter_predictors(filter, keep)
            trained_predictor = predictor.train(shortcut=shortcut, apply_log = apply_log,
                                                weightsFunc = weightFunc, show_plot=False)
            trained_predictor.out_dir = "out/models/"
            trained_predictor.draw_Feature_importances(show_plot=False)
            for validation_file in validation_files:
                # trained_predictor.validate(validation_file + contact_type + suffix, show_plot=False,
                #                            transformation=decorate_oe2obs(oe2obs, expected_folder=expected_folder, cell_type=cell_type))
                                           #transformation=xxxx)
                for h in range(100, 101):
                    if apply_log:
                        trained_predictor.validate(validation_file + contact_type + ".gz.12.1500000.50001.1191309.25000.txt", show_plot=False,
                                                   validators=[trained_predictor.decorate_scc(trained_predictor.scc, h=h, loop_file="input/Loops/GM12878/GM12878.25000.loops"),
                                                               trained_predictor.plot_juicebox])
                    else:
                        trained_predictor.validate(validation_file + contact_type + ".gz.12.1500000.50001.1191309.25000.txt" , show_plot = False,
                                               transformation=decorate_oe2obs(oe2obs, expected_folder=expected_folder, cell_type=cell_type),
                                               validators=[trained_predictor.decorate_scc(trained_predictor.scc, h=h, loop_file="input/Loops/GM12878/GM12878.25000.loops"),
                                                           trained_predictor.plot_juicebox])
                                               # trained_predictor.r2score,trained_predictor.plot_matrix, trained_predictor.scc,
                                               #         trained_predictor.plot_juicebox])
                                           #transformation=decorate_oe2obs(oe2obs, expected_folder=expected_folder, cell_type=cell_type) )#trained_predictor.plot_juicebox,
                #my_plot_matrix = partial(trained_predictor.plot_matrix,juicebox=True)


