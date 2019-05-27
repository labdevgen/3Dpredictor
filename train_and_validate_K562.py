import logging
from Predictor import Predictor
from Weight_funcs_modul import *
import numpy as np
from functools import partial
from shared import decorate_oe2obs, oe2obs, decorate_return_coordinates_after_deletion, return_coordinates_after_deletion
import pandas as pd
from shared import Interval

logging.basicConfig(format='%(asctime)s %(message)s', datefmt='%I:%M:%S', level=logging.DEBUG)

expected_folder = "/mnt/scratch/ws/psbelokopytova/201905031108polinaB/3DPredictor/input/expected/K562/5KB/"
cell_type = 'K562'
# contact_type = ["oe"]
# not_cor_predictors_file = "/mnt/scratch/ws/psbelokopytova/201905031108polinaB/3DPredictor/out/Hepat/validating_chrms_1/not_cor_predictors0.8"
# suffix = ".gz.8.1500000.50001.25000.25000.txt"

# suffix_val = ".gz.1000000.50001.500000.25000.txt"#".gz.12.1000000.50001.1117748.25000.txt"
# deletion =Interval("chr10",60755595, 60761099)
deletion = "no_deletion"
training_files = [
# "/mnt/scratch/ws/psbelokopytova/201905031108polinaB/3DPredictor/out/K562/5KB/training.RandOnchr10contacts.gz.11.1500000.50001.100000.5000.txt",
#     "/mnt/scratch/ws/psbelokopytova/201905031108polinaB/3DPredictor/out/K562/5KB/training.RandOnchr5contacts.gz.11.1500000.50001.1000000.5000.txt",
    "/mnt/scratch/ws/psbelokopytova/201905031108polinaB/3DPredictor/out/K562/5KB/training.RandOnchr10oe.gz.11.1500000.50001.100000.5000.txt",
#     "/mnt/scratch/ws/psbelokopytova/201905031108polinaB/3DPredictor/out/K562/5KB/training.RandOnchr10contacts.gz.True.11.1500000.10001.100000.5000.txt",
#         "/mnt/scratch/ws/psbelokopytova/201905031108polinaB/3DPredictor/out/K562/5KB/training.RandOnchr10contacts.gz.False.11.1500000.10001.100000.5000.txt",
#     "/mnt/scratch/ws/psbelokopytova/201905031108polinaB/3DPredictor/out/K562/5KB/1.1withCTCF_wotraining.RandOnchr10contacts.gz.False.11.1500000.10001.100000.5000.txt",
#     "/mnt/scratch/ws/psbelokopytova/201905031108polinaB/3DPredictor/out/K562/5KB/1.1_1bins_withCTCF_wotraining.RandOnchr10contacts.gz.False.11.1500000.10001.100000.5000.txt",
# "/mnt/scratch/ws/psbelokopytova/201905031108polinaB/3DPredictor/out/K562/5KB/1.2_withCTCF_wotraining.RandOnchr10contacts.gz.False.11.1500000.10001.100000.5000.txt"
        ]

# predictors_data=pd.read_csv(not_cor_predictors_file,sep="\t")
# predictors=str("|".join(list(predictors_data["predictors"])))
# print(predictors)
# i=1
# assert i>1
for training_file in training_files:
    validation_files = [
# "/mnt/scratch/ws/psbelokopytova/201905031108polinaB/3DPredictor/out/K562/5KB/Interval_chr2_118000000_129000000validatingOrient.contacts.gz.11.1500000.50001.1000000.5000.txt",
# "/mnt/scratch/ws/psbelokopytova/201905031108polinaB/3DPredictor/out/K562/5KB/Interval_chr4_86000000_114000000validatingOrient.contacts.gz.11.1500000.50001.1000000.5000.txt",
#         "/mnt/scratch/ws/psbelokopytova/201905031108polinaB/3DPredictor/out/K562/5KB/Interval_chr2_118000000_129000000validatingOrient.contacts.gz.True.11.1500000.10001.14448381.5000.all_contacts.txt",
#         "/mnt/scratch/ws/psbelokopytova/201905031108polinaB/3DPredictor/out/K562/5KB/Interval_chr2_118000000_129000000validatingOrient.contacts.gz.False.11.1500000.10001.24241897.5000.txt",
        "/mnt/scratch/ws/psbelokopytova/201905031108polinaB/3DPredictor/out/K562/5KB/Interval_chr2_118000000_129000000validatingOrient.oe.gz.False.11.1500000.10001.24241897.5000.txt"
    ]
    #Some examples of predictors filtering:
    #predictor.filter_predictors(".*CTCF.*|.*RNA.*", keep=False) #-- discard CTCF AND RNA
    #predictor.filter_predictors(".*CTCF.*", keep=False) #-- discard CTCF
    #predictor.filter_predictors(".*contact_dist.*|.*CTCF_W.*", keep=True) #-- keep only distance and CTCF in window

    for contact_type,apply_log in zip(["oe"],[False]):
    #for contact_type,apply_log in zip(["contacts"],[False]):
        for (filter,keep),shortcut in zip(zip(["CTCF_|contact_dist|RNA"], [True]),
                                          ["CTCF_|contact_dist|RNA"]):#, "CTCF", "CTCF|contact_dist"]#"Loop","Loop|E1"]#,"Loop|E1|contact_dist","Loop|E1|contact_dist","Loop|E1|contact_dist|CTCF_L|CTCF_W|CTCF_R"] \
                # [True],#, True, True]),#,False,True,True]),
                #                                ["all"])):#,"CTCF_only", "CTCF and contact_dist" ]):#, "no loop,no E1,no dist", \
                #                                 #"loop,E1,dist", "loop,E1,dist,LRW_CTCF,Nbl"]):
            # filter2, keep2 = zip(["OrientBlock|plus_orientation|minus_orientation"], [False])
            if contact_type == "oe":
                weightFuncs = [abs_log]
                               #abs_log, decorate_mult_abs_log(mult_abs_log,100),decorate_overweight_loops(overweight_loops,10),
                    # decorate_overweight_loops(overweight_loops,1000), \
                                #decorateContactWeither(contactWeitherFunction, coeff=5),decorateContactWeither(contactWeitherFunction, power=3),
                                # decorateContactWeither(contactWeitherFunction, power=3, abs=True), abs_log, decorateContactWeither(contactWeitherFunction, coeff=5, piecing=True), \
                                #  decorateContactWeither(contactWeitherFunction, threshold=1.2, coeff=5, piecing=True, asymmetric=1)] \
                                # decorateContactWeither(contactWeitherFunction, power=3,asymmetric=1)]
                                # [ones_like, array, abs_log, decorate_mult_abs_log(mult_abs_log,100), decorate_overweight_loops(overweight_loops,100)] \
                                # [decorateContactWeither(contactWeitherFunction, coeff=5),decorateContactWeither(contactWeitherFunction, power=3), \
                                # decorateContactWeither(contactWeitherFunction, power=3, abs=True), decorateContactWeither(contactWeitherFunction, coeff=5, piecing=True), \
                                # decorateContactWeither(contactWeitherFunction, threshold=1.2, coeff=5, piecing=True, asymmetric=1), \
                                # decorateContactWeither(contactWeitherFunction, power=3,asymmetric=1)]
                                # #decorate_mult_abs_log(mult_abs_log,10000)]
            else:
                weightFuncs = [ones_like]
            for h in range(2, 3):
                for weightFunc in weightFuncs:
                    predictor = Predictor()
                    predictor.read_data_predictors(training_file)# + contact_type + suffix)
                    predictor.filter_predictors(filter, keep)
                    # for (filter2, keep2), shortcut in zip(
                    #         zip(["CTCF_ConvergentPair"], [False]),
                    #         ["CTCF|contact_dist|RNA"]):
                    #     predictor.filter_predictors(filter2, keep2)
                    trained_predictor = predictor.train(shortcut=shortcut, apply_log = apply_log,
                                                        weightsFunc = weightFunc, show_plot=False)
                    trained_predictor.out_dir = "/mnt/scratch/ws/psbelokopytova/201905031108polinaB/3DPredictor/out/models/"
                    trained_predictor.draw_Feature_importances(show_plot=False)
                    for validation_file in validation_files:
                        if apply_log:
                            trained_predictor.validate(validation_file, show_plot = False,cell_type=cell_type,
                            #                            transformation=
                            # [decorate_return_coordinates_after_deletion(return_coordinates_after_deletion, interval=deletion)],
                                                       validators=[trained_predictor.plot_juicebox,
                                                           trained_predictor.decorate_scc(
                                                           trained_predictor.scc, h=h, scc_file="scc.r",cell_type=cell_type,
                                                          ),
                                                           ])

                                                   # validators=[trained_predictor.r2score, trained_predictor.scc,
                                                   #             trained_predictor.plot_matrix, trained_predictor.plot_juicebox])
                        else:
                            trained_predictor.validate(validation_file, show_plot=False,cell_type=cell_type,
                                                       transformation=[decorate_oe2obs(oe2obs,
                                                                                      expected_folder=expected_folder,
                                                                                      cell_type=cell_type,
                                                                                      coeff_fname="coefficient." + cell_type + ".5KB.txt")],
                                                       validators=[trained_predictor.plot_juicebox,
                                                                   # trained_predictor.decorate_scc(
                                                           # trained_predictor.scc, h=h, scc_file="scc.r",cell_type=cell_type,
                                                           #         ),
                                                           ])
                                                       # validators=[trained_predictor.r2score, trained_predictor.scc,
                                                       #             trained_predictor.plot_matrix, trained_predictor.plot_juicebox],
                                                       # )
                        #my_plot_matrix = partial(trained_predictor.plot_matrix,juicebox=True)


