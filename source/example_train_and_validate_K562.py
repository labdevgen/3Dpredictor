import logging
from Predictor import Predictor
from Weight_funcs_modul import *
import numpy as np
from functools import partial
from shared import decorate_oe2obs, oe2obs, decorate_return_coordinates_after_deletion, return_coordinates_after_deletion
import pandas as pd
from shared import Interval

logging.basicConfig(format='%(asctime)s %(message)s', datefmt='%I:%M:%S', level=logging.DEBUG)

cell_type = 'K562'
out_dir = "output/models/"
# If you use OE values for training and validating:
expected_folder = "input/expected/K562/5KB/"
deletion = 'no deletion'

# set all training files
training_files = [
    "output/K562/chr1_chr3_chr5_chr7_chr9training.RandOncontacts.gz.False.11.1500000.50001.1.1.cont_with_CTCF666406.25000",
        ]

for training_file in training_files:
    # set all validation files
    validation_files = [
        "output/K562/Interval_chr18_0_78000000validatingOrient.contacts.gz.False.11.1500000.50001.505429.all_cont.25000.txt",
    ]
    #Some examples of predictors filtering:
    #predictor.filter_predictors(".*CTCF.*|.*RNA.*", keep=False) #-- discard CTCF AND RNA
    #predictor.filter_predictors(".*CTCF.*", keep=False) #-- discard CTCF
    #predictor.filter_predictors(".*contact_dist.*|.*CTCF_W.*", keep=True) #-- keep only distance and CTCF in window

    #set contact_type. If "oe" - apply_log = False, elif "contacts" - apply_log=True
    for contact_type,apply_log in zip(["contacts"],[True]):
        #you can choose predictors for training and validating here:
        for (filter,keep),shortcut in zip(zip(["CTCF_|contact_dist|RNA"], [True]),
                                          ["CTCF_|contact_dist|RNA"]):
            # filter2, keep2 = zip(["OrientBlock|plus_orientation|minus_orientation"], [False])
            if contact_type == "oe":
                # You can choose weight function from the Weight_funcs_modul
                weightFuncs = [ones_like]
            else:
                weightFuncs = [ones_like]
            for h in range(2, 3): # h is a coefficient for matrix smoothing. Used for SCC estimation
                for weightFunc in weightFuncs:
                    predictor = Predictor()
                    predictor.read_data_predictors(training_file)
                    predictor.filter_predictors(filter, keep)
                    # for (filter2, keep2), shortcut in zip(
                    #         zip(["CTCF_ConvergentPair"], [False]),
                    #         ["CTCF|contact_dist|RNA"]):
                    #     predictor.filter_predictors(filter2, keep2)
                    trained_predictor = predictor.train(shortcut=shortcut, apply_log = apply_log,
                                                        weightsFunc = weightFunc, show_plot=False)
                    trained_predictor.out_dir = out_dir
                    trained_predictor.draw_Feature_importances(show_plot=False)
                    for validation_file in validation_files:
                        if apply_log:
                            trained_predictor.validate(validation_file, show_plot=False, cell_type=cell_type,
                                                       #                            transformation=
                                                       # [decorate_return_coordinates_after_deletion(return_coordinates_after_deletion, interval=deletion)],
                                                       validators=[trained_predictor.plot_juicebox,
                                                                   trained_predictor.decorate_scc(
                                                                       trained_predictor.scc, h=h, scc_file="scc.r",
                                                                       cell_type=cell_type),
                                                                   ], out_dir=out_dir)

                        # validators=[trained_predictor.r2score, trained_predictor.scc,
                        #             trained_predictor.plot_matrix, trained_predictor.plot_juicebox])
                        else:
                            trained_predictor.validate(validation_file, show_plot=False, cell_type=cell_type,
                                                       transformation=[decorate_oe2obs(oe2obs,
                                                                                       expected_folder=expected_folder,
                                                                                       cell_type=cell_type,
                                                                                       coeff_fname="coefficient." + cell_type + ".5KB.txt")],
                                                       validators=[trained_predictor.plot_juicebox,
                                                                   trained_predictor.decorate_scc(
                                                                       trained_predictor.scc, h=h, scc_file="scc.r",
                                                                       cell_type=cell_type,
                                                                       out_dir=out_dir + "scc/" + cell_type + "/"
                                                                   ),
                                                                   ])




