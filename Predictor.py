# Written by Minja, 2018-09
# Class Predictor, to train and validate models
# Basic functionality:
# 1. Loading files with predictors
# 2. Drop different sets of predictors
# 3. Train and validate models
# 4. Perform estimation and visualization of results (TODO)

import logging,os,re,pickle,sklearn
import numpy as np
import pandas as pd
from numpy import int64,float32
from sklearn import ensemble
from shared import str2hash,Interval,write_XML
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matrix_plotter import MatrixPlotter
from collections import OrderedDict
from matplot2hic import MatPlot2HiC
import subprocess
from functools import partial
from add_loop import add_loop
import numpy as np
import math
import xgboost
#from catboost import CatBoostRegressor
from shared import Interval, get_bin_size
import datetime

def ones_like(contacts,*args):
    return [1]*len(contacts)

def equal(contacts, data,**kwargs):
    if kwargs['data_type']=='predicted':
        return contacts
    elif kwargs['data_type'] == 'validation':
            return data

class Predictor(object):
    def __setattr__(self, key, value):
        if key == "predictors":
            self.trained_model = None
        self.__dict__[key] = value

    def toXMLDict(self):
        try:
            self.apply_log
        except:
            raise Exception("Please read validation data firts")
        result = OrderedDict()
        result["shortcut"] = self.shortcut
        result["input_file"] = self.input_file
        result["predictors"] = ".".join(self.predictors)
        result["algrorithm"] = str(self.alg.__class__.__name__)
        result["algorithm_params"] = str(self.alg_params)
        result["apply_log"] = str(self.apply_log)
        result["weights_func"] = str(self.weightsFunc.__name__)
        return result

    def __represent_validation__(self):
        try:
            self.validation_file
        except:
            raise Exception("Please read validation data firts")
        return "model"+str(self)+".validation."+"."+str(self.transformation_for_validation_data)+"."+self.cell_type+"."+\
                        str2hash(os.path.basename(self.validation_file))

    def __repr__(self): # Representation should reflect all paramteres
                        # so that 2 models having same repr are identical
        try:
            self.apply_log
        except:
            raise Exception("Please train model first")
        properties = self.toXMLDict()
        representation = ".".join(list(properties.values()))
        return str2hash(representation)

    def __init__(self):
        self.predictors = None
        self.contacts = None
        self.input_file = None
        self.apply_log = None
        self.shortcut = "model"
        self.alg = None
        self.trained_model = None
        self.constant_nonpredictors = ["chr", "contact_st", "contact_en", "contact_count"]

    def draw_Feature_importances(self, show_plot=True):
        if self.trained_model is None:
            logging.getLogger(__name__).error("Please train model first")
            return
        dump_path = os.path.join(self.out_dir,self.representation)
        try:
            self.trained_model.feature_importances_
        except:
            logging.getLogger(__name__).warning(
                "Features importance is not avaliable for the model " + str(self.trained_model.__class__.__name__))
            return 0
        # create file with feature importances
        importances = pd.Series(self.trained_model.feature_importances_, index = self.predictors ).sort_values(ascending=False)
        #print(importances)
        importances.to_csv(dump_path + ".featureImportance.txt", sep='\t')
        plt.plot(range(len(self.trained_model.feature_importances_)),
                 self.trained_model.feature_importances_, marker="o")
        plt.xticks(range(len(self.trained_model.feature_importances_)),
                   self.predictors, rotation='vertical')
        plt.setp(plt.gca().get_xticklabels(), rotation='vertical', fontsize=7)
        plt.grid(which='both',axis="x",ls=":")
        plt.xlabel("Feature")
        plt.ylabel("Importance")

        # these are matplotlib.patch.Patch properties
        props = dict(boxstyle='round', facecolor='wheat', alpha=0.5)
        # place a text box in upper left in axes coords
        xml = self.toXMLDict()
        textstr = "\n".join(key + " " +val for key,val in xml.items() if key != "predictors" )
        plt.gca().text(0.05, 0.95, textstr , transform=plt.gca().transAxes, fontsize=8,
                verticalalignment='top', bbox=props)

        plt.title("Features importance for model " + self.representation)
        plt.savefig(dump_path + ".FeatureImportance.png", dpi=300)
        if show_plot:
            plt.show()
        plt.clf()

    # Train model
    # if dump = True, save it to file dump_file
    # show_plot = True/False show features importance plot
    # returns class instance with trained_model object
    def train(self,
              alg = xgboost.XGBRegressor(n_estimators=100,max_depth=9,subsample=0.7
              # alg = sklearn.ensemble.GradientBoostingRegressor(n_estimators=100, max_depth=9, subsample=0.7
              # alg=CatBoostRegressor(n_estimators=100, max_depth=9, subsample=0.7
              # alg=sklearn.ensemble.RandomForestRegressor(n_estimators=100, max_depth=9
                                   ),
              shortcut = "model", apply_log = True,
              dump = True, out_dir = "output/models/",
              weightsFunc = ones_like,
              show_plot = True,
              *args, **kwargs):

        # Check that we have got data file
        try:
            self.input_file
            self.predictors
        except:
            raise Exception("Please read the data and set predictors first")

        # Save paramters to be able to hash model name
        self.predictors = sorted(self.predictors)
        self.alg = alg
        self.alg_params = ".".join([str(k)+"_"+str(v) for k,v in sorted(alg.get_params().items())])
        self.shortcut = shortcut
        self.apply_log = apply_log
        self.weightsFunc = weightsFunc
        self.out_dit = out_dir


        # remove validation data since we are going to dump instance and
        # do not want file to be large
        try:
            del self.validation_file
            del self.validation_data
        except:
            pass

        # First try to load model dump
        self.representation = str(self)
        dump_path = os.path.join(out_dir,self.representation)
        if os.path.exists(dump_path):
            logging.info("Found dump for model " + dump_path)
            return pickle.load(open(dump_path,"rb"))
        else:
            # read data
            self.input_data = self.read_file(self.input_file)
            self.train_chrms = set(self.input_data["chr"].values)
            print("!!!!!!!!1", self.train_chrms)
            self.input_data.fillna(value=0, inplace=True)
            self.contacts = np.array(self.input_data["contact_count"].values)
            print("train contacts", self.contacts)

            # fit new model
            if apply_log:
                self.contacts = np.log(self.contacts)
            logging.getLogger(__name__).info("Fitting model")
            alg.fit(self.input_data[self.predictors],self.contacts,
                    sample_weight=self.weightsFunc(self.contacts,self.input_data))
            self.trained_model = alg
            if dump:
                logging.getLogger(__name__).info("Saving to file "+dump_path)
                #Remove large variables before dump
                del self.contacts
                del self.input_data
                pickle.dump(self,open(dump_path,"wb"))
                write_XML(self.toXMLDict(),str(self),dump_path+".xml")

        logging.getLogger(__name__).info("Done")
        return self


    def r2score(self,validation_data,predicted,out_dir,**kwargs):
        real = validation_data["contact_count"].values
        r2 = sklearn.metrics.r2_score(predicted, real)

        # Plot r2
        subset = max(len(predicted) // 5000, 1)
        if "contact_dist" in validation_data.columns.get_values():
            plt.scatter(predicted[::subset], real[::subset],
                        c=validation_data["contact_dist"].values[::subset])
        else:
            plt.scatter(predicted[::subset], real[::subset])
        plt.title("Test: Predicted vs real, r^2 score = " + str(r2) + "\n" + self.__represent_validation__())
        plt.xlabel("Predicted Contact")
        plt.ylabel("Real Contact")
        print("Saving file " + os.path.join(out_dir,self.__represent_validation__()) + ".r2scatter.png")
        plt.savefig(os.path.join(out_dir,self.__represent_validation__()) + ".r2scatter.png", dpi=300)
        if not ("show_plot" in kwargs) or kwargs["show_plot"]:
            plt.show()
        plt.clf()

    def plot_matrix(self,validation_data,predicted,out_dir,**kwargs):
        predicted_data = validation_data.copy(deep=True)
        predicted_data["contact_count"] = predicted
        mp = MatrixPlotter()
        mp.set_data(validation_data)
        mp.set_control(predicted_data)
        matrix = mp.getMatrix4plot(Interval(validation_data["chr"].iloc[0],
                                            min(validation_data["contact_st"].values),
                                            max(validation_data["contact_en"].values)))
        #if not self.apply_log:
        matrix = np.log(matrix)

        tick_pos, tick_labels = mp.get_bins_strart_labels(maxTicksNumber=15)
        plt.xticks(tick_pos, tick_labels, rotation=45)
        plt.imshow(matrix, cmap="OrRd")
        plt.title(self.__represent_validation__())
        # these are matplotlib.patch.Patch properties
        props = dict(boxstyle='round', facecolor='wheat', alpha=0.5)
        # place a text box in upper left in axes coords
        xml = self.toXMLDict()
        textstr = "\n".join(key + " " + val for key, val in xml.items() if key != "predictors")
        plt.gca().text(0.05, 0.95, textstr, transform=plt.gca().transAxes, fontsize=6,
                       verticalalignment='top', bbox=props)

        plt.imsave(os.path.join(out_dir,self.__represent_validation__()) + ".matrix.png", matrix,
                    cmap="OrRd", dpi=600)
        if not ("show_plot" in kwargs) or kwargs["show_plot"]:
            plt.show()
        plt.clf()

    def scc(self,validation_data,predicted,out_dir,**kwargs):
        # if self.apply_log:
        #     print(validation_data["contact_count"])
        #     print(predicted)
        #     validation_data["contact_count"] = validation_data["contact_count"].apply(lambda x: math.exp(x))
        #     predicted = np.exp(np.array(predicted))
        #     print(validation_data["contact_count"])
        #     print(predicted)
        # print(validation_data["chr"])
        chromosome = str(validation_data["chr"][1])
        binsize = str(get_bin_size(validation_data))
        # print("chromosome", chromosome)
        if "h" not in kwargs:
            kwargs["h"] = 2
        else:
            logging.info("for scc using h = " + str(kwargs["h"]))
        if "loop_file" not in kwargs:
            d = pd.concat([validation_data["contact_st"],validation_data["contact_en"],validation_data["contact_count"],pd.DataFrame(predicted)], axis=1)
        else:
            add_loop(validation_data, kwargs["loop_file"])
            d = pd.concat([validation_data["contact_st"],validation_data["contact_en"],validation_data["contact_count"],pd.DataFrame(predicted),validation_data["IsLoop"]], axis=1)
        out_fname = os.path.join(out_dir+"scc/",chromosome + "." + binsize+"."+ self.__represent_validation__()) + ".scc"
        pd.DataFrame.to_csv(d, out_fname, sep=" ", index=False)
        logging.info(datetime.datetime.now())
        if "p_file" not in kwargs or "e_file" not in kwargs:
            if "interact_pr_en"  not in kwargs:
                out = subprocess.check_output(["Rscript", kwargs["scc_file"], out_fname, str(kwargs["h"]), chromosome])
            else:
                out = subprocess.check_output(["Rscript", kwargs["scc_file"], out_fname,str(kwargs["h"]), chromosome, kwargs["interact_pr_en"]])
        else:
            out = subprocess.check_output(["Rscript", kwargs["scc_file"], out_fname, str(kwargs["h"]), chromosome, kwargs["p_file"], kwargs["e_file"]])
        print(str(out))


    def decorate_scc(self, func, h, scc_file,cell_type,**kwargs):
        if "loop_file" in kwargs:
            if "p_file" in kwargs and "e_file" in kwargs:
                result = partial(func, h=h, scc_file=scc_file, loop_file=kwargs["loop_file"], p_file=kwargs["p_file"], e_file=["e_file"])
            else:
                result = partial(func, h=h, scc_file=scc_file, loop_file=kwargs["loop_file"])
        elif "loop_file" not in kwargs:
            if "p_file" in kwargs and "e_file" in kwargs:
                # print(kwargs["p_file"])
                result = partial(func, h=h, scc_file=scc_file, p_file=kwargs["p_file"], e_file=kwargs["e_file"])
            elif "interact_pr_en" in kwargs:
                result = partial(func, h=h, scc_file=scc_file, interact_pr_en=kwargs["interact_pr_en"])
            else:
                result = partial(func, h=h, scc_file=scc_file)
        self.h_for_scc = "h="+str(h)
        self.cell_type = cell_type
        result.__name__ = str(h) + cell_type+ func.__name__
        return result

    def plot_juicebox(self,validation_data,predicted,out_dir,**kwargs):
        print("predicted", predicted)
        print("validation cont_count", validation_data["contact_count"])
        predicted_data = validation_data.copy(deep=True)
        predicted_data["contact_count"] = predicted
        # print(predicted_data.query("contact_count >10000")['contact_count'])
        out_dir = "output/hic_files"
        mp = MatrixPlotter()
        mp.set_data(predicted_data)
        mp.set_control(validation_data)
        chromosome = str(validation_data["chr"][1])
        #mp.set_apply_log(self.apply_log)
        MatPlot2HiC(mp, chromosome + "." + self.__represent_validation__(), out_dir)

    def return_predicted(self, validation_data,predicted,out_dir,**kwargs):
        return [validation_data, predicted]
    # Validate model
    def validate(self, validation_file,
                 out_dir = "output/",
                 validators = None,
                 transformation = [equal],
                 cell_type=None,
                 **kwargs):
        # validation_file - file with validation data
        # out_dir - directory to save output produced during validation
        # validators - list of functions used to validate each region, i.e. funcs to calc r2score. scc and etc.
        # transformation - function to apply to contact values before running validation funcs
        # i.e. if using o/e values it can transform it back to contacts based on expected values
        # kwargs will be passed to validation functions

        validators = validators if validators is not None else [self.r2score,self.plot_matrix,self.scc]
        self.cell_type = cell_type
        self.validation_file = validation_file
        self.validation_data = self.read_file(validation_file)
        self.validation_data.fillna(value=0, inplace=True)
        # check that train chrms not in validate
        validate_chrms = set(self.validation_data["chr"].values)
        assert [chr not in self.train_chrms for chr in validate_chrms]
        self.transformation_for_validation_data = ""
        self.predicted = self.trained_model.predict(self.validation_data[self.predictors])
        print("!!!!!!!!!!!predicted", self.predicted)
        for transformation_function in transformation:
            print(transformation_function.__name__)
            print(self.validation_data)
            self.transformation_for_validation_data+=transformation_function.__name__
            self.predicted = transformation_function(self.predicted,
                                        data=self.validation_data, data_type="predicted")
            self.validation_data = transformation_function(self.validation_data["contact_count"].values,
                                                               data=self.validation_data, data_type="validation")

        print(self.validation_data)
        #do this for validation with observed contacts
        if self.apply_log:
            self.predicted = np.exp(self.predicted)

        for validataion_function in validators:
            validataion_function(self.validation_data.copy(),self.predicted.copy(),
                                 out_dir = out_dir, **kwargs)


    # Read header of predictors file, get list of avaliable predictors
    # Returns list of predictors
    def get_avaliable_predictors(self,fname):
        predictors = open(fname).readline().strip().split("\t")  # read avaliable features
        return predictors

    # Filter list of available predictors according to the filter rule
    # Rule should be a valid re expression
    def filter_predictors(self, rule, keep):
        try:
            self.predictors
        except:
            raise Exception("Please read data first")
        expression = re.compile(rule)
        self.predictors = sorted([p for p in self.predictors if (expression.search(p) != None) == keep ])
        logging.info("Using following predictors: "+" ".join(self.predictors))


    def read_file(self,inp_file):
        # Pandas reads file much faster if it knows datatypes
        # Actually, we know that all predictors are float except chr (str)
        # and contact start/end/dist, which is int
        dtypes = {"chr": str, "contact_st": int64, "contact_en": int64,
                  "contact_dist": int64, "contact_count": float32}
        header = open(inp_file, "r").readline().strip().split("\t")

        fixed_dtypes_count = len(dtypes.keys())
        for i in header[fixed_dtypes_count:]:
            dtypes[i] = float32

        logging.getLogger(__name__).info("Reading file "+inp_file)
        input_data = pd.read_csv(inp_file, delimiter="\t", dtype=dtypes,
                                 header=0, names=header)
        # print(input_data.keys())
        # print(len(input_data.keys()))
        input_data.fillna(value=0, inplace=True) # Filling N/A values TODO check why N/A appear
        return input_data

    # Read available predictors, drop non-predictors
    # Fill self.predictors
    def read_data_predictors(self,inp_file):
        header = self.get_avaliable_predictors(inp_file)
        self.predictors = [h for h in header if not h in self.constant_nonpredictors]
        # print(self.predictors)
        self.input_file = inp_file