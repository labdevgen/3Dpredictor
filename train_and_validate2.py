import sklearn
import math
import pandas as pd
from sklearn import linear_model,ensemble
import matplotlib.pyplot as plt
from matrix_plotter import MatrixPlotter
import pickle
import os
from numpy import int64,float32
import numpy as np
from shared import Interval
import logging

logging.basicConfig(format='%(asctime)s %(message)s', datefmt='%I:%M:%S', level=logging.DEBUG)
logFunc = True

def train(model,inp_file,featuresSubset):
    dtypes = {"chr":str,"contact_st":int64,"contact_en":int64,"contact_dist":int64,"contact_count":float32}
    header = open(inp_file,"r").readline().strip().split("\t")

    fixed_dtypes_count = len(dtypes.keys())
    for i in header[fixed_dtypes_count:]:
#        if i.startswith("CTCF") or i.startswith("E1"):
            dtypes[i] = float32

    logging.info("Reading data")
    input_data = pd.read_csv(inp_file,delimiter="\t",dtype=dtypes,header=1,names=header)

    if logFunc:
        results = input_data["contact_count"].apply(math.log)
    else:
        results = input_data["contact_count"]

    input_data.drop(["contact_count","chr"],axis = 1,inplace=True)
    input_data = input_data[[i for i in input_data.columns.get_values() if i in featuresSubset]] #preserve order of columns

    logging.info("Fitting model")
    model.fit(input_data,results)
    return model

def validate(model,inp_file,featuresSubset,prefix):
    input_data = pd.read_csv(inp_file,delimiter="\t")
    if logFunc:
        results = input_data["contact_count"].apply(math.log)
    else:
        results = input_data["contact_count"]

    input_data.drop(["contact_count","chr"],axis = 1,inplace=True)
    input_data = input_data[[i for i in input_data.columns.get_values() if i in featuresSubset]] #preserve order of columns

    predicted = model.predict(input_data)
    r2 = sklearn.metrics.r2_score(predicted, results)

    #plt.scatter(predicted,results,c=(input_data["contact_en"]-input_data["contact_st"]).values)
    plt.title("Test: Predicted vs real, r^2 score = "+str(r2)+"\nModel: "+str(model.__class__))
    plt.xlabel("Log(Predicted Contact)")
    plt.ylabel("Log(Real Contact)")
    plt.savefig(inp_file+".scatter"+prefix+".png",dpi=300)
    #plt.show()
    plt.clf()

    plt.plot(range(len( model.feature_importances_)), model.feature_importances_, marker = "o")
    print (input_data.columns.names)
    plt.xticks(range(len( model.feature_importances_)), input_data.columns.get_values(), rotation='vertical')
    plt.xlabel("Feature")
    plt.ylabel("Importance")
    plt.title("Feature importances")
    plt.savefig(validation_file+".featImportance."+str(logFunc)+prefix+".png",dpi=300)
    plt.show()
    plt.clf()


    mp = MatrixPlotter()
    input_data = pd.read_csv(inp_file,delimiter="\t") #read again to get chr and contacts count
    predicted_data = input_data.copy(deep=True)
    input_data["contact_count"] = results #just to account for log that might have been applied
    predicted_data["contact_count"] = predicted
    mp.set_data(input_data)
    mp.set_control(predicted_data)
    matrix = mp.getMatrix4plot(Interval(input_data["chr"].iloc[0],
                                        min(input_data["contact_st"].values),
                                        max(input_data["contact_en"].values)))
    if not logFunc:
        matrix = np.log(matrix)
    plt.imshow(matrix,cmap="OrRd")
    plt.show()
    plt.savefig(inp_file+".matrix."+prefix+".png",dpi=300)


def trainAndValidate(lm, training_file=None, validation_file=None, featuresSubset="all",order="keep",rewriteModel=False):
    if training_file != None:
        features = open(training_file).readline().strip().split("\t") #read avaliable features
    elif validation_file != None:
        features = open(validation_file).readline().strip().split("\t")  # read avaliable features
    else:
        raise Exception("Please provide either training or validation file")

    prefix = str(lm.__class__.__name__)
    if featuresSubset == "all":
        prefix += order+"."+featuresSubset+str(len(features))
        featuresSubset = features
    else:
        assert all([f in features for f in featuresSubset])
        prefix += order + "." + ".".join(featuresSubset)
        if order == "keep":
            pass
        if order == "drop":
            featuresSubset = [f for f in features if not (f in featuresSubset)]

    if "contact_count" in features: #contacts count is always dropped
        del features[features.index("contact_count")]
    if "chr" in features: #chr is always dropped
        del features[features.index("chr")]

    prefix += ".log"+str(logFunc)
    model_file = training_file + ".model"+prefix+".dump"

    logging.info("Using following features: "+" ".join(featuresSubset))

    if os.path.exists(model_file) or rewriteModel:
        model = pickle.load(open(model_file,"rb"))
    else:
        model = train(lm,training_file,featuresSubset=featuresSubset)
        pickle.dump(model,open(model_file,"wb"))

    logging.info("Validating model")
    validate(model,validation_file,featuresSubset=featuresSubset,prefix=prefix)

#lm = linear_model.LinearRegression()
#lm = linear_model.Lasso(alpha=0.2)
#lm = linear_model.SGDRegressor()
#lm = linear_model.TheilSenRegressor()
#lm = linear_model.HuberRegressor()
#lm = ensemble.AdaBoostRegressor()
#lm = ensemble.RandomForestRegressor()
#training_file = "training.RandOnChr11000000.50001.1000000.5000.100000.txt"
#training_file = "training.RandOnChr13000000.50001.3000000.10000.500000.txt"
training_file = "data/2018-08-20-trainingSmall.RandOnChr1.20000.contacts.3000000.50001.500000.25000.txt"
#validation_file = "validating.38Mb_58MbOnChr21000000.50001.1000000.5000.100000.txt"
#validation_file = "validating.38Mb_58MbOnChr23000000.50001.3000000.10000.500000.txt"
#validation_file = "validatingSmall.38Mb_58MbOnChr2.20000.contacts.3000000.50001.500000.25000.txt"
#validation_file = "Interval_chr10_59000000_62000000validatingSmall.20000.contacts.3000000.50001.500000.25000.txt"
#validation_file = "Interval_chr2_47900000_53900000validatingSmall.noChr2.20000.contacts.3000000.50001.500000.25000.txt"
validation_file = "data/Interval_chr2_47900000_53900000validatingSmall.20000.contacts.3000000.50001.500000.25000.txt"
#validation_file = "Interval_chr2_47900000_53900000validatingSmall.noCTCF.20000.contacts.3000000.50001.500000.25000.txt"

#drop = ["chr", "contact_count","contact_st","contact_en"]
keep = ["all",["CTCF_W","contact_dist"]]
keep = [["CTCF_W","contact_dist","CTCF_L","CTCF_R","CTCF_LDist_0","CTCF_LDist_4","CTCF_RDist_0","CTCF_RDist_4"]]
for featuresSubset in keep:
    trainAndValidate(lm=ensemble.GradientBoostingRegressor(),
                     training_file = training_file, validation_file = validation_file,
                     featuresSubset = featuresSubset, order= "keep")

# for drops in [[],
#               ["E1_L","E1_R"],
#               ["CTCF_L","CTCF_R"],
#               ["CTCF_W"]]:
#     drop += drops
#     prefix = ".".join(drop[1:])