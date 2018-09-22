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
from shared import Interval, str2hash
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

    logging.getLogger(__name__).info("Reading data")
    input_data = pd.read_csv(inp_file,delimiter="\t",dtype=dtypes,header=1,names=header)
    input_data.fillna(value=0, inplace=True)

    if logFunc:
        results = input_data["contact_count"].apply(math.log)
    else:
        results = input_data["contact_count"]

    input_data.drop(["contact_count","chr", "contact_st", "contact_en"],axis = 1,inplace=True)
    input_data = input_data[[i for i in input_data.columns.get_values() if i in featuresSubset]] #preserve order of columns

    logging.getLogger(__name__).info("Fitting model")
    model.fit(input_data,results)
    return model

def validate(model,inp_file,featuresSubset,prefix):
    input_data = pd.read_csv(inp_file,delimiter="\t")
    input_data.fillna(value=0, inplace=True)
    if logFunc:
        results = input_data["contact_count"].apply(math.log)
    else:
        results = input_data["contact_count"]

    input_data.drop(["contact_count","chr", "contact_st", "contact_en"],axis = 1,inplace=True)
    input_data = input_data[[i for i in input_data.columns.get_values() if i in featuresSubset]] #preserve order of columns

    predicted = model.predict(input_data)
    r2 = sklearn.metrics.r2_score(predicted, results)

    #Plot r2
    subset = max(len(predicted) // 5000, 1)
    if "contact_dist" in input_data.columns.get_values():
        plt.scatter(predicted[::subset],results[::subset],
                c=input_data["contact_dist"].values[::subset])
    else:
        plt.scatter(predicted[::subset],results[::subset])
    plt.title("Test: Predicted vs real, r^2 score = "+str(r2)+"\n"+prefix)
    plt.xlabel("Log(Predicted Contact)")
    plt.ylabel("Log(Real Contact)")
    print ("Saveing file " + inp_file+".scatter"+str2hash(prefix)+".png")
    plt.savefig(inp_file+".scatter"+str2hash(prefix)+".png",dpi=300)
    plt.show()
    plt.clf()

    try:
        plt.plot(range(len( model.feature_importances_)), model.feature_importances_, marker = "o")
        print (input_data.columns.names)
        plt.xticks(range(len( model.feature_importances_)), input_data.columns.get_values(), rotation='vertical')
        plt.xlabel("Feature")
        plt.ylabel("Importance")
        plt.title("Feature importances")
        plt.savefig(validation_file+".featImportance."+str(logFunc)+str2hash(prefix)+".png",dpi=300)
        plt.show()
        plt.clf()
    except:
        logging.getLogger(__name__).warning("Feature importances is not avaliable for the model "+str(model.__class__.__name__))

    mp = MatrixPlotter()
    input_data = pd.read_csv(inp_file,delimiter="\t") #read again to get chr and contacts count
    predicted_data = input_data.copy(deep=True)
    input_data["contact_count"] = results #just to account for log that might have been applied
    predicted_data["contact_count"] = predicted
    mp.set_data(input_data)
    mp.set_control(predicted_data)
    pickle.dump(mp, open(inp_file + '.mpobject.dump', 'wb') )
    matrix = mp.getMatrix4plot(Interval(input_data["chr"].iloc[0],
                                        min(input_data["contact_st"].values),
                                        max(input_data["contact_en"].values)))
    if not logFunc:
        matrix = np.log(matrix)

    tick_pos, tick_labels = mp.get_bins_strart_labels(maxTicksNumber=15)
    plt.xticks(tick_pos,tick_labels,rotation=45)
    plt.imshow(matrix,cmap="OrRd")
    plt.title(prefix)
    plt.imsave(inp_file+".matrix."+str2hash(prefix)+".png",matrix,cmap="OrRd",dpi=600)
    plt.show()

def get_avaliable_predictors(file):
    predictors = open(file).readline().strip().split("\t")  # read avaliable features
    return predictors

def trainAndValidate(lm, training_file=None, validation_file=None, featuresSubset="all",order="keep",rewriteModel=False):
    if training_file != None:
        features = get_avaliable_predictors(training_file) #read avaliable features
    elif validation_file != None:
        features = get_avaliable_predictors(validation_file)  # read avaliable features
    else:
        raise Exception("Please provide either training or validation file")

    prefix = str(lm.__class__.__name__)
    if featuresSubset == "all":
        prefix += order+"."+featuresSubset+str(len(features))
        featuresSubset = features
    else:
        found_features = [f in features for f in featuresSubset]
        if not all(found_features):
            logging.error("Some requested features not found in dataset")
            logging.error("\t".join(np.array(featuresSubset)[np.logical_not(found_features)]))
            raise Exception("Feature not found")
        prefix += order + "." + ".".join(featuresSubset)
        if order == "keep":
            pass
        if order == "drop":
            featuresSubset = [f for f in features if not (f in featuresSubset)]

    if "contact_count" in features: #contacts count is always dropped
        del features[features.index("contact_count")]
    if "chr" in features: #chr is always dropped
        del features[features.index("chr")]
    if "contact_st" in features:
        del features[features.index("contact_st")]
    if "contact_en" in features:
        del features[features.index("contact_en")]

    prefix += ".log"+str(logFunc)

    model_file = training_file + ".model"+str2hash(prefix)+".dump"

    logging.getLogger(__name__).info("Using following features: "+" ".join(featuresSubset))

    if os.path.exists(model_file) or rewriteModel:
        model = pickle.load(open(model_file,"rb"))
    else:
        model = train(lm,training_file,featuresSubset=featuresSubset)
        pickle.dump(model,open(model_file,"wb"))

    logging.getLogger(__name__).info("Validating model")
    validate(model,validation_file,featuresSubset=featuresSubset,prefix=prefix)

#lm = linear_model.LinearRegression()
#lm = linear_model.Lasso(alpha=0.2)
#lm = linear_model.SGDRegressor()
#lm = linear_model.TheilSenRegressor()
#lm = linear_model.HuberRegressor()
#lm = ensemble.AdaBoostRegressor()
#lm = ensemble.RandomForestRegressor()
lm=ensemble.GradientBoostingRegressor()

#training_file = "training.RandOnChr11000000.50001.1000000.5000.100000.txt"
#training_file = "training.RandOnChr13000000.50001.3000000.10000.500000.txt"
#training_file = "Data/2018-08-20-trainingSmall.RandOnChr1.20000.contacts.3000000.50001.500000.25000.txt"
training_file = "2018-09-17-trainingOrient.RandOnChr1.20000.contacts.3000000.50001.500000.25000.txt"

validation_files = ["Interval_chr1_100000000_110000000validatingOrient.20000.contacts.3000000.50001.500000.25000.txt",
                    "Interval_chr2_47900000_53900000validatingOrient.20000.contacts.3000000.50001.500000.25000.txt",
                    "Interval_chr2_85000000_92500000validatingOrient.20000.contacts.3000000.50001.500000.25000.txt",
                    "Interval_chr10_59000000_62000000validatingOrient.20000.contacts.3000000.50001.500000.25000.txt"]
# keep = ["all",["CTCF_W","contact_dist"],
#         ["CTCF_W","contact_dist","CTCF_L","CTCF_R","CTCF_LDist_0","CTCF_LDist_2","CTCF_RDist_0","CTCF_RDist_2"]]
predictors = get_avaliable_predictors(training_file)
keep = []
#all except E1
keep+=[[p for p in predictors if p.find("E1")==-1 and p.find("contact")==-1]]

# #Contact distance + all non-CTCF chipSeqs
# keep += ["contact_dist"] + [p for p in predictors if p.find("CTCF")==-1 and p.find("E1")==-1]
#
# #All except non-CTCF chipSeq
# keep = [["contact_dist"] + [p for p in predictors if p.find("CTCF")!=-1 or p.find("E1")!=-1]]
#
# #Contacts_dist + CTCF
# keep += [["contact_dist"] + [p for p in predictors if p.find("CTCF")!=-1]]

for featuresSubset in keep:
    for validation_file in validation_files:
        trainAndValidate(lm=lm,
                     training_file = training_file, validation_file = validation_file,
                     featuresSubset = featuresSubset, order= "keep")

# for drops in [[],
#               ["E1_L","E1_R"],
#               ["CTCF_L","CTCF_R"],
#               ["CTCF_W"]]:
#     drop += drops
#     prefix = ".".join(drop[1:])