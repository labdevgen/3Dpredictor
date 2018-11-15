from Predictor import Predictor, ones_like, equal
import xgboost
import pandas as pd
import logging
from sklearn.metrics import f1_score, recall_score, precision_score
from Weight_funcs_modul import decorate_overweight_loops_for_classifier

class LoopsPredictor(Predictor):

    def __init__(self):
        self.classes_ratio = None
        super(LoopsPredictor, self).__init__()

    def toXMLDict(self):
        result = super(LoopsPredictor, self).toXMLDict()
        result["classes_ratio"] = str(self.classes_ratio)
        return result


    def train(self, alg = xgboost.XGBClassifier(n_estimators=100), #TODO,
              shortcut = "model", apply_log = False,
              dump = True, out_dir = "out/models/loops/",
              weightsFunc = ones_like,
              show_plot = True,
              classes_ratio = None,
              *args, **kwargs):

        self.classes_ratio = classes_ratio

        return super(LoopsPredictor, self).train(alg,shortcut,apply_log,dump,out_dir,weightsFunc,show_plot,
                                          *args,**kwargs)


    # Read file
    # Rename contact_count to isLoop
    def read_file(self,inp_file):
        input_data = super(LoopsPredictor,self).read_file(inp_file)
        input_data["contact_count"] = input_data["IsLoop"].astype(int)
        input_data.drop(["IsLoop"],axis="columns",inplace=True)
        logging.info("Equalizing classes")
        try:
            del self.predictors[self.predictors.index("IsLoop")]
        except:
            logging.getLogger(__name__).warning("IsLoop not found in predictors header")
        return input_data

    # Loops are rare, so drop some of non-looping contacts so that N of loops
    # is ~same as N of not-loops
    def equalize_classes(self, input_data):
        try:
            self.classes_ratio
        except:
            raise Exception("Please set classes ratio")

        ratio = input_data["contact_count"].value_counts(normalize=True)
        logging.getLogger(__name__).info(" Classes ratio: \n"+str(ratio))
        logging.getLogger(__name__).info(" Classes ratio: \n"+ \
                                         str(input_data["contact_count"].value_counts(normalize=False)))

        classes_counts = input_data["contact_count"].value_counts(normalize=False)
        assert classes_counts[0] >= classes_counts[1]
        to_keep = (self.classes_ratio[0]/self.classes_ratio[1])*classes_counts[1]

        loops = input_data.query("contact_count == 1")
        nonloops = input_data.query("contact_count == 0")
        nonloops = nonloops.sample(frac=min(1,to_keep/len(nonloops)))
        input_data = pd.concat([loops,nonloops])

        ratio = input_data["contact_count"].value_counts(normalize=True)
        logging.getLogger(__name__).info(" Classes ratio: \n"+str(ratio))
        logging.getLogger(__name__).info(" Classes ratio: \n" + \
                                         str(input_data["contact_count"].value_counts(normalize=False)))
        return input_data

    def check_loops_avaliable(self):
        if not "IsLoop" in self.predictors:
            raise Exception ("IsLoop fieled not found in predictors \n"+str(self.predictors))

    def read_data_predictors(self,inp_file):
        super(LoopsPredictor, self).read_data_predictors(inp_file)
        self.check_loops_avaliable()

    def f1score(self,validation_data,predicted,out_dir,**kwargs):
        f1 = f1_score(y_true=validation_data["contact_count"].values,
                                            y_pred=predicted)
        p = precision_score(y_true=validation_data["contact_count"].values,
                                            y_pred=predicted)
        r = recall_score (y_true=validation_data["contact_count"].values,
                                            y_pred=predicted)
        logging.getLogger(__name__).info("F1 score: "+ str(f1))
        logging.getLogger(__name__).info("Percision score: "+ str(p))
        logging.getLogger(__name__).info("Recall score: " + str(r))
        logging.getLogger(__name__).info("N Positive: "+ str(sum(predicted)))

    def validate(self, validation_file,
                 out_dir="out/pics/",
                 validators=None,
                 transformation=equal,
                 **kwargs):
        super(LoopsPredictor, self).validate(validation_file,out_dir,[self.f1score],transformation,**kwargs)

logging.getLogger(__name__).setLevel(logging.INFO)
train_file = "out/Hepat/2018-09-25-training.RandOnchr2contacts.gz.1000000.50001.500000.25000.txt"
validation_file = "out/Hepat/2018-09-25-training.RandOnchr10contacts.gz.1000000.50001.500000.25000.txt"

#for z in [1,2,5,10,20,50]:
for coeff in [2]:#,5,10,20,30,35]:
    logging.getLogger(__name__).info(str(coeff))
    for z in [50]:
        logging.getLogger(__name__).info("1:"+str(z))
        loops_predictor = LoopsPredictor()
        loops_predictor.read_data_predictors(train_file)
        trained_predictor = loops_predictor.train(shortcut="TestLoops", apply_log=False,
                                            classes_ratio = {1:1,0:z},
                                            show_plot=True,
                                            weightsFunc=decorate_overweight_loops_for_classifier(coeff=coeff))
        trained_predictor.out_dir = "out/models/loops/"
        trained_predictor.validate(validation_file,show_plot=True)
