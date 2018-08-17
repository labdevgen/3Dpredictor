import pandas as pd
import sklearn
from sklearn.linear_model import LinearRegression
from sklearn import linear_model,ensemble
import matplotlib.pyplot as plt
import math

input_file = "out.txt"
input_data = pd.read_csv(input_file,delimiter="\t")

results = input_data["contact_count"].apply(math.log)
drop = ["chr","contact_count"]
#drop += ["contact.st","contact.en",
#                 "window_start","window_end","contact_relative_start",
#                 "contact_relative_end"]
input_data.drop(drop,axis = 1,inplace=True)

X_train,X_test,Y_train,Y_test = sklearn.model_selection.train_test_split(
    input_data,results,test_size=0.15,random_state=5)
#lm = linear_model.LinearRegression()
#lm = linear_model.Lasso(alpha=0.2)
#lm = linear_model.SGDRegressor()
#lm = linear_model.TheilSenRegressor()
#lm = linear_model.HuberRegressor()
#lm = ensemble.AdaBoostRegressor()
#lm = ensemble.RandomForestRegressor()
lm = ensemble.GradientBoostingRegressor()
lm.fit(X_train,Y_train)

#plt.plot(lm.coef_[:200])
#plt.show()
#plt.clf()

predicted = lm.predict(X_test)
r2 = sklearn.metrics.r2_score(predicted, Y_test)
plt.scatter(predicted,Y_test,c=(X_test["contact_en"]-X_test["contact_st"]).values)
plt.title("Test: Predicted vs real, r^2 score = "+str(r2)+"\nModel: "+str(lm.__class__))
plt.xlabel("Log(Predicted Contact)")
plt.ylabel("Log(Real Contact)")
plt.savefig(input_file+".scatter.png",dpi=300)
plt.show()
plt.clf()

print ("----")

predicted = lm.predict(X_train)
plt.scatter(predicted,Y_train)
plt.title("Train: Predicted vs real")
plt.show()
plt.clf()