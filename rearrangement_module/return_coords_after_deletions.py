import pandas as pd
from shared import Interval
import pickle
import swifter
"Z:/scratch/201903101031polina/3DPredictor/out/pics/scc/model4227386657.validation..equal.h=2.NPC.Interval_chr10_3100000_130575000validatingOrient.contacts.gz.8.1500000.50001.809529.25000.txt.scc"
"model4931696814.validation..equal.h=2.deletion_unc5b.3062033198.scc"
"model4931696814.validation..equal.h=2.deletion_unc5b.1869607811.scc"
"model4931696814.validation..equal.h=2.deletion_unc5b.Interval_chr10_3100000_130575000validatingOrient.contacts.gz.6.1500000.50001.809529.25000.txt.scc"
"model4931696814.validation..equal.h=2.deletion_unc5b.3972806291.scc"
"model4931696814.validation..equal.h=2.deletion_unc5b.1869607811.scc"
"model4931696814.validation..equal.h=2.deletion_unc5b.Interval_chr10_3100000_130575000validatingOrient.contacts.gz.7.1500000.50001.809529.25000.txt.scc"
"chr1.5000.model2224079568.validation..equal.Dactyly.chr1validatingOrient.contacts.gz.False.11.1500000.10001.MusDelB..5.all_cont.5000.txt.scc"
data=pd.read_csv("MusDelB/5KB/chr1.5000.model2224079568.validation..equal.Dactyly.chr1validatingOrient.contacts.gz.False.11.1500000.10001.MusDelB..5.all_cont.5000.txt.scc", sep=" ")
print(data.keys())
# interval = Interval("chr10", 60755595, 60761099)
interval = Interval("chr1", 76392403, 78064264)
print("pandas applying")
data["contact_st"] = data["contact_st"].apply(lambda x: x if x < interval.start else x+interval.len)
data["contact_en"] = data["contact_en"].apply(lambda x: x if x < interval.start else x + interval.len)
data = data[["contact_st", "contact_en", "contact_count"]]
data.rename(columns={"contact_count":"count", "contact_st":"st", "contact_en":"end"}, inplace=True)
# print(data.keys())
# print("done")
# print("dumping_data")
# with open("del_data_returned_coords.pickle", 'wb') as f:
#     pickle.dump(data, f)
# print("writing_data")
data.to_csv("MusDelB/5KB/MusDelB_predicted.txt", sep="\t", index=False)
print("done")