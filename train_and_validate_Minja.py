import logging
from Predictor import Predictor

logging.basicConfig(format='%(asctime)s %(message)s', datefmt='%I:%M:%S', level=logging.DEBUG)
training_file = "out/2018-09-22-trainingOrient.RandOnChr1.contacts.gz.1000000.50001.500000.25000.txt"
validation_files = [
    "out/Interval_chr2_36000000_41000000validatingOrient.contacts.gz.1000000.50001.500000.25000.txt",
    "out/Interval_chr2_47900000_53900000validatingOrient.contacts.gz.1000000.50001.500000.25000.txt",
    "out/Interval_chr2_85000000_92500000validatingOrient.contacts.gz.1000000.50001.500000.25000.txt"
            ]
predictor = Predictor()
predictor.read_data(training_file)
predictor.train(shortcut="all")
for validation_file in validation_files:
    predictor.validate(validation_file)
