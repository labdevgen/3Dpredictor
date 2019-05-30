**3DPredictor**

**How to use the code**

The project contains 2 major modules:

**1. Data Generation module**

The module _GenerateData_K562.py_ loads some external files 
(i.e. file with known contact frequencies, ChipSeq, E1 and etc.)
and builds a dataset with predictor values for each contact.

The arhitecture of the data generation is following:

1. Reader object is responsible for parsing of specific data (i.e. ChiP-seq) and store it in pandas dataframe format.

2. Predictor_generator object uses reader object as a proxy to access data and based on epigenetic data provided by reader and coordinates of pair of loci generates specific predictors (in other words, performs parametrization)

3. DataGenerator object uses one contact_reader object and list of predictor_generator objects (each linked to its reader object)  to generate predictors for each contacts accessible by contact_reader object.

The _GenerateData_K562.py_ simply wraps DataGenerators object, providing specific file names.

Basic usage:
Set variables:

    '''python
    params.window_size = 25000 #region around contact to be binned for predictors
    params.small_window_size = 12500 #region  around contact ancors to be considered as cis
    params.mindist = 50001 #minimum distance between contacting regions
    params.maxdist = params.window_size #max distance between contacting regions
    params.maxdist = 3000000
    params.binsize = 20000 #when binning regions with predictors, use this binsize
    params.sample_size = 500000 #how many contacts write to file
    params.conttype = "contacts"
    '''

as well as filenames and genomic intervals of interest 
(see the code), and run. Data files sholud be located in 
 
    '''python
    input_folder = "set/path/to/this/variable"
The data files currently used could be downloaded from 
http://genedev.bionet.nsc.ru/site/hic_out/3DPredictor/

Note that predictors generation takes ~3h for 500 000 contacts.
One may change parallelization options by tweak code in DataGenerator.py:


    '''python
    n_cpus = multiprocessing.cpu_count()
    '''


There is an example of generating data for K562 cells

**Rearrangements**

If you want to predict contacts after rearrangement you should generate predictors for contacts with rearrangement.You should apply special rearrangement function for EVERY 'reader'. There are functions for duplication, deletion and inversion. For example:

    '''python
    params.contacts_reader.delete_region(Interval("chr22", 16064000, 16075000))
    '''

Output of this module are file with features (predictors) for all contacts and xml file with definition of dataset. Output format: contact_st--contact_end--contact_dist--other features...

**2. Training and validation module**

Run module _train_and_validate_K562.py_ to train model.
Set up following variables before running:

training_files - files with data for training

validation_files - a list of files with data for validation

contact type - contacts or OE values:
 
    '''python
    for contact_type,apply_log in zip(["contacts"],[True]): 
    '''

keep - a list with predictors to use. Each entery should 
contain a list of predictors, or a single string _"all"_
(for using all avaliable predictors). 
E.g. with
 
    '''python
    keep = ["CTCF_W","contact_dist"] 
    '''

 only CTCF sites inbetween contacting loci and distance between them will be used.
 
learning algorithm - you can change this in the _Predictor.py_ module:
 
    '''python
    def train(self,alg = xgboost.XGBRegressor(n_estimators=100,max_depth=9,subsample=0.7 
    '''

Also set folders name for validating results in _Predictor.py_ module

    '''python
    dump = True, out_dir = "/mnt/scratch/ws/psbelokopytova/201905031108polinaB/3DPredictor/out/models/" 
    out_dir = "/mnt/scratch/ws/psbelokopytova/201905031108polinaB/3DPredictor/out/pics/",
    '''
 Please note that fitting model is time-consuming so fitted model is saved to the file with the name 
 representing model parameters (predictors and algorithm). 
 It is automatically determined wheather such file exists and
 model can be simply loaded without fitting again 
 (to change this behaviour pass _rewriteModel=True_ when calling 
 _trainAndValidate_ function)
 
Choose validators to estimate accuracy of the algorithm. There are 2 main validators: SCC metric and plot_juicebox module for visual assessment of prediction
Output files:

| file name | definition |
| --------- | ---------------- |
| file.xml | model definition |
| featureImportances.png | histogram of feature importances |
| featureImportances.txt | list of feature importances scores |
| file.ssc.out | file with standard metrics and SCC for prediction |
| file.scc | pre-file for SCC counting which looks like contact_st--contact_end--real_cont_count--predicted_cont_count |
| data.hic | predicted heatmap |
| control.hic | heatmap with experimental data |

**Rearrangements**

For heatmap prediction after rearrangement use corresponding validating file and choose transformation option (now it works only for deletion):

    '''python
    trained_predictor.validate(validation_file, show_plot = False,cell_type=cell_type,
                            #use transformation option if you want to return coordinates after rearrangement
                            #                            transformation=
                            # [decorate_return_coordinates_after_deletion(return_coordinates_after_deletion, interval=deletion)], 
    '''
