**How to migrate to gitLab**

1.(optional, for PyCharm users):

--Commit and push current version to GitHub.
--Install GitLab Project Plugin (https://plugins.jetbrains.com/plugin/7975-gitlab-projects)
--Add GitLab repo in remotes (VCS->Git->Remotes)
--Fetch, Pull, Push (all under VCS control)

2.To migrate with your local copy of GitHub repo, see 
https://stackoverflow.com/questions/20359936/import-an-existing-git-project-into-gitlab


**How to use the code**

The project contains 2 major modules:

**1. Data Generation module**

The module _GenerateData_new.py_ loads some external files 
(i.e. file with contacts frequencies, ChipSeq, E1 and etc.) 
and builds a dataset with predictor values for each contact
It mainly wraps  DataGenerators classes with specific file name,
so read comments and code in DataGenerators.py to get idea of
how it works.
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
http://genedev.bionet.nsc.ru/hic_out/3DPredircor

Note that predictors generation takes ~3h for 500 000 contacts.

**2. Training and validation module**

Run module train_and_validate2.py to train model.
Set up following variables before running:

lm - learning algorithm

training_file - file with data for training

validation_files - a list of files with data for validation

keep - a list with predictors to use. Each entery should 
contain a list of predictors, or a single string _"all"_
(for using all avaliable predictors). 
E.g. with
 
    '''python
    keep = ["CTCF_W","contact_dist"] 
    '''
 only CTCF sites inbetween contacting loci and distance between them will be used.
 Please note that fitting model is time-consuming so fitted model is saved to the file with the name 
 representing model parameters (predictors and algorithm). 
 It is automatically determined wheather such file exists and
 model can be simply loaded without fitting again 
 (to change this behaviour pass _rewriteModel=True_ when calling 
 _trainAndValidate_ function)