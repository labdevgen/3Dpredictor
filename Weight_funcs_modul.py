import numpy as np
from functools import partial

def ones_like(contacts, *args):
    return np.ones_like(contacts)

def array(contacts,*args):
    return np.array(contacts)

def abs_log(contacts,predictors):
    return np.abs(np.log(contacts))

def mult_abs_log(contacts,predictors,coeff):
    return np.abs(np.log(contacts))*coeff

def decorate_mult_abs_log(func,coeff):
    result = partial(func,coeff=coeff)
    result.__name__ = str(coeff) + func.__name__
    return result

def overweight_loops(contacts,predictors,coeff): #contacts is np.array of contact_count, predictors
    idx_loop= np.flatnonzero(predictors['IsLoop'])
    result = np.array(contacts)
    result[[idx_loop]] = result[[idx_loop]]*coeff
    return result

    # idx = predictors.columns.get_loc("isloop")
    # for i in range(len(contacts)):
    #     if predictors.iloc[i, idx] == 1:
    #         result[i] =
def decorate_overweight_loops(func,coeff):
    result = partial(func,coeff=coeff)
    result.__name__ = str(coeff) + func.__name__
    return result

#threshold - threshold counts of contacts, by default threshold is symmetric
#weigth = ((contacts/threshold)**power)*coeff
#abs -  boolean only, if True - only absolute value of contacts will be calculated
#piecing - boolean only, if True - weith function has piecewising form 
#asymmetric - -1,0,1
#if 1 - only 'big' (observed/expected > treshold) contacts will be reweighted
#if -1 - only 'small' (observed/expected < 1/treshold) contacts will be reweighted 
#if 0 - all contacts will are reweighted
def contactWeitherFunction(contacts,predictors, threshold,power,coeff,abs,piecing,asymmetric):
    log_con = np.log2(contacts)
    sign = np.sign(log_con)
    if asymmetric != 0: sign = np.trunc(asymmetric*sign+1)/2
    if threshold == 1 or threshold <= 0:
        print ('threshold = 1 or <= 0, returned contact weigths = 1')
        return contacts*0+1
    else:
        if abs == True: sign **= 2
        result = ((np.abs(log_con)/np.abs(np.log2(threshold)))**power)*sign
        if piecing == True: result = np.sign(np.trunc(result))
        result *= coeff
        nulls = np.abs(np.sign(np.trunc(result)))
        result *= nulls
        nulls = (nulls+1) % 2
        result += nulls
        return result

def decorateContactWeither(func, **kwargs):
    try:kwargs['threshold']
    except KeyError: threshold = 2.
    else: threshold = kwargs['threshold']
    try: kwargs['power']
    except KeyError: power = 1.
    else: power = kwargs['power']
    try: kwargs['coeff']
    except KeyError: coeff = 1.
    else: coeff = kwargs['coeff']
    try: kwargs['abs']
    except KeyError: abs = False
    else: abs = kwargs['abs']
    try: kwargs['piecing']
    except KeyError: piecing = False
    else:piecing = kwargs['piecing']
    try: kwargs['asymmetric']
    except KeyError: asymmetric = 0
    else: asymmetric = np.sign(kwargs['asymmetric'])
    name = ''
    for key in kwargs: name += ( key[:3]+str(kwargs[key]) )
    #result = contactWeitherFunction(contacts,treshold=treshold,power=power,coeff=coeff,abs=abs,piecing=piecing,asymmetric=asymmetric)
    result = partial(func,threshold=threshold,power=power,coeff=coeff,abs=abs,piecing=piecing,asymmetric=asymmetric)
    result.__name__ = name + func.__name__
    return result

