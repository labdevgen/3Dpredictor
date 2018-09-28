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

def overweight_loops(contacts,predictors):
    pass
