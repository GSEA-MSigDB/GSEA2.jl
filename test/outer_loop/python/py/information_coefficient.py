from numpy import asarray, exp, finfo, isnan, log, log2, nan, sign, sqrt, unique, power, e, pi


#from rpy2.robjects.numpy2ri import numpy2ri
# import rpy2.robjects.numpy2ri as numpy2ri
from rpy2.robjects.packages import importr
import rpy2.robjects as ro
from rpy2.robjects import r, pandas2ri

import numpy as np
from scipy.stats import pearsonr

eps = finfo(float).eps

# ro.conversion.py2ri = numpy2ri

pandas2ri.activate()

mass = importr("MASS")
    
def information_coefficient(x, y, n_grid=24, delta = 1.0):


    try:
        pearson_correlation = pearsonr(x, y)[0]
        
    except BaseException as err:
        print('Exception={} x={} y={}'.format(err, x, y))

    #print('pearson={}'.format(pearson_correlation))
    
    if isnan(pearson_correlation) or unique(x).size == 1 or unique(y).size == 1:

        return nan

    else:

        pearson_correlation_abs = abs(pearson_correlation)

       # xr = pandas2ri.py2ri(x)
       # yr = pandas2ri.py2ri(y)

        #print('IC: sum x:{} sum y:{}'.format(np.sum(x), np.sum(y)))
        
        #print('bandwidth x: {}'.format(mass.bcv(x)[0]))
        bandwidth_x = delta * mass.bcv(x)[0] * (1 - pearson_correlation_abs * 0.75)
        #bandwidth_x = delta * mass.bcv(asarray(list(x)))[0] * (1 - pearson_correlation_abs * 0.75)
        
        #print('bandwidth y: {}'.format(mass.bcv(y)[0]))
        bandwidth_y = delta * mass.bcv(y)[0] * (1 - pearson_correlation_abs * 0.75)
        #bandwidth_y = delta * mass.bcv(asarray(list(y)))[0] * (1 - pearson_correlation_abs * 0.75)
        
        #print('bandwidth... done')
        
        fxy = (asarray(mass.kde2d(x, y, asarray((bandwidth_x, bandwidth_y)), n=asarray((n_grid,)))[2]) + eps)
        #fxy = (asarray(mass.kde2d(asarray(list(x)), asarray(list(y)),
        #                              asarray((bandwidth_x, bandwidth_y)), n=asarray((n_grid,)))[2]) + eps)        

        #dx = (x.max() - x.min()) / (n_grid - 1)

        #dy = (y.max() - y.min()) / (n_grid - 1)

        pxy = fxy / (fxy.sum())

        px = pxy.sum(axis=1)
        px = px/px.sum()

        py = pxy.sum(axis=0)
        py = py/py.sum()

        #mi = (pxy * log2(pxy / (asarray((px,)).T * asarray((py,))))).sum()

        hxy = - (pxy * log2(pxy)).sum() 
        hx = -(px * log2(px)).sum() 
        hy = -(py * log2(py)).sum() 
        mi = hx + hy - hxy
        if mi < 0:
            mi = 0
            
        #print('MI={}'.format(mi))

        # The mutual information is normalized as in  Linfoot's Information Coefficient (https://www.sciencedirect.com/science/article/pii/S001999585790116X)
        # using the mutual information for a Gaussian distribution (see e.g. Example 8.5.1 in Elements of Information Theory 2nd ed - T. Cover, J. Thomas Wiley, 2006)

        IC = sign(pearson_correlation) * sqrt(1 - power(2.0, -2.0 * mi))

        return IC
             
