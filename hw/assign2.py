# Author: Hoang NT
# Created: 2016/05/28
# Measure the convergen rate of Laplace Problem numerical
# solution, compare with the exact solution in diferent
# configuration. Also implement higher order fdm scheme.

from __future__ import division
import numpy as np
import math
from matplotlib import pyplot, cm
from mpl_toolkits.mplot3d import Axes3D

# Import previous solution
import assign1 as a1

# Redefine exact solution
def exact(nx=32, ny=32, inf=101, draw=False) :
    '''
    nx : number of points for x-axis
    ny : number of points for y-axis
    inf : pseudo infinity for computation
    draw : plot the final result
    '''
    pi = math.pi
    x = np.zeros(nx)
    y = np
