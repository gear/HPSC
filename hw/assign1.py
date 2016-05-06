# Author: Hoang NT
# Created: 2016/05/04
# Measure the convergence rate of FDM step09.py
# compare with the exact solution from BEM step02.py

from __future__ import division
import numpy as np
import math
from matplotlib import pyplot, cm
from mpl_toolkits.mplot3d import Axes3D

# Discrete solution for 2D Laplace's solution
def fdm(nx=11, ny=11, nit=100, draw=False) :
  '''
  nx : x-axis spartial resolution
  ny : y-axis spartial resolution
  nit : number of time step
  draw : plot the result
  '''
  dx = 2./(nx-1)
  dy = 1./(ny-1)
  x = np.linspace(0,2,nx)
  y = np.linspace(0,1,ny)
  X, Y = np.meshgrid(x,y)
  p = np.zeros((ny,nx))
  
  for it in range(nit):
    for i in range(0,nx):
      p[0,i] = p[1,i]
      p[ny-1,i] = p[ny-2,i]
    for j in range(0,ny):
      p[j,0] = 0
      p[j,nx-1] = y[j]
    pn = p.copy()
    for i in range(1,nx-1):
      for j in range(1,ny-1):
        p[j,i] = (dy**2*(pn[j,i+1]+pn[j,i-1])+dx**2*(pn[j+1,i]+pn[j-1,i]))/(2*(dx**2+dy**2))
      if draw:
        fig = pyplot.figure(figsize=(11,7), dpi=100)
        ax = fig.gca(projection='3d')
        ax.plot_surface(X, Y, p, rstride=1, cstride=1, cmap=cm.coolwarm)
        ax.set_zlim3d(0,1)
        ax.view_init(elev=50., azim=-130.)
        pyplot.pause(0.05)
        pyplot.clf()
  if draw:
    pyplot.show()
  return p

# Exact boundary solution for 2D Laplace's equation
def exact(n=64, inf=101, draw=False) :
  '''
  n : number of boundary points
  int : pseudo infinity (101 is maximum for sinh func)
  draw : plot the final result
  '''
  pi = math.pi
  x = np.zeros(n)
  y = np.zeros(n)
  for i in range(0,n):
    if i < (n/4) :
      x[i] = i*8./n
      y[i] = 0
    elif i < n/2 + 1 :
      x[i] = 2
      y[i] = (i-n/4)*4./n
    elif i < 3*n/4:
      x[i] = (3*n/4-i)*8./n
      y[i] = 1
    else :
      x[i] = 0
      y[i] = (n-i)*4./n

  ux = np.zeros(n)
  for i in range(0,n):
    uxi = 0
    for j in range(1,inf,2):
      uxi += 1/(j*pi)**2/math.sinh(2*j*pi)*math.sinh(j*pi*x[i])*math.cos(j*pi*y[i])
      ux[i] = x[i] / 4 - 4*uxi
  
  if draw:
    fig = pyplot.figure(figsize=(11,7), dpi=100)
    ax = fig.gca(projection='3d')
    #ax.scatter(x,y,u,c='b')
    ax.scatter(x,y,ux,c='r')
    ax.set_zlim3d(0,1)
    ax.view_init(elev=40., azim=-130.)
    pyplot.show()
  return x,y,ux

# Plot 2 result p, q in the same figure
def scatter_plot(x,y,p,q):
  '''
  x : x-axis
  y : y-axis
  p : 1D boundary (exact) result
  q : 1D boundary (approximate) result
  '''
  fig = pyplot.figure(figsize=(11,7), dpi=100)
  ax = fig.gca(projection='3d')
  ax.scatter(x,y,p,c='b')
  ax.scatter(x,y,q,c='r')
  ax.set_zlim3d(0,1)
  ax.view_init(elev=40., azim=-130.)
  pyplot.show()

# Get border result from 2D fdm method
def get_border(a):
  '''
  a : 2D result from fdm method
  '''
  size = a.shape
  length = size[0]
  size = 2*(size[0]+size[1])
  ret = np.zeros(size)
  ret[0:length] = a[:,0]
  ret[length:2*length] = a[length-1,:]
  temp = a[:,length-1]
  ret[2*length:3*length] = temp[::-1]
  temp = a[0,:]
  ret[3*length:] = temp[::-1]
  return ret[::-1]

# Compute error between p and q
def error(p,q,tor=1e-6):
  '''
  p : 1D boudnary exact result
  q : 1D boundary approximate result
  tor : pseudo zero
  '''
  assert len(p) == len(q)
  length = len(p)
  error = 0.0
  for i in range(len(p)):
    d = p[i] - q[i] 
    r = 0.0
    if q[i] < tor:
      r = (d**2) / tor
    else:
      r = (d**2) / (p[i]**2)
    error += r
  return math.sqrt(error) / math.sqrt(length)

# Compute alternative relative square error
def error_rel(p,q):
  '''
  p : 1D boundary exact result
  q : 1D boundary approximate result
  '''
  assert len(p) == len(q)
  length = len(p)
  error = 0.0
  for i in range(len(p)):
    d = p[i] - q[i]
    r = 0.0
    if d != 0:
      r = (d**2) / ((p[i] + q[i])**2)
    error += r
  return math.sqrt(error) / math.sqrt(length)

# Experimental: Compute sum error
def error_rat(p,q,tor=1e-8):
  assert len(p) == len(q)
  sum1 = 0.0
  sum2 = 0.0
  for i in range(len(p)):
    d = p[i] - q[i]
    sum1 += d**2
    sum2 += p[i]**2
  if sum2 == 0:
    sum2 = tor
  return math.sqrt(sum1/sum2)
