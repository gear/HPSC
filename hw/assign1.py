# Author: Hoang NT
# Created: 2016/05/04
# Measure the convergence rate of FDM step09.py
# compare with the exact solution from BEM step02.py

import numpy as np
from matplotlib import pyplot, cm
from mpl_toolkits.mplot3d import Axes3D

def fmd(nx=11, ny=11, nit=100) :
  dx = 2./(nx-1)
  dy = 1./(ny-1)
  x = np.linspace(0,2,nx)
  y = np.linspace(0,1,ny)
  X, Y = np.meshgrid(x,y)
  p = np.zeros((ny,nx))
  fig = pyplot.figure(figsize=(11,7), dpi=100)
  
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
      ax = fig.gca(projection='3d')
      ax.plot_surface(X, Y, p, rstride=1, cstride=1, cmap=cm.coolwarm)
      ax.set_zlim3d(0,1)
      ax.view_init(elev=50., azim=-130.)
      pyplot.pause(0.05)
      pyplot.clf()
  pyplot.show()


