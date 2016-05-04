# Author: Hoang NT
# Created: 2016/05/04
# Measure the convergence rate of FDM step09.py
# compare with the exact solution from BEM step02.py

import numpy as np
import math
from matplotlib import pyplot, cm
from mpl_toolkits.mplot3d import Axes3D

def fem(nx=11, ny=11, nit=100, draw=False) :
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
      if draw:
        ax = fig.gca(projection='3d')
        ax.plot_surface(X, Y, p, rstride=1, cstride=1, cmap=cm.coolwarm)
        ax.set_zlim3d(0,1)
        ax.view_init(elev=50., azim=-130.)
        pyplot.pause(0.05)
        pyplot.clf()
  if draw:
    pyplot.show()
  return p

def bem(n=64, inf=101, draw=False) :
  pi = math.pi
  x = np.zeros(n)
  y = np.zeros(n)
  u = np.zeros(n)
  bc = np.zeros(n)
  for i in range(0,n):
    if i < (n/4) :
      x[i] = i*8./n
      y[i] = 0
      u[i] = 0
      bc[i] = 1
    elif i < n/2 + 1 :
      x[i] = 2
      y[i] = (i-n/4)*4./n
      u[i] = (i-n/4)*4./n
      bc[i] = 0
    elif i < 3*n/4:
      x[i] = (3*n/4-i)*8./n
      y[i] = 1
      u[i] = 0
      bc[i] = 1
    else :
      x[i] = 0
      y[i] = (n-i)*4./n
      u[i] = 0
      bc[i] = 0
  ip1 = np.arange(n)
  ip1 += 1
  ip1[n-1] = 0

  xm = 0.5 * (x+x[ip1])
  ym = 0.5 * (y+y[ip1])
  dx = x[ip1] - x
  dy = y[ip1] - y
  d = np.zeros(n)
  
  for i in range(0,n):
    d[i] = math.sqrt(dx[i]*dx[i] + dy[i]*dy[i])

  G = np.zeros((n,n))
  H = np.zeros((n,n))
  for i in range(0,n):
    for j in range(0,n):
      if i != j:
        rx = xm[i] - xm[j]
        ry = ym[i] - ym[j]
        r = math.sqrt(rx*rx + ry*ry)
        G[i,j] = -math.log(r) * d[j]/2/pi
        H[i,j] = (rx*dy[j] - ry*dx[j])/r/r/2/pi
    G[i,i] = d[i] * (1-math.log(d[i]/2))/2/pi
    H[i,i] = 0.5

  for i in range(0,n):
    if bc[i] == 1:
      tmp = G[:,i]
      G[:,i] = -H[:,i]
      H[:,i] = -tmp

  b = H.dot(u)

  un = np.linalg.solve(G,b)

  for i in range(0,n):
    if bc[i] == 1:
      tmp = u[i]
      u[i] = un[i]
      un[i] = tmp

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


def scatter_plot(x,y,p,q):
  fig = pyplot.figure(figsize=(11,7), dpi=100)
  ax = fig.gca(projection='3d')
  ax.scatter(x,y,p,c='b')
  ax.scatter(x,y,q,c='r')
  ax.set_zlim3d(0,1)
  ax.view_init(elev=40., azim=-130.)
  pyplot.show()

def get_border(a):
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
