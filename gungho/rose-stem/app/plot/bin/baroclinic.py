#!/usr/bin/env python2.7
# -*- coding: utf-8 -*-
##############################################################################
# (C) Crown copyright 2017 Met Office. All rights reserved.
# For further details please refer to the file LICENCE which you should have
# received as part of this distribution.
##############################################################################
'''
Python script to plot xz slices along the y=0 and xy slices on a specified level 

This version takes nodal format output files and
interpolates onto a regular grid.

Filename hardcoded.

Levels are determined from the data

This version stitches together a directory of files
and extracts all levels so it can work in the serial
case where there is one file or the parallel case where
there is a file for each processor.

'''

import numpy as np
# Need to set a non-interactive backend for suites
import matplotlib
matplotlib.use('Agg')

import matplotlib.pyplot as plt
import matplotlib.cm as cm

from scipy.interpolate import griddata

import sys

from read_data import read_nodal_data

levels = None
data = None

       
def make_figure(plotpath, field, component, timestep, levels):

  val_col = 'c' + str(component)

   # get min and max of x,y data for plot axes

  min_lev = min(levels)

  deltaz = data.loc[data['level'] == (min_lev+1)][val_col].min() - data.loc[data['level'] == (min_lev)][val_col].min()
  xmin = data.loc[data['level'] == min_lev]['x'].min()
  xmax = data.loc[data['level'] == min_lev]['x'].max()
  ymin = data.loc[data['level'] == min_lev]['y'].min()
  ymax = data.loc[data['level'] == min_lev]['y'].max()

  zmin = min(levels)*1000.0
  zmax = max(levels)*1000.0

  r2d = 180.0/np.pi;
  nx,ny,nz = 360,180,len(levels)

  #create 2D plot
  x2d = np.linspace(xmin, xmax, nx)
  z2d = np.linspace(zmin, zmax, nz)
  y2d = np.linspace(ymin, ymax, ny)
  zi = np.zeros([ny,nx,len(levels)])
  xi, yi = np.meshgrid(x2d, y2d) 
  for p in xrange(len(levels)):
    p_data = data.loc[data['level'] == levels[p]]
    zi[:,:,p] = griddata((p_data['x'].values, p_data['y'].values), p_data[val_col].values, (xi, yi), method='linear')


  cc = np.linspace(np.amin(p_data[val_col].values),np.amax(p_data[val_col].values),13)

  # Pressure plot
  plotlevel = 1
  if field == 'theta':
    cc = np.linspace(220,330,12)
  elif field == 'exner':
    rd = 287.05
    p0 = 100000.0
    kappa = rd/1005.0
    zi = 0.01*zi**(1.0/kappa) * p0
    cc = np.linspace(900,1000,11)
    plotlevel = 0

  fig = plt.figure(figsize=(10,5))
  dz = zi[:,:,plotlevel]
  c_map = cm.summer
  cf = plt.contourf(xi *r2d, yi * r2d, dz, cc, cmap=c_map)
  plt.colorbar(cf,  cmap=c_map)
  cl = plt.contour(xi * r2d, yi * r2d, dz, cc, linewidths=0.5, colors='k')
  plt.xlabel('Longitude')
  plt.ylabel('Latitude')
  out_file_name = plotpath + "/" "slice_xy_" + field + "_" + timestep +  ".png"
  plt.savefig(out_file_name , bbox_inches='tight')
 
if __name__ == "__main__":

  try:
    config, datapath, fields, timesteps, plotpath = sys.argv[1:6]
  except ValueError:
    print("Usage: {0} <config> <datapath> <fields_list> <timestep_list> <plotpath>".format(sys.argv[0]))
    exit(1)

  # Split out the list of fields
  field_list = fields.split(':')

  # Split out the list of timesteps
  ts_list = timesteps.split(':')


  for ts in ts_list:
    for field in field_list:

      filestem =  datapath + "/" + config + "_nodal_" + field + "_" + ts + "*"      

      data = read_nodal_data(filestem, 1, 1)
      levels = data.level.unique()

      # Only try to plot if we found some files for this timestep
      if len(levels) > 0:
        make_figure(plotpath, field, 1, ts, levels)


