# -*- coding: utf-8 -*-
"""
Created on Sun Mar 31 19:04:32 2019

@author: ash
"""

#!python
#!/usr/bin/env python
import tables
file = tables.openFile('test.mat')
lon = file.root.lon[:]
lat = file.root.lat[:]
# Alternate syntax if the variable name is in a string
varname = 'lon'
lon = file.getNode('/' + varname)[:]