# -*- coding: utf-8 -*-
"""
Created on Wed Sep 25 19:05:07 2019

@author: ash
"""
import math
# sphinx_gallery_thumbnail_number = 3
import matplotlib.pyplot as plt
import numpy as np
#import parametre

from tkinter import filedialog
from tkinter import *

file = open ("recordday", "r")
temps = []
raquette = []
balle =[]
for e in file :
    a = (e.split(","))
    temps.append(int (a[0]))
    raquette.append(int(a[2]))
    balle.append(int(a[1]))
file.close()

plt.figure(1)
plt.figure( figsize=(8, 6))
#plt.gcf().subplots_adjust(wspace = 0, hspace = 4)

plt.plot(balle,'k' )
plt.show()
plt.figure(1)
plt.figure( figsize=(8, 6))
#plt.gcf().subplots_adjust(wspace = 0, hspace = 4)

plt.plot(raquette,'k' )
plt.show()