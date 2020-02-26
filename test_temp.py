# -*- coding: utf-8 -*-
"""
Created on Wed Oct 30 18:35:57 2019

@author: ash
"""
from pathlib import Path

p = Path('.')

file_to_open = p.absolute() / "calib_raquette.txt"
f= open(file_to_open,"r")
ex_val = 0
print(file_to_open)
print (f.readline(1000))
#

for line in f:
    val = int (line.split(",")[0])
    inter_temp = val-ex_val
    if inter_temp > 1100 : print (inter_temp)
    ex_val = val
    
f.close()