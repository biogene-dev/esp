# -*- coding: utf-8 -*-
"""
Created on Thu Mar 21 18:00:03 2019

@author: ash
"""

import serial

f = open("myfile.txt", "w") 

ser = serial.Serial('com4', 115200)

ser.reset_input_buffer()

while 1 :
    #ser.reset_input_buffer()
    data = ser.readline().decode("utf-8")
    f.write(data)
    f.flush()

