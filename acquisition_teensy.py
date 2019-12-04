# -*- coding: utf-8 -*-
"""
Created on Wed Sep 25 17:41:08 2019

@author: ash
"""
import serial
import time
import math
# sphinx_gallery_thumbnail_number = 3
import matplotlib.pyplot as plt

ser = serial.Serial('com4', 115200)
file = open ("impact.txt", "w")
log = []
for i in range (1000):
        ser.readline().decode("utf-8")
temp_acquisition = 200
time_start = time.time()

for i in range (1000):
        ser.readline().decode("utf-8")
        
while time.time() < time_start+temp_acquisition :
    data = ser.readline().decode("utf-8")
    
    log.append(int(data[:-1].split(",")[1]))
    file.write(data[:-1])
    file.flush()
ser.close()
print (log)
plt.figure(1)
plt.figure( figsize=(8, 6))
plt.plot( log,'k' )
plt.show()


