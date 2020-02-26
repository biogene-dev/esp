# -*- coding: utf-8 -*-
"""
Created on Thu Mar 21 18:00:03 2019

@author: ash
"""

import serial

f = open("impact.txt", "w") 

ser = serial.Serial('com3', baudrate=2000000, bytesize=8, parity='N', stopbits=1, timeout=100, xonxoff=1, rtscts=0)
ser.set_buffer_size(rx_size = 1280000, tx_size = 1280000)

for i in range (1000):
    data = ser.readline().decode("utf-8")
try:
    while 1 :
        #ser.reset_input_buffer()
        data = ser.read(500).decode()
    #    print(data)
        f.write(data)
    #    f.flush()
except (KeyboardInterrupt, SystemExit):
    f.close()
    ser.close()
    print ("je quitte")

