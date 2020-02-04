# -*- coding: utf-8 -*-
"""
Created on Thu Mar 21 18:00:03 2019

@author: ash
"""

import serial
import time
tps = time.localtime()
file_name =str(tps.tm_year)+"_"+str(tps.tm_mon)+"_"+str(tps.tm_mday)+"_" + str(tps.tm_hour)+"_"+str(tps.tm_min)+ ".txt"
f = open(file_name, "w") 

ser = serial.Serial('com3', baudrate=2000000, bytesize=8, parity='N', stopbits=1, timeout=100, xonxoff=1, rtscts=0)
ser.set_buffer_size(rx_size = 1280000, tx_size = 1280000)

for i in range (1000):
    data = ser.readline().decode("utf-8")
print ("ACQUISITION EN COURS")
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
    print ("Fin acquisition")

print ("traitement des donnees")
import numpy as np

global temps
global raquette
global balle
global accelero
temps= np.array([])
balle= np.array([])
raquette= np.array([])
accelero= np.array([])

f = open(file_name, "r") 
(f.read(1000)) 
val = 0.0 
print ("Formatage des datas")
rating = 0
for line in f:
    rating +=1
    if line :
        try :
#            print (line)
            
            to_test =line.split(',')
            val1 = float(to_test[0])
            val2 = float (to_test[2]) * -1
            val3 = float (to_test[3])
            val4 = float (to_test[4])
#            print ("val =",val)
        except :
            pass
            print("bad")
#            print (line)
        else :
#            print ("good" , val)
            # ----------------------------
            # conversion binaire, valeur reelle
    
            # fonction de transformation bit -> degre
            if ( len(to_test) == 6 ) and (val2 < 10000) and (val3 < 10000) and (val4 < 10000) :
                if int(to_test[0])< 9999999999:
                    temps = np.append(temps , int(int(val1)/1000))
                    raquette = np.append(raquette , (val2*360/(8192))*0.8 )
                    balle = np.append(balle , (val3*360/(8192))*0.8 )
                    accelero = np.append(accelero , val4)
#                val2=(val2*360/(8192))*0.8; # le 0.8 est du à la bande de 10# -> 90# du capteur
##            transAng=(val*360/((2^13)*0.8)); # le 0.8 est du à la bande de 10# -> 90# du capteur
#                if transAng < 10000 : a =np.append(a, transAng)
#                if int(to_test[0])< 9999999999 : Temps.append(int(int(to_test[0])/1000))

#            print (len(a))
        
print ("len de rec =" , len(balle))



