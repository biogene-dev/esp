# -*- coding: utf-8 -*-
"""
Created on Thu Mar 21 20:38:57 2019

@author: ash
"""
import math
# sphinx_gallery_thumbnail_number = 3
import matplotlib.pyplot as plt
import numpy as np
import parametre

from tkinter import filedialog
from tkinter import *

data_bal=[]
data_raq=[]
data_3=[]
temps = []
data_int_1 =[]
data_int_2=[]

vit_raq=[]
vit_bal=[]

size_int =[]
v_raq=0
v_bal=0
'''
root = Tk()
root.filename =  filedialog.askopenfilename(initialdir = "C:/Users/ash/Documents/Thorbot",title = "Select file",filetypes = (("jpeg files","*.jpg"),("all files","*.*")))
print (root.filename)
fic = open(root.filename)
'''
fic = open("C:/Users/ash/Documents/Thorbot/1_1khz.txt")
t_zero =[0,0,0,0]
decouvert = False
t_impact = 0

data_moins_un=0

def E_cin (m , l , V_phi):
    return  (0.5 * m * math.pow(l,2) * math.pow(V_phi,2))

for i, line in enumerate (fic):
    #print (line)
    data= line.split(',')
    #print (len(data))
    if len(data)==4 :
        if (data_moins_un !=0) and (int (data[2]) > 7000) and (int(data[2])<data_moins_un) :
            t_impact  = i
            break
        data_moins_un = int(data[2])
            
#print ("temps impact = ", t_impact)
fic.seek(0)

t_zero=0
for i , line in enumerate(fic,1):
    #print (line)
    line = line[:-1]
    #print ("line",line)
    
    data= line.split(',')
    #print (len(data))
    if len(data)==4 :
        data = list(map(int, data))
        brut_bal = data[1]
        brut_raq = data[2]
        data_bal.append(brut_bal)
        data_raq.append(brut_raq)
        data_3.append(data[3])
        temps.append(data[0] -t_impact )
        if ( i - t_impact > -50 ) and ( i - t_impact < 50) :
            data_int_1.append(brut_bal)          
            data_int_2.append(brut_raq)
            size_int.append(i - t_impact)
        #print ("i-timpact" , i- t_impact)
        if ( i - t_impact > -50 and i - t_impact < -20 ):
            print ("raq", data[2] , " / ", t_zero[2], " / " ,  brut_raq - t_zero[2])
            v_raq += brut_raq - t_zero[2]
        if ( i - t_impact > 0 and i - t_impact < 200) :
            print ("bal", data[1] , " / ", t_zero[1], " / ", t_zero[1] - data[1])
            v_bal +=   t_zero[1] - data[1]
        t_zero =  data
print ("masse raquette ",parametre.m_tot_tennis , ", l raquette =",parametre.x_tot_tennis  )
print ("masse balle ",parametre.m_tot_pendule_balle , ", l balle =",parametre.x_tot_pendule_balle  )

print (" v raqu = ", v_raq/15  , "v bal = ", v_bal)

#du coté ucontroleur l'adc :  ref =3.3v pour 8191 donc chaque niveau = 0.000403 V
# pour le capteur 5V pour 360 degres , donc pour chaque degres on a = 0.0139 V
# on a donc pour la raquette 
# sur 5 ms un deplacmemnt de 69 donc 69 * 0.000403 = 0.0278V
#soit 0.0278 /0.0139=2degres , donc 2 degres pour 5ms soit une periode de 0.9
# saxchant que v  = 2PIr / T(periode)
# pour la balle
# sur 5 ms on a un deplacement de 34 soit 0,0137mV

#donc v = 2*3.141*0.6/0.9
#v_raq_radseconde : 2 degres = 0,0349066 rad
# 2 degres pour 5 ms font 400degres /s  
v_raq_kmh =(2*math.pi*0.6)/0.9
print ("vitesse km/h  = ", v_raq_kmh*3.6) 
E_raq = E_cin(parametre.m_tot_tennis, parametre.x_tot_tennis ,12 )
E_bal = E_cin(parametre.m_tot_pendule_balle , parametre.x_tot_pendule_balle , v_bal/5)
#E_raq = E_cin(parametre.m_tot_tennis, 1 , v_raq/5)
#E_bal = E_cin(parametre.m_tot_pendule_balle , 1 , v_bal/5)
print ("E raq = ", E_raq , ", E balle = ",E_bal , ", ratio =" , E_bal / E_raq)

#root.destroy()
print (len(data_int_1))
plt.figure(1)
plt.subplot(411)
plt.ylabel('teta')

plt.plot(temps, data_bal, 'bo', temps, data_bal, 'k')
plt.subplot(412)
plt.xlabel('temps')
plt.ylabel('degres/s')
plt.plot(temps, data_raq, 'bo', temps, data_raq, 'k')

plt.subplot(413)
plt.xlabel('temps')
plt.ylabel('degres/s')
plt.plot(size_int, data_int_1, 'bo', size_int, data_int_1, 'k')

plt.subplot(414)
plt.xlabel('temps')
plt.ylabel('degres/s')
plt.plot(size_int, data_int_2, 'bo', size_int, data_int_2, 'k')

plt.show()