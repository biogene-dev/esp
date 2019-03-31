# -*- coding: utf-8 -*-
"""
Created on Mon Mar 18 18:13:05 2019

@author: ash
"""
import math
# sphinx_gallery_thumbnail_number = 3
import matplotlib.pyplot as plt
import numpy as np
import parametre

phi_zero = math.pi/2#angle de lancement
g = 9.81
l = parametre.x_tot_tennis # longueur de l'axe en metres

x=[]
y1=[]
v1 =[]
E=[]
phi_t1 = phi_zero
puls0 = 2*math.pi*math.sqrt(l/g)

E_impact = 0
m = parametre.m_tot_tennis
print ("moment inerti" , parametre.x_tot_tennis)
E_pot = m * g * (l) * (1-math.cos(phi_zero))

def E_cin (V_phi):
    return  (0.5 * m * math.pow(l,2) * math.pow(V_phi,2))


#for i in range (phi_zero,0,-1):
for i in range (0,int(math.pi*100.0),1):
    t=i/100
    print (t)
    phi= phi_zero * (math.cos (puls0 * t))
    #print (phi)
    
    if phi < 0 : break
    x.append(t)
    y1.append(phi)
    print ("phi = ",phi, " ; phi-1 =", phi_t1)
    delta_phi = abs(phi - phi_t1)*800
    v1.append((delta_phi))
    e_cin = E_cin(delta_phi)
    E.append(e_cin)
    #print ("v1 = ",v1)
    phi_t1 = phi
    #y2.append(phi)

print ("E potentiel a t=0 : ", E_pot, "Joules")
print ("E cinetique a tmax : ", e_cin, "Joules")

plt.figure(1)
plt.subplot(211)
plt.ylabel('teta')

plt.plot(x, y1, 'bo', x, y1, 'k')

plt.subplot(212)
plt.xlabel('temps')
plt.ylabel('degres/s')
plt.plot(x, E, 'r--')

plt.show()


'''
plt.plot(x, y1, label='wo')
#plt.plot(x, y2, label='w1')

plt.xlabel('x label')
plt.ylabel('y label')

plt.title("Simple Plot")

plt.legend()

plt.show()

vit.plot ( x , v1, label="vit")
vit.show()
'''