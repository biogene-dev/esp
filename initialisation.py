# -*- coding: utf-8 -*-
"""
Created on Wed May  8 19:36:26 2019

@author: ash
"""
# permet de redémarrer le programme et supprimer les variables de la mémoire
"""
clear all
close all
fclose all
clc
"""
global nomPortCOM, nomRaquette, g, optionsOde

# ------------------------------------------------------------
# nom du port COM sur lequel est connecté la carte, ex : 'COM3'
nomPortCOM='COM3';
# ------------------------------------------------------------

# variables globales :
g=9.81; # gravité

# options pour les interpolations
"""
optionsOde=odeset('AbsTol',1e-8,'RelTol',1e-5);
"""