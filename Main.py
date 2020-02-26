# -*- coding: utf-8 -*-
"""
Created on Wed May  8 19:19:26 2019

@author: ash
"""
from pathlib import Path


# programme principal à lancer en premier
# ------------------------------------------------------------
#                   INITIALISATION
# ------------------------------------------------------------
import initialisation 
global temps
global raquette
global balle
global accelero

# ====================================
#           variables
# ====================================

# secondes (temps nécessaire entre 2 mesures, avec système au repos)
tpsEntreMesures = 3; 

# temps (en seconde) nécessaire avant le début d'une série d'acquisition
# (système au repos)
tpsRepos = 1; 

# nb d'oscillations minimum à imposer lors du calibrage
nbOscillations = 7; 

##Rajouté par Pierrick
#PositionBrute = input('Distance (en m) entre l''axe du codeur angulaire 
#raquette et l''axe du haut de la clé allen qui fixe le début du 
#manche ? : ');
##exemple de valeur = 0.261

# ----------------------------
## variables géométrique de conception


## Pendule Balle

# Masse tige pendule (kg)
mTigePendule = 1.212; # PB et GR le 17/11/17

# mTigePendule = 1.198;  # Ancien data (Capsule)


# Distance du CG de la tige sur le pendule Balle avec l'axe de rotation
# (mètre)
posCdGTigePendule = 0.385; # PB et GR le 17/11/17

# posCdGTigePendule = 0.391; # Ancien data (Capsule)


# Masse de la boule qui représente la balle (kg)
mBoule = 0.290; # changé par PB et GR le 17/11/17

# mBoule = 0.346; # Ancien data (Capsule)


# Distance du CG de la boule sur le pendule Balle avec l'axe de rotation
# (mètre)
posCdGBoule = 0.777; # PB et GR le 17/11/17

# posCdGBoule = 0.78; # Ancien data (Capsule)




## Pendule Raquette

# Poids du pendule Raquette (avec 26g de rondelles au niveau des vis qui 
# tiennent le manche, rajoutées par PB le 26/02/2014)
#mBras = 4.892; # PB et GR le 06/02/18

# mBras = 3.95; # Gauvain Touch padel

# Autres possibilités
# mBras = 3.790; # poids initial du pendule tennis sans l'aimant
mBras = 4.163; # poids du pendule tennis avec l'aimant installé dessus

# mBras = 4.563; # Poids du pendule padel 





# Distance du CG du pendule raquette avec l'axe de rotation

# avec le module tennis
posCdGBras = 0.132; # PB et GR le 06/02/18

# avec le module Padel
# posCdGBras = 0.260; # PB et GR le 06/02/18


# Autres possibilités
# posCdGBras = 0.121; # avec module tennis et sans l'aimant fixé
# posCdGBras = 0.247; # avec le module padel





## Position de la raquette par rapport à l'axe de rotation

# Distance entre l'axe de rotation du codeur angulaire et le haut de 
# la clé allen qui fixe la butée du manche (mètre). 
posRaquette = 0.269; # raquettes de longueur 685mm
posRaquette = 0.26; # raquettes de longueur 685mm


# Autres possibilités
# posRaquette = 0.286; # raquettes de longueur 660mm
# 
# posRaquette = 0.254; # raquettes de longueur 700mm
# posRaquette = 0.210; #Position Haute (partie haute du tamis)
# posRaquette = 0.250; #Position Standart (partie milieu du tamis)
# posRaquette = 0.290; #Position Basse (partie basse du tamis)
# posRaquette = 0.440; # raquettes de padel de longueur 460mm



# E2 Lab 2015

# Position spéciale pour les raquettes de taille différentes
# posRaquette=0.207; #Position Haute Pure drive +
# posRaquette=0.237; #Position Standart Pure drive +
# posRaquette=0.267; #Position Basse Pure drive +
# posRaquette=0.240; #Position Haute Pure drive Courte
# posRaquette=0.270; #Position Standart Pure drive Courte
# posRaquette=0.300; #Position Basse Pure drive Courte
# posRaquette=0.280; #Position Basse Pure Drive (que 3cm de débattement)
# posRaquette=0.220; #Position Haute Pure Drive (que 3cm de débattement)



# Babolat Schelle 2017

# # Pour la raquette PD Lite (longueur 685 mm)
# posRaquette=0.269; 

# # Pour la raquette PD 110, Hammer head et Isometric (longueur 660 mm)
# posRaquette=0.286; 




## On détermine les chemins et/ou current directory où on écrit les données
data_folder = Path("source_data/text_files/")

file_to_open = data_folder / "raw_data.txt"

# Celui où on stocke les données brutes
data_folder = '\Validation_Energie_Balle\Donnees_brutes';

# Celui où sont stockés les scripts 
script_folder = '\Validation_Energie_Balle\esp_thorbot_install\python';






## ---------------------------
global typeTest 


print('*******************************************************')
print('Choisir le type de test à effectuer')
print('1     calibrage côté balle');
print('2     calibrage côté raquette');
print('3     mesure de l''impact');
typeTest=input ('type de test =======>  ');



##
#switch typeTest ;# type de test choisi par l'utilisateur
        
## -----------------------------------------------------------
#               TEST CALIBRAGE BALLE
# ------------------------------------------------------------
# lancer la balle seule en pendule libre, sûr de petites oscillations 
# angle initial ~ 30-50 degrés
print (" typeTest ===",typeTest)
typeTest = int (typeTest)
if typeTest == 1 :
    print ("===TEST CALIBRAGE BALLE===")
    exec(open("acqui.py").read())
    exec(open("analyseCourbesCali_nd.py").read())

#    mesures = Acquisition(100,typeTest,current);   # 100 secondes de test par default pour calibrage balle
    # fait l'acquisition et save data dans calBalle-date-.csv
    
#    analyseCourbesCali   # fait l'analyse de cette courbe

    # en output il sort:
#         mLI_balle.mat, rendement_balle.mat



## -----------------------------------------------------------
#               TEST CALIBRAGE RAQUETTE
# ------------------------------------------------------------
if typeTest == 2 :  # raquette
        
        # lancer la raquette en pendule libre, sur de petites oscillations 
        # angle initial ~ 30-50 degrés
    print ("===TEST CALIBRAGE RAQUETTE===")
    exec(open("acqui.py").read())
    exec(open("analyseCourbesCali_nd.py").read())
    
#        nomRaquette=input('nom de la raquette ? ==> ','s'); # demande nom de la raquette
#        save('nomRaquette','nomRaquette'); # enregistre en dur sous fichier .mat
#        mesures = Acquisition(80,typeTest,current); # 80sec de temps d'acquisition
        # fait l'acquisition et save data dans calibBalle-nomRaquette-date.csv

#        analyseCourbesCali   # fait l'analyse de cette courbe
        
        # en output il sort:
#         mLI_raquette.mat, rendement_raquette.mat, donneesGeomRaquette.mat
        
        

## -----------------------------------------------------------
#               ANALYSE COURBE REPONSE IMPACT
# ------------------------------------------------------------
if typeTest == 3 :
    print ("===IMPACT===")

# va utiliser les données renseignées par la calibration côté balle et côté
# raquette
# type test 1 et type test 2 donc doivent être effectués avant
    import os
#        mesures = Acquisition(150,typeTest,current); # output :1x3 cell avec les mesures angles + l'accelero
    exec(open("acqui.py").read())
    exec(open("analyseCourbesImpactBalle_nd.py").read())  # output: Energie_balle (pour n tests)
                  
    exec(open("analyseCourbesImpactRaquette_nd.py").read()) # outpu: Energie_raquette (pour n tests)
    exec(open("analyseCourbesAccelero_nd.py").read()) # outpu: Energie_raquette (pour n tests)
    print ("Energie_balle_cplt = ", Energie_balle_cplt)
#    testAcc= 1;     # Ajouté par Pierrick
#    # testAcc=input('faire l''analyse sur l''accéléro ? 1=oui (0=non) ===> '); # pour pas faire automatiquement
#    # l'analyse sur l'accéléromètre (ex s'il est pas branché...)
#    if testAcc :
#        analyseCourbesAccelero
#    end
        
## -----------------------------------------------------------
#               CALCUL COEFF RESTITUTION
# ------------------------------------------------------------
# utilise les données (Energies) sorties par les routines analyseCourbesImpactBalle
# et analyseCourbesImpactRaquette
#
## enlever les resultats des tests non concluants (ceux non validés par
## l'utilisateur)
#        Energie_balle(sort(indAEnlever))=[];
#        Energie_raquette(sort(indAEnlever))=[];
#        Energie_balle_cplt(:,sort(indAEnlever))=[];
#        Energie_raquette_cplt(:,sort(indAEnlever))=[];
#
#
#        print('*******************************************************')
#        print('Les coefficients de restitution (#) trouvés sont :');
#        restitution=Energie_balle./Energie_raquette*100   # rapport des 2 énergies
#        #==> ce sont les coefficients de restitution énergétique (et non de
#        #vitesse)
#
#        # possibilité de sauvegarder en dur les résultats
#        testSauve=1;        # ajoutée par Pierrick
##         testSauve=input('Voulez-vous sauvegarder les coefficients de restitution ? (oui => 1 / non par default) ');
#        switch testSauve
#            case 1
#                fid=fopen('restitutionResultats.csv','a');
#                maDate=datestr(now);
#                fprintf(fid,'#s\n',maDate);
#                load('nomRaquette.mat');
#                fprintf(fid,'#s\n',nomRaquette);
#                str='';
#                for ii=1:length(restitution)
#                    str=[str '#4.2f '];
#                end
#                str=[str '\n'];
#                fprintf(fid,str,restitution');
#                fclose(fid);
#
#            otherwise # par default, on n'enregistre rien
#        end
#
### -----------------------------------------------------------
##               Autre choix pas possible
#    otherwise
#        print('Mauvaise valeur du type de test, choisir entre 1 et 3')
#end
#
####Pierrick a crée les toutes les lignes suivantes. 
#switch typeTest # type de test choisi par l'utilisateur
#case 3
##On crée une version pour copier-coller sur Excel
#zz1_Restitution = restitution';
#zz2_AmpliPP = amplitudePP';
#zz3_EnergieRaquette = Energie_raquette';
#zz4_VitesseRaquette = (vitesse*3.6)';
#zz5_EnergieBalle = Energie_balle';
#zz6_VitesseBalle = vitesseBallekmh';
#zz7_NRJ_Totale = valEnergieFreqTout';
#zz8_NRJ_60a80Hz = valEnergieFreqBande';
#zzz(:,1) = zz1_Restitution;
#zzz(:,2) = zz2_AmpliPP;
#zzz(:,3) = zz3_EnergieRaquette;
#zzz(:,4) = zz4_VitesseRaquette;
#zzz(:,5) = zz5_EnergieBalle;
#zzz(:,6) = zz6_VitesseBalle;
#zzz(:,7) = zz7_NRJ_Totale;
#zzz(:,8) = zz8_NRJ_60a80Hz;
#
##On enregistre la courbe de l'angle de la raquette
#testSauve=1;        # ajoutée par Pierrick
## testSauve=input('Voulez-vous sauvegarder toutes les datas de cette raquette ? (oui => 1 / non par default) ');
#cd(current.data) # Je change mon current directory
#save('nomRaquette.mat','nomRaquette') #On a besoin du nom de la raquette
#        switch testSauve
#            case 1
#                load('nomRaquette.mat');
#                save(nomRaquette)
#            otherwise # par default, on n'enregistre rien
#        end
#end
#cd(current.script) # je restaure mon ancien current directory