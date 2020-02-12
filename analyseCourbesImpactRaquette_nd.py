# -*- coding: utf-8 -*-
# coding: utf-8

print ("start")

# ------------------------------------------------------------
# ------------------------------------------------------------
#                      ESP THORBOT
#                     decembre 2012
# ------------------------------------------------------------
# ------------------------------------------------------------
# auteur : benjamin Chouvion - CAPSULE
# contact : bchouvioncapsule-ea.fr
#           04 84 25 07 00
# ------------------------------------------------------------
# ------------------------------------------------------------

# analyse de la courbe de reponse cote raquette apres impact


## ========================================================================
# Recupere les donnees enregistrees lors du calibrage
# =========================================================================
print ("start ENERGIE RAQUETTE")
"""
if exist('rendement_raquette.mat','file')
    disp('*******************************************************')
    error('faire calibrage raquette d''abord');
end
if exist('mLI_raquette.mat','file')
    error('faire calibrage raquette d''abord');
end
load('mLI_raquette');
load('rendement_raquette');
load('donneesGeomRaquette');
""" 
import os.path
import numpy as np

import math
# sphinx_gallery_thumbnail_number = 3
import matplotlib.pyplot as plt

if os.path.isfile('mLI_Raquette.py'):# recupere precedentes valeurs
    from mLI_raquette import mLI_r 
else :
    print ('faire calibrage Raquette dabord')
#    exit()
if os.path.isfile('rendement_raquette.py'):# recupere precedentes valeurs
    from rendement_raquette import rendement_r 
else :
    print ('faire calibrage Raquette dabord')
#    exit()
    
import time

## ========================================================================
# cherche le nombre de tests et les intervales d'interet
# =========================================================================
import scipy.io
# open file
from pathlib import Path
import numpy as np
a = np.array([])
t = np.array([])
#data_folder = Path("source_data/text_files/")
#p = Path('.')
#print (p.absolute())
#file_to_open = p.absolute() / "impact_nd.txt"
#
#print ("opening file")
#print(file_to_open)
#f= open(file_to_open,"r")
#(f.read(1000)) 
#val = 0.0 
#print ("Formatage des datas")
#rating = 0
#for line in f:
#    rating +=1
#    if line :
#        try :
##            print (line)
#            to_test =line.split(',')
#            val = float(to_test[2])*-1
##            print ("val =",val)
#        except :
#            pass
#            print("bad")
##            print (line)
#        else :
##            print ("good" , val)
#            # ----------------------------
#            # conversion binaire, valeur reelle
#    
#            # fonction de transformation bit -> degre
#            if len(to_test) == 6 :
#                transAng=(val*360/(8192*0.8)); # le 0.8 est du à la bande de 10# -> 90# du capteur
##            transAng=(val*360/((2^13)*0.8)); # le 0.8 est du à la bande de 10# -> 90# du capteur
#                if transAng < 10000 : a =np.append(a, transAng)
#                if int(to_test[0])< 9999999999 : t= np.append(t ,  (float(to_test[0])/1000000) ) 
#
##            print (len(a))

global temps
global raquette
global balle
global accelero

a = raquette
t = temps
print ("len de rec =" , len(a))
print (a)
t = t -t[0]
#ang = a; # angle en radian
ang=np.radians(-a); # angle en radian

# fait la moyenne des n premieres valeurs pour trouver celle au repos
#ang0=mean(ang(t<(tpsRepos),1));
ang0 = (np.mean(ang[0:1000]))
print ("angle 0 = ", ang0)
print ("ang0:10n=" , ang[0:10])
# soustrait a toutes les autres pour mettre le repos sur angle=0
#ang=ang-ang0;
ang = np.subtract(ang, ang0)
print ("ang0:10n=" , ang[0:10])
time.sleep(3)

angMin = np.radians(-30); # angle en radian

angMax = np.radians(-20); # angle en radian
print("max admis", angMax)
print ("max = ",max(ang))
print("min admis", angMin)
print ("min = ",min(ang))

a = 0
if max(ang)< angMax:
    print('changer les valeurs de angMin et angMax car je vais calculer lenergie sur des points apres le premier pic (non recommande)')

# ------ pour chercher les pics d'impact
# cherche les endroits ou les angles sont sup a angMin deg et inf a angMax
#angMin=deg2rad(-30); # rad         
#angMax=deg2rad(-20); # rad
#Pierrick a change Ang Max pour être bien avant l'impact (en se mettant de -30 a -20)

#[indPic,valMaxPicsRaq]=peakfinder(ang,0.2,deg2rad(0.5)); # isole les pics
###Pierrick : changer l'angle de 0.5 si besoin. C'est l'angle de detection
###que fait la raquette avec la verticale apres l'impact. 
    ###02/04/2015 : Pierrick vient de changer cet angle : anciennement 1, il
    ###vient d'ecrire 0.5 car la raquette en cours n'atteint pas 1.


# selectionne que les pics reels (et non les rebonds)
#ind=[1;indPic]; # rajoute l'indice 0 en debut

# cherche ceux pour lesquels il n'y a pas de pics juste avant
#indSub=ind(2:end)-ind(1:end-1);
    
from detect_peaks import detect_peaks
indexes = detect_peaks(ang, mph=0.3 , mpd=1000)

ind=[]
from detect_peaks import detect_peaks
indPic = detect_peaks(ang, mph=0.02 , mpd=1000)
compt = 1
ind = indPic 
np.insert(ind , 0 ,0 )

#% cherche ceux pour lesquels il n'y a pas de pics juste avant
#indSub=ind(2:end)-ind(1:end-1);
    
indSub = ind[1:] - ind[:-1]
print ("lzen indsub " , len(indSub))    
    
#for i in indSub:
#    if i>1 :
#        print ('igrand',i)
print ("indSub " , indSub) 


# pas de pics si le nb de points entre 2 max consecutifs est grand
#periode=indSub*(t(2,1)-t(1,1)); # diff des indices *  tps increment => periode en us
#indPic=ind(find((periode)>tpsEntreMesures)+1);
#
#nbTestRaq=length(indPic);
#if nbTests~=nbTestRaq; error('pas le meme nb dimpact trouves cote raquette et cote Raquette'); end
periode = indSub*(t[2]-t[1])
tpsEntreMesures = 3 
indPic = ind[np.where((periode)>tpsEntreMesures)]

nbTestRaq=len(indPic);
print ("le nombre de test raq =",nbTestRaq)

plt.figure(1)
plt.figure( figsize=(8, 6))
#plt.gcf().subplots_adjust(wspace = 0, hspace = 4)
#    plot(t(indTestDebutB(ii):indTestFin(ii)),rad2deg(ang(indTestDebutB(ii):indTestFin(ii))),'.r');

plt.ylabel('teta')
plt.plot(t , ang,'k' )
plt.plot(t[indPic] , ang[indPic],'ro')
plt.show()

indTestDebut = []
indTestFin = []
compt  = 0
for i in range (nbTestRaq) :
    print ("pour le test n= ",i)
    indTestFin.append(indPic[i])
    while ang[indTestFin[i]] > angMax : 
        indTestFin[i] = indTestFin[i] - 1 
    indTestDebut.append(indTestFin[i])
    while ang[indTestDebut[i]] > angMin :
        indTestDebut[i] = indTestDebut[i] -1 
        
print ("nbrde test = " ,len(indTestDebut) )
print ("temps debuts impact",indTestDebut)

plt.figure(1)
plt.figure( figsize=(8, 6))
plt.plot(t , ang,'k' )

for ii in range(nbTestRaq):
    plt.plot(t[int(indTestDebut[ii]):indTestFin[ii]],ang[int(indTestDebut[ii]):indTestFin[ii]],'.r');
plt.show()

# part de ces pics et recule dans les increments pour recuperer les indTestDebut et indTestFin juste avant

#for ii=1:nbTestRaq
#    indTestFin(ii,1)=indPic(ii); # initialisation au pic
#    compteurLoc=0;
#    while ang(indTestFin(ii))>angMax
#        indTestFin(ii)=indTestFin(ii)-1;
#    end
#    indTestDebut(ii,1)=indTestFin(ii); # initialisation
#    while ang(indTestDebut(ii))>angMin
#        indTestDebut(ii)=indTestDebut(ii)-1;
#    end
#end



## ========================================================================
#  etudie chaque test independamment
# =========================================================================
#Energie_raquette=zeros(1,nbTests); # initialisation
#Energie_raquette_cplt=zeros(size(mLI_raquette,1),nbTests); # initialisation
Energie_raquette=np.zeros(nbTestRaq) # initialisation
Energie_raquette_cplt=np.zeros([nbTestRaq]) # initialisation
#
#for n=1:nbTests
#    t_test=t(indTestDebut(n):indTestFin(n));
#    if isempty(t_test)
#        continue  # je saute ce cas si il est empty
#    end
#    
#    ang_test=ang(indTestDebut(n):indTestFin(n)); 
#    # recale de telle maniere que le premier temps soit=0
#    t_test=t_test-t_test(1);
#    
#    #----- trace
#    figure
#    plot(t_test,rad2deg(ang_test),'.')
#    xlabel('temps (s)')
#    ylabel('angle (deg)')    
#    title(['test numero: ' int2str(n)])
#    
#    #----- check si l'utilisateur est OK avec ce test
#    testValid= 1;   #ajoute par Pierrick
#    # testValid=input('doit obtenir une ligne croissante. Test OK ? (oui ==> default / non ==> 0) ===> ');
#    if testValid==0
#        disp('on oublie ce test et passe a celui d''apres'); 
#        indAEnlever=[indAEnlever,n];
#        continue 
#    end
#    
#    #----- interpolation lineaire de la mesure sur ce petit intervalle
#    p=polyfit(t_test,ang_test,1);
#    mesRegress=p(1)*t_test+p(2);
#    hold on
#    plot(t_test,rad2deg(mesRegress),'-c');
#    


#    omega=p(1); # vitesse en rad/s
#    vitesse(n)=omega*0.78;    # en lineaire (m/s) pour information
#    
#    theta=mean([angMin,angMax]); # en radian
#
#    disp(['Test numero ' int2str(n) ',      vitesse raquette a ' num2str(round2(rad2deg(theta),0)) ' deg avant impact : ' num2str(round2(vitesse(n)*3.6,1)) ' km/h']);
#    
#    #----- calcul energie raquette
#    for ii=1:size(mLI_raquette,1) # pour chaque valeur de I qui viennt des tests de calibrage
#        mL=mLI_raquette(ii,1);
#        I=mLI_raquette(ii,2);
#        Energie_raquette_cplt(ii,n)=1/2*I*omega^2+mL*g*(1-cos(theta));
#    end
#    
#end

vitesseRaquettekmh=[]

for n in range (nbTestRaq):
    print("n=,",n)
    print("indTestFin=,",indTestFin[n])
    print("indTestDebutB=,", indTestDebut[n])
    t_test=t[int(indTestDebut[n]):int(indTestFin[n])];
    ang_test=ang[int(indTestDebut[n]):int(indTestFin[n])]; 
    # recale de telle maniere que le premier temps soit=0
    print ("t_test = ",t_test)
    t_test=t_test-t_test[0]
    
#    #----- interpolation lineaire de la mesure sur ce petit intervalle
#    p=polyfit(t_test,ang_test,1)
    p=np.polyfit(t_test,ang_test,1)
    mesRegress = p[0]*t_test+p[1]
    print('recherche accent')
#    mesRegress=p(1)*t_test+p(2);
#    hold on
#    plot(t_test,rad2deg(mesRegress),'-c'); # trace l'interpolation lineaire
#    
#    # la vitesse est egale au coefficient directeur de la droite (provenant de l'interpolation lineaire)
#    omega=p(1); # vitesse en rad/s
#    vitesse=omega*0.78;    # en lineaire (m/s) pour information
#    vitesseBallekmh(n)=vitesse*3.6;
    omega =p[0]
    vitesse = omega*0.78
    vitesseRaquettekmh.append(vitesse*3.6)
    theta = np.mean(ang_test)
#    theta=mean([ang_test(1),ang_test(end)]); # en radian
#    
#    disp(['Test numero ' int2str(n) ',      vitesse balle a ' num2str(round2(rad2deg(theta),0)) ' deg apres impact : ' num2str(round2(vitesse*3.6,1)) ' km/h']);
#    
#    #----- energie apres impact    
#    for ii=1:size(mLI_balle,1) # pour chaque valeurs de I qui viennent des tests de calibrage
#        mL=mLI_balle(ii,1);
#        I=mLI_balle(ii,2);
#        Energie_balle_cplt(ii,n)=1/2*I*omega^2+mL*g*(1-cos(theta)); # E=1/2*I om^2+mgh (cinetique + potentielle)
#        # energie sans considerer le rendement => energie reelle
#    end
#    
#end
#
## close(hf);   # ferme toutes les figures (les droites)
#
## moyenne des energies sur les x valeurs de mLI
#Energie_balle=mean(Energie_balle_cplt,1);
#disp('*******************************************************')
#disp('Energie (J) transmise dans la balle (sans perte) est :')
## disp('j''eme colonne pour test d''impact numero j')
#Energie_balle
    print ("Test numero =", n , " / vitesse balle a = " , theta ," /degres apres impact=", vitesse*3.6 , "km/h = " )
#    #----- energie apres impact    
    for ii in range (len(mLI_r)):
        g = 9.81
        mL=mLI_r[ii][0]
        I=mLI_r[ii][1]
        print ("ml = ", mL , " I = ", I)
        Energie_raquette_cplt[n] = (1/2*I*math.pow(omega,2)+mL*g*(1-math.cos(theta)))
#Energie_balle
#
Energie_raquette=np.mean(Energie_raquette_cplt);
print('*******************************************************')
print('Energie (J) transmise dans la raquette (sans perte) est :')
print('eme colonne pour test d''impact numero j')
print(Energie_raquette)        
#plt.figure(1)
#plt.figure( figsize=(8, 6))
##plt.gcf().subplots_adjust(wspace = 0, hspace = 4)
##    plot(t(indTestDebutB(ii):indTestFin(ii)),rad2deg(ang(indTestDebutB(ii):indTestFin(ii))),'.r');
#
#plt.subplot(311)
#plt.ylabel('teta')
#plt.plot(t , ang,'k' )
#for ii in range(nbTestRaq):
#    plt.plot(t[int(indTestDebutB[ii]):indTestFin[ii]],ang[int(indTestDebutB[ii]):indTestFin[ii]],'.r');
#plt.plot(t[indTestDebutB],ang[indTestDebutB],'.g')
#    
#plt.subplot(312)
#plt.plot(t[indic] , ang[indic],'r' )
#plt.subplot(313)
#plt.plot(t[indic] , ang[indic],'r' )
#
#plt.show()
