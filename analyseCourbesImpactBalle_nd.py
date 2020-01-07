# ------------------------------------------------------------
# ------------------------------------------------------------
#                      ESP THORBOT
#                     décembre 2012
# ------------------------------------------------------------
# ------------------------------------------------------------
# auteur : benjamin Chouvion - CAPSULE
# contact : bchouvion@capsule-ea.fr
#           04 84 25 07 00
# ------------------------------------------------------------
# ------------------------------------------------------------

# analyse de la courbe de réponse coté balle après impact de la raquette


## ========================================================================
# Recupère les données enregistrées lors du calibrage
# =========================================================================
"""
if ~exist('rendement_balle.mat') # un calibrage doit avoir été fait
    error('faire calibrage balle d''abord');
end
if ~exist('mLI_balle.mat')
    error('faire calibrage balle d''abord');
end
load('rendement_balle.mat');
load('mLI_balle');
"""
import os.path
import numpy as np

import math
# sphinx_gallery_thumbnail_number = 3
import matplotlib.pyplot as plt

if os.path.isfile('mLI_balle.py'):# recupere précédentes valeurs
    from mLI_balle import mLI_b 
else :
    print ('faire calibrage balle d''abord')
    exit()
if os.path.isfile('rendement_balle.py'):# recupere précédentes valeurs
    from rendement_balle import rendement_b 
else :
    print ('faire calibrage balle dabord')
    exit()
    
import scipy.io
# open file
from pathlib import Path
import numpy as np
a = np.array([])
t = np.array([])
data_folder = Path("source_data/text_files/")
p = Path('.')
print (p.absolute())
file_to_open = p.absolute() / "impact_nd.txt"

## ========================================================================
# cherche le nombre de tests et les intervales d'intérêt
# =========================================================================

#t = mesureBalle[][0]*1e-6 # en secondes

print ("opening file")
print(file_to_open)
f= open(file_to_open,"r")
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
            val = float(to_test[3])
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
            if len(to_test) == 6 :
                transAng=(val*360/(8192*0.8)); # le 0.8 est du à la bande de 10# -> 90# du capteur
#            transAng=(val*360/((2^13)*0.8)); # le 0.8 est du à la bande de 10# -> 90# du capteur
                if transAng < 10000 : a =np.append(a, transAng)
                if int(to_test[0])< 9999999999 : t= np.append(t , int (int(to_test[0])/1000) ) 

#            print (len(a))
        
print ("len de rec =" , len(a))
print (a)

#ang = a; # angle en radian
ang=np.radians(-a); # angle en radian

# fait la moyenne des n premieres valeurs pour trouver celle au repos
#ang0=mean(ang(t<(tpsRepos),1));
ang0 = (np.mean(ang[0:1000]))

print ("angle 0 = ", ang0)
print ("ang0:10n=" , ang[0:10])
# soustrait à toutes les autres pour mettre le repos sur angle=0
#ang=ang-ang0;
ang = np.subtract(ang, ang0)
print ("ang0:10n=" , ang[0:10])

# cherche les pics de la courbe d'oscillations libres amorties
#ind=peakfinder(ang,0.3);
#ind=[1;ind]; # rajoute l'indice 0 en début


#
#t=mesureBalle(:,1)*1e-6; # en secondes
#ang=deg2rad(-mesureBalle(:,2));  # en rad et opposé (signe '-') pour avoir avec v0>0
#

## fait la moyenne des premieres valeurs pour trouver la valeur au repos
#ang0=mean(ang(t<(tpsRepos),1));
## soustrait à toutes les autres pour mettre le repos sur angle=0
#ang=ang-ang0;
#
## enleve les outliers (fonction effectuées plusieurs fois)
#for ii=1:6
#    ang=RemoveOutlier(ang,deg2rad(1));
#end
#
#
angMin = np.radians(30); # angle en radian

angMax = np.radians(150); # angle en radian
print("max admis", angMax)
print ("max = ",max(ang))
print("min admis", angMin)
print ("min = ",min(ang))

a = 0
if max(ang)< angMax:
    print('changer les valeurs de angMin et angMax car je vais calculer l''energie sur des points après le premier pic (non recommandé)')

## ------ pour chercher les pics d'impact
## cherche les endroits ou les angles sont sup à angMin deg et inf à angMax
#angMin=deg2rad(3); # deg
#angMax=deg2rad(16); # deg
## 
#if max(ang)<max(angMax)
#        warning('changer les valeurs de angMin et angMax car je vais calculer l''energie sur des points après le premier pic (non recommandé)')
#end
## 

ind=[]

compt = 1 
ind.append(0)

for i  in ang : 
    if (i>angMin and i<angMax) : 
        ind.append(compt)
#        print (ind)
    compt+=1
print ("ind " , t[ind])    
print ("lzen ind " , len(ind), "sur ", len(ang))    
#ind=[1;find(and(ang>angMin,ang<angMax))];
#
#indSub=ind(2:end)-ind(1:end-1); # temps entre 2 mesures
indic = np.asarray(ind)
    
indSub = indic[1:] - indic[:-1]
print ("lzen indsub " , len(indSub))    

for i in indSub:
    if i>1 :
        print ('igrand',i)
print ("ind " , indSub)    

## correspond à un debut d'impact (nouveau test) si le nombre entre 2 mesures est grand
#indTestDebutB=ind(find(indSub>(tpsEntreMesures/(t(2)-t(1))))+1); # 
#mat["tpsEntreMesures"]

tpsEntreMesures = 3 
indTestDebutB = []
compt  = 0
for i in range (len(indic[:-1])) :
    if indSub[compt] > (tpsEntreMesures/(t[2]-t[1])) :
        indTestDebutB.append(indic[i+1])
        print ("le gaggnat",i)
    compt+=1
   
nbrTests = len (indTestDebutB)
print ("nbrde test = " ,len(indTestDebutB) )
print ("temps debuts impact",indTestDebutB)
#nbTests=length(indTestDebutB);
#
## garde n points apres chaque debut tels que la valeur finale soit inferieure à angMax
#for ii=1:nbTests
#    indTestFin(ii,1)=indTestDebutB(ii);
#    compteurLoc=0;
#    while (compteurLoc<1000) && (ang(indTestFin(ii))<angMax) # max 1000 points
#        indTestFin(ii)=indTestFin(ii)+1;
#        compteurLoc=compteurLoc+1; 
#    end
#end
#
vitesseBallekmh=[]
indTestFin= np.arange(nbrTests)
for ii in range (nbrTests) :
#    print("ii=,",ii)
#    print("indTestFin=,",indTestFin[ii])
#    print("indTestDebutB=,", indTestDebutB[ii])

    indTestFin[ii] = indTestDebutB[ii]
    compteurLoc=0
    while (compteurLoc<100) and (ang[indTestFin[ii]]<angMax): # max 1000 points
        indTestFin[ii]+=1;
        compteurLoc=compteurLoc+1; 

## # trace la figure avec toutes les mesures et les zones d'impact étudiées
#figure
#plot(t,rad2deg(ang),'.b');
#hold on 
#xlabel('temps (s)')
#ylabel('angle (deg)') 
#    
#for ii=1:nbTests
#    # plot les sélections sur la figure principale
#    plot(t(indTestDebutB(ii):indTestFin(ii)),rad2deg(ang(indTestDebutB(ii):indTestFin(ii))),'.r');
#end
#
#
###
##========================================================================
## étudie chaque test indépendamment
##=========================================================================
#indAEnlever=[];
#Energie_balle=zeros(1,nbTests); # initialisation
#Energie_balle_cplt=zeros(size(mLI_balle,1),nbTests); # initialisation
#
#for n=1:nbTests
#    t_test=t(indTestDebutB(n):indTestFin(n));
#    ang_test=ang(indTestDebutB(n):indTestFin(n)); 
#    # recale de telle manière que le premier temps soit=0
#    t_test=t_test-t_test(1);
#    
#    #----- trace
#    figure # si je veux que ca soit sur des nouveaux graphes
##     ha=gca; # recupère le current axis handle
#    plot(t_test,rad2deg(ang_test),'.r')
#    xlabel('Temps (s)')
#    ylabel('Angle (deg)')    
#    title(['Test numero: ' int2str(n)])
#    hf(n)=gcf; # recupère le current figure handle
#
#    
#    #----- check si l'utilisateur est OK avec ce test
#    testValid=1;        # ajouté par Pierrick
#    # testValid=input('Doit obtenir une ligne croissante. Test OK ? (oui ==> default / non ==> 0) ===> ');
#    if testValid==0
#        disp('on oublie ce test et passe à celui d''après'); 
#        indAEnlever=[indAEnlever,n];
#        continue 
#    end
Energie_balle=np.zeros(nbrTests) # initialisation
Energie_balle_cplt=np.zeros([nbrTests]) # initialisation
for n in range (nbrTests):
    print("n=,",n)
    print("indTestFin=,",indTestFin[n])
    print("indTestDebutB=,", indTestDebutB[n])
    t_test=t[int(indTestDebutB[n]):int(indTestFin[n])];
    ang_test=ang[int(indTestDebutB[n]):int(indTestFin[n])]; 
#    # recale de telle manière que le premier temps soit=0
    t_test=t_test-t_test[1];
    
#    #----- interpolation lineaire de la mesure sur ce petit intervalle
#    p=polyfit(t_test,ang_test,1)
    p=np.polyfit(t_test,ang_test,1)
    mesRegress = p[0]*t_test+p[1]
#    mesRegress=p(1)*t_test+p(2);
#    hold on
#    plot(t_test,rad2deg(mesRegress),'-c'); # trace l'interpolation linéaire
#    
#    # la vitesse est egale au coefficient directeur de la droite (provenant de l'interpolation linéaire)
#    omega=p(1); # vitesse en rad/s
#    vitesse=omega*0.78;    # en lineaire (m/s) pour information
#    vitesseBallekmh(n)=vitesse*3.6;
    omega =p[0]
    vitesse = omega*0.78
    vitesseBallekmh.append(vitesse*3.6)
    theta = np.mean(ang_test)
#    theta=mean([ang_test(1),ang_test(end)]); # en radian
#    
#    disp(['Test numero ' int2str(n) ',      vitesse balle à ' num2str(round2(rad2deg(theta),0)) ' deg après impact : ' num2str(round2(vitesse*3.6,1)) ' km/h']);
#    
#    #----- énergie après impact    
#    for ii=1:size(mLI_balle,1) # pour chaque valeurs de I qui viennent des tests de calibrage
#        mL=mLI_balle(ii,1);
#        I=mLI_balle(ii,2);
#        Energie_balle_cplt(ii,n)=1/2*I*omega^2+mL*g*(1-cos(theta)); # E=1/2*I om^2+mgh (cinétique + potentielle)
#        # énergie sans considérer le rendement => énergie réelle
#    end
#    
#end
#
## close(hf);   # ferme toutes les figures (les droites)
#
## moyenne des énergies sur les x valeurs de mLI
#Energie_balle=mean(Energie_balle_cplt,1);
#disp('*******************************************************')
#disp('Energie (J) transmise dans la balle (sans perte) est :')
## disp('j''eme colonne pour test d''impact numero j')
#Energie_balle
    print ("Test numero =", n , " / vitesse balle a = " , theta ," /degres apres impact=", vitesse*3.6 , "km/h = " )
#    #----- énergie après impact    
    for ii in range (len(mLI_b)):
        g = 9.81
        mL=mLI_b[ii][0]
        I=mLI_b[ii][1]
        print ("ml = ", mL , " I = ", I)
        Energie_balle_cplt[n] = (1/2*I*math.pow(omega,2)+mL*g*(1-math.cos(theta)))
#Energie_balle
#
Energie_balle=np.mean(Energie_balle_cplt);
print('*******************************************************')
print('Energie (J) transmise dans la balle (sans perte) est :')
print('eme colonne pour test d''impact numero j')
print(Energie_balle)        
plt.figure(1)
plt.figure( figsize=(18, 16))
#plt.gcf().subplots_adjust(wspace = 0, hspace = 4)
#    plot(t(indTestDebutB(ii):indTestFin(ii)),rad2deg(ang(indTestDebutB(ii):indTestFin(ii))),'.r');

plt.subplot(311)
plt.ylabel('teta')
plt.plot(t , ang,'k' )
for ii in range(nbrTests):
    plt.plot(t[int(indTestDebutB[ii]):indTestFin[ii]],ang[int(indTestDebutB[ii]):indTestFin[ii]],'.r');
plt.plot(t[indTestDebutB],ang[indTestDebutB],'.g')
    
plt.subplot(312)
plt.plot(t_test,mesRegress,'-c');
plt.subplot(313)
plt.plot(t[indic] , ang[indic],'r' )

plt.show()