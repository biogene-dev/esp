# -*- coding: utf-8 -*-
"""
Created on Sat May 11 15:35:53 2019

@author: ash
"""
#import Main
tpsEntreMesures = 2; 
# nb d'oscillations minimum à imposer lors du calibrage
nbOscillations = 7; 

import math
# sphinx_gallery_thumbnail_number = 3
import matplotlib.pyplot as plt
import numpy as np
#import parametre

from tkinter import filedialog
from tkinter import *
import tkinter
#from Main import *

from matplotlib.backends.backend_tkagg import (
    FigureCanvasTkAgg, NavigationToolbar2Tk)
# Implement the default Matplotlib key bindings.
from matplotlib.backend_bases import key_press_handler
from matplotlib.figure import Figure

import numpy as np

# script utilisé lors de calibrage (raquette ou balle) pour trouver les propriétés
# de masse, inertie, et rendement des 2 axes
# les données sont dans la variable 'mesures' qui sort de la fonction 'Acquisition'
def plot_simple (to_plot):
    plt.figure(1)
    plt.figure( figsize=(8, 6))
    plt.plot(to_plot,'k' )
    plt.show()

# si calibrage raquette, demande autres informations

print ("=========calibration=======")
global typeTest
if typeTest==2 :
    print('*******************************************************')
    mRaquette=input('Masse de la raquette ? (en g) ==> ') 
    posCdGRaquette=input('Position du centre de gravité de la raquette ? (en mm) ==> ')
    print('*******************************************************')
#end
# open file
from pathlib import Path
import numpy as np
a = np.array([])
Temps=[]
data_folder = Path("source_data/text_files/")
p = Path('.')
print (p.absolute())
file_to_open = p.absolute() / "calib.txt"

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
            val = float(to_test[2])
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
                transAng=(val*360/(8192))*0.8; # le 0.8 est du à la bande de 10# -> 90# du capteur
#            transAng=(val*360/((2^13)*0.8)); # le 0.8 est du à la bande de 10# -> 90# du capteur
                if transAng < 10000 : a =np.append(a, transAng)
                if int(to_test[0])< 9999999999 : Temps.append(int(int(to_test[0])/1000))

#            print (len(a))
#global temps
#global raquette
#global balle
#global accelero
        
print ("len de rec =" , len(a))
print (a)
ang = np.radians(a)
#Temps = ((temps))
## ========================================================================
# cherche le nombre de tests et les intervales intéressants
# =========================================================================
#t=mesures(:,1)*1e-6 # temps en s
print ("=====passage au radioans")
#if typeTest==2 :
#    ang=np.radians(raquette); # angle en radian
#else : 
#    ang=np.radians(balle); # angle en radian
    
# fait la moyenne des n premieres valeurs pour trouver celle au repos
#ang0=mean(ang(t<(tpsRepos),1));
ang0 = (np.mean(ang[0:1000]))
print ("angle 0 = ", ang0)
# soustrait à toutes les autres pour mettre le repos sur angle=0
#ang=ang-ang0;
ang = np.subtract(ang, ang0)
# cherche les pics de la courbe d'oscillations libres amorties
#ind=peakfinder(ang,0.3);
#ind=[1;ind]; # rajoute l'indice 0 en début
plot_simple(ang)
from detect_peaks import detect_peaks
indexes = detect_peaks(ang, mph=0.15 , mpd=1000)

""" ancienne solution
import peakutils
indexes = peakutils.indexes(ang, thres=0.7, min_dist=1000)
"""
pic = np.array([])

print ("indexe des pics =")
print ((indexes))
print (len(indexes))
#indexes = indexes[1:]
indexes = np.append(0, indexes)
print ("fin")
for i in indexes:
    #print (i , ang[i])
    pic = np.append(pic, ang[i])
# cherche ceux pour lesquels il n'y a pas de pics juste avant
#indSub=ind(2:end)-ind(1:end-1); # difference des indices où les pics sont trouvés*
indSub = indexes[1:]-indexes[0:-1]

print ("indSub (difrence des indices = )", indSub)

time = np.arange(0.0, len(a), 1)

plt.figure(1)
plt.figure( figsize=(8, 6))

plt.plot(time , ang,'k', indexes, pic , 'ro',)

# pas de pics si le debut entre 2 max consecutifs est grand
#periode=indSub*(t(2,1)-t(1,1)); # diff des indices *  tps increment => periode en us

#indTest1erPic=ind(find((periode)>tpsEntreMesures)+1); # indices des débuts de tests
indTest1erpic = np.array([])
total_time=0;
for indexe in indSub:
    total_time +=indexe 
    if indexe > (tpsEntreMesures +1) *1000:
        print("kiki", (tpsEntreMesures +1) *1000, total_time)
        indTest1erpic = np.append(indTest1erpic , total_time)

nbTests=len(indTest1erpic) ; # nombre de tests
print("nombre de test = ", nbTests)
print ("premier pics a = ", indTest1erpic)


""
# ---
# figure
# plot(t,ang,'.');
# hold on
# plot(t(indTest1erPic),ang(indTest1erPic),'.r');


# je choisis de commencer l'acquisition non pas sur le premier pic, mais
# sur le premier 0 d'après pour éviter certains problèmes de conditions limites
#indMes0=find(ang<0);
print ("premier 0 dapres")
indMes0 =np.array([])
for counter, angl in enumerate(ang):
    if  ang[counter]  < 0 :
        #print (counter  + indexes[1], " ::  " , ang[counter + indexes[1]])
        indMes0 = np.append(indMes0 , counter)
        #break
#print ("indMEs0 = ", indMes0 ,indexes[1], ang[indMes0])

#indMes0 = np.array([])
#indMes0 = np.where(ang < 0)
print ("indMEs0 = ", len(indMes0) )
print ('indicateur =',indMes0[2])
#for ii=1:nbTests # pour chaque test, je cherche le 0 juste après le indTest1erPic
#    val=indMes0(find(indMes0>indTest1erPic(ii)));
#    ind0=[ind0;val(1)]; # la premiere valeur où ca croise l'axe est la valeur interessante
#end
ind0=np.array([])
val_ind0=np.array([])
for i in range (nbTests):
    val = np.array([])
    print ('pour le tesst n=', i, indTest1erpic[i])
    for indi in(indMes0):
#        print ('indi = ' , indi)
        if indi >  indTest1erpic[i]:
            val = np.append(val, indi)
            #print("find indi sup a indtest1piec", indTest1erpic[i])
    print ("val de 1 = ",val[1])
    ind0 = np.append(ind0, val[1])    
    val_ind0 = np.append(val_ind0 , ang[int(val[1])])

print ("ind0 ==" , ind0)


# regarde si la valeur juste avant le ind0 (soit juste avant l'axe n'est pas mieux)
#indChoix=[ind0,ind0-1];
#[aa,bb]=min(abs(ang(indChoix)),[],2);

#x.argmin(axis=0)

#ind0OK=[];
#for ii=1:nbTests
#    ind0OK=[ind0OK;indChoix(ii,bb(ii))];  # nouvelle valeur de ind0
#end
    
indTestFin = np.array([]) 
valTestFin = np.array([])
for n in range (nbTests):
    print ("for nb de tzest = ", nbTests)
    for i , ind in enumerate (indexes):
        if ind == indTest1erpic[n] : 
            print ("indice et indtest1pic ", ind , indTest1erpic[n])
            indTestFin = np.append(indTestFin , indexes[i + nbOscillations])
            valTestFin = np.append(valTestFin , ang[indexes[i + nbOscillations]])
    
"""
clear ii
ind0=ind0OK;
clear ind0OK

# je garde un certain nombre (nbOscillations) d'oscillations
# recupère les indice de fin indTestFin
try
    for n=1:nbTests
        indTestFin(n)=ind(find(ind==indTest1erPic(n))+nbOscillations);
    end
catch # si ca plante, c'est qu'il n'y a pas eu assez d'oscillations
    nbTests=nbTests-1; # j'enleve 1 sur nombre de tests, et je recommence
    indTest1erPic(end)=[];
    for n=1:nbTests
        indTestFin(n)=ind(find(ind==indTest1erPic(n))+nbOscillations);
    end
end



## ========================================================================
#  etudie chaque test indépendamment
# =========================================================================
indAEnlever=[]; # initialisation
mLI=zeros(nbTests,2);       # initialisation donnée de masse, longueur, et inertie
rendement=zeros(nbTests,1); # initialisation rendement
"""
indAEnlever = np.array([])
mLI = np.zeros((nbTests,2),)
rendement = np.zeros((nbTests,1))
#ang_test = [None] 
#t_test = [None]
for i in range (nbTests):
    print ("pour inbrdetest donne indice zero et fin = ", i, ind0[i] , indTestFin [i] )
    t_test = np.arange(int(ind0[i]), int(indTestFin[i]),0.001  )
#    Temps = Temps / 1000000
#    t_test = Temps[ int(ind0[i]): int(indTestFin[i]) ] 
    temp = ((int(indTestFin[i]) - int(ind0[i])) ) / 1000
#    temp = 
    temp_test = np.arange(0,temp , step = 1 / 1000 ) 
    ang_test =   - ang [ int (ind0[i]) : int (indTestFin [i]) ]
#-----nouvelle version avec temp reel
    t_test = Temps[ int(ind0[i]): int(indTestFin[i]) ] 
    temp =  (t_test[-1]-t_test[0])/1000
    t_test = np.asarray(t_test)
    temp_test= t_test - t_test[0]
    temp_test = temp_test 
#------fin nouvelleversion
#    ang_test =   ang [ int (ind0[i]) : int (indTestFin [i]) ]  
    print ("fin de la zone test")
    print ("len de time test = ", len(t_test) , t_test)
    print (" len de and_test =" , len (ang_test))
    print ("len de temp test = ", len(temp_test) , temp_test)
    print ()
    print ('le temps calculé est de = ' ,  temp)
    
    
    """
    for n=1:nbTests  # pour chaque test identifié
      
        # test si il y a bien des valeurs dedans (et donc s'il n'y a pas eu un
        # pb lors du choix des ind0 et indTestFin
        t_test=t(ind0(n):indTestFin(n));
        if isempty(t_test)
            continue  # je saute ce cas si il est empty
        end
             
        ang_test=ang(ind0(n):indTestFin(n));
        
        ang_test=-ang_test; # j'inverse pour avoir la premiere oscillation qui monte
            # bon dans tous les cas (balle ou raquette)
        
        # recale de telle manière que le premier temps soit=0
        t_test=t_test-t_test(1);
    """


    """   
    
    # -------- filtrage / simplification    
    # méthode moyenne flottante
    p=1/100; # coefficient
    angMoy=zeros(length(ang_test),1);
    angMoy(1)=ang_test(1);
    for i=2:length(ang_test)
        angMoy(i)=angMoy(i-1)*(1-p)+p*ang_test(i);
    end    
    
    # Trace les valeurs moyennées
    figure
    plot(t_test,angMoy,'.','color',rgb('lightblue'))  
    
    # je prends qu'une valeur sur nSimpl (pour simplifier au niveau des interpolations et ne pas tout prendre en compte
    # car le temps d'échantillonnage est très faible)
    nSimpl=20;
    t_test=t_test(1:nSimpl:end);
    ang_test=ang_test(1:nSimpl:end);
        
    
    # Trace les mesures (simplifiées par nSimpl)
    hold on
    plot(t_test,ang_test,'.') # en bleu marine
    xlabel('Temps (sec)')
    ylabel('Angle (rad)')    
    title(['Test de calibrage numero ' int2str(n)])
       
    # check si l'utilisateur est OK avec ce test
    testValid= 1;       # ajouté par Pierrick
    #     testValid=input('Doit obtenir de belles oscillations. OK ? 1=oui (0=non) ? ===> ');
    if testValid==0
        print('On oublie ce test et passe à celui d''après'); 
        indAEnlever=[indAEnlever,n];
        continue 
    end
    # ferme la figure (utiliser le close(gcf)) ou pas (alors tout se trace sur la même)
    # close(gcf);
    
    
    # -------- Calcul le rendement sur les maxima
    """
    
    ind_n =detect_peaks(ang_test, mph=0.3 , mpd=1000) ;
    print (" nex indexes = ", ind_n)
    pic_n = np.array([])
    for i in ind_n :
        pic_n = np.append ( pic_n ,ang_test[i])
    print (" et les pics = ", pic_n)
    
    #indexes = detect_peaks(ang, mph=0.3 , mpd=1000)
    rendement[n] = math.pow( ((1-math.cos(pic_n[-1])))/(1-math.cos(pic_n[1])) , 1/(4*(nbOscillations - 1)) )
    print ( " le rendemlent calculé = ", rendement[n])
    diffT = np.array([])
    diffT =  t_test[ind_n[2:]] -t_test[ind_n[1:-1]]
    print ("les peiodes sont = ", diffT)
    T = np.mean(diffT)
    print ( " le moyenne des periodes = ", T)

    """            
    # cherche les pics
    [ind,pic]=peakfinder(ang_test);
    # => rendement calculé pour 1/4 d'oscillations (de haut en bas 1 fois) en suivant une
    # interpolation lineaire sur les courbes d'oscillations) : les pertes sont approximées linéairement
    rendement(n)=((1-cos(pic(end)))/(1-cos(pic(1))))^(1/(4*(nbOscillations-1)));

    # périodes de chaque oscillation
    diffT=t_test(ind(2:end))-t_test(ind(1:end-1));
    T=mean(diffT);  # moyenne des périodes
    """

    """     
    
 	# -------- Calcul courbe de regression
    
    # valeurs initiales et bornes pour l'algorithme d'optimisation
    switch typeTest
        case 1 # balle
        #             mLinit=0.7;
            mLinit=mTigePendule*posCdGTigePendule+mBoule*posCdGBoule; # mL est imposé
            mL=mLinit;
            kv1init=0.01;
            Iinit=0.55; 
            cmax=[mL; 0.02; 0.6;   0.01;  3.5];         # bornes supérieures pour le vecteur [mL,kv1,I,pos0,vit0]
            cmin=[mL; 0;    0.5;  -0.01;  1];     
        case 2
        #             mLinit=0.6;     
            mLinit=mRaquette*1e-3*(posCdGRaquette*1e-3+posRaquette)+mBras*posCdGBras; # mL est imposé
            mL=mLinit;
            kv1init=0.035;                               
            Iinit=T^2*mLinit*9.81/(4*pi^2);      # T période moyenne          
            cmax=[mL; 0.04;    0.3;   0.01;  4];        # bornes supérieurs 
            cmin=[mL; 0.03;    0.2;  -0.01;  2];        # bornes inférieures
    end
    """
    if typeTest == 1 :
        mLinit=mTigePendule*(posCdGTigePendule+mBoule*posCdGBoule); # mL est imposé
        mL = mLinit 
        kvlinit = 0.01
        Iinit = 0.55
        cmin = [mL, 0,    0.4,  -0.01,  1]
        cmax = [mL, 0.02, 0.7  ,   0.01,  5]
        cmin = [mL, 0,    0.5,  -0.01,  1]
        cmax = [mL, 0.02, 0.6,   0.01,  3.5]
    if typeTest == 2 :
        mLinit=float(mRaquette)*1e-3*(float(posCdGRaquette)*1e-3+float(posRaquette))+float(mBras)*float(posCdGBras); # mL est imposé
        print ("mlinit calcule", mLinit)
        mL = mLinit 
        mL = 0.74
        kvlinit = 0.035
        Iinit=float(math.pow(T,2))*float(mLinit)*9.81/(4*math.pow(math.pi,2));      # T période moyenne          
        Iinit=0.30;      # T période moyenne          
        cmin = [mL, 0.03, 0.2, -0.01,  2]
        cmax = [mL, 0.04, 0.3,  0.01,  5]
#    posRaquette = 0.269; # raquettes de longueur 685mm
#    mBras = 4.892; # PB et GR le 06/02/18
#    posCdGBras = 0.179; # PB et GR le 06/02/18
#    mRaquette = 500 
#    posCdGRaquette = 25 
#    mLinit=mRaquette*1e-3*(posCdGRaquette*1e-3+posRaquette)+mBras*posCdGBras; # mL est imposé
#    mL = mLinit 
#    mL =0.238
#    Iinit = 1.05
#    print ("mLinit = " , mL)
#    kvlinit = 0.035
#    #☺kvlinit = 0.6
#    Iinit = math.sqrt (T) * mLinit * 9.81 / (4*math.sqrt ( math.pi) )
#    cmax = [ ml , 0.04 , 0.3 , 0.01 , 4]
#    cmin = [ ml , 0.03 , 0.2 , -0.01 , 2 ]
    #Iinit = 0.238
    pos0init=0;
    vit0init = (ang_test[10] - ang_test[1]) / (t_test[10]/1000 - t_test[1]/1000)
#    vit0init = vit0init * 10000
    print ("calcule vit0=", ang_test[10] , ang_test[1], t_test[10]/1000,t_test[1]/1000)
    
    vit0init = 3.05

    print ("vit0)init  = ", vit0init)  
#--------------------------------
#    def EDCoulomb(t,x,param):
#        g = 9.81
#        mL = param[0]
#        kv1 = param[1]
#        I = param[2]
#        
#        theta = x[0] # angle
#        dtheta = x[1] # vitesse
#        ddtheta=-mL * g / I * math.sin(theta) - kv1 / I * dtheta # acceleration
#
#        return (dtheta , ddtheta)
#    def Myfun2Minimize(c,donnees):
#        param = (mL , kvlinit , Iinit)
#        t0 = 0.0 
#        tmax = 12.0
#        from scipy import integrate
#        solED = integrate.solve_ivp(lambda t, y:EDCoulomb(t,y,param) , (t0, tmax) , (pos0init, vit0init)  , dense_output = True, max_step = 0.0008, rtol = 1e-8,atol = 1e-10)    
##        print( solED)
#        fxi  = solED.y[0][:len (ang_test)]
#        res  = np.sum(fxi-ang_test)*math.ceil(math.log2(abs(2)))  
#        print ("res =¨" , res)
#        return res
#    """
#        # utilisation de fminsearchbnd
#    options=optimset('fminsearch');
#    options.TolFun=1e-2; options.TolX=1e-2; options.TolCon=1e-2; # options de convergence
#    [xsol,fval]=fminsearchbnd(@(x)Myfun2Minimize(x,[t_test,ang_test]),[mLinit,kv1init,Iinit,pos0init,vit0init],cmin,cmax,options);
#    # si je ne veux pas lancer l'optimisation mais utiliser des valeurs approchées
#    """
#
#    from scipy import optimize
#    
##    (xsol , fval) = optimize.fmin(lambda x : Myfun2Minimize ( x, [t_test , ang_test]) ,[mLinit,kvlinit,Iinit,pos0init,vit0init], ftol = 1e-2 , xtol = 1e-2  )    
#    xsol = optimize.fmin(lambda x : Myfun2Minimize ( x, [t_test , ang_test]) ,[mLinit,kvlinit,Iinit,pos0init,vit0init], ftol = 1e-2 , xtol = 1e-2  )    
#    print ("totut = ", xsol) 
#------------------------------------------
    def EDCoulomb(t,x,param):
        g = 9.81
        mL = param[0]
        kv1 = param[1]
        I = param[2]
        
        theta = x[0] # angle
        dtheta = x[1] # vitesse
        ddtheta=-mL * g / I * math.sin(theta) - kv1 / I * dtheta # acceleration

        return (dtheta , ddtheta)
    def Myfun2Minimize(c,donnees):
#        param = (mL , kvlinit , Iinit)
#        print ("donnes =" , c)
        param = (c[0], c[1] , c[2])
        t0 = 0.0 
        tmax = 14.0
        from scipy import integrate
        #from scipy.integrate import odeint
        solED = integrate.solve_ivp(lambda t, y:EDCoulomb(t,y,param) , (t0, tmax) , (c[3], c[4])  , max_step = 0.001,dense_output=True, rtol = 1e-5,atol = 1e-8)    
#        print( solED)
#        fxi  = solED.y[0][:len (ang_test)]
#        res  = np.sum(fxi-ang_test)*math.ceil(math.log2(abs(2)))  
        from scipy.integrate import odeint
#        XsolED = odeint(solED, 0, len (donnees[1]))
        
        fxi  = solED.y[0][:len (donnees[1])]
#        res  = np.sum(fxi-donnees[1])*math.ceil(math.log2(abs(2))) 
        res = 0
        for i in range(len (fxi)):
            res += math.pow((fxi[i]-donnees[1][i]) , 2)
        print ("res =¨" , res)
        return res
    """
        # utilisation de fminsearchbnd
    options=optimset('fminsearch');
    options.TolFun=1e-2; options.TolX=1e-2; options.TolCon=1e-2; # options de convergence
    [xsol,fval]=fminsearchbnd(@(x)Myfun2Minimize(x,[t_test,ang_test]),[mLinit,kv1init,Iinit,pos0init,vit0init],cmin,cmax,options);
    # si je ne veux pas lancer l'optimisation mais utiliser des valeurs approchées
    """
    from scipy.optimize import Bounds
#    bounds = Bounds([mL, 0.03,    0.2,  -0.01,  2], [mL, 0.04,    0.3,   0.01,  4])
    bounds = Bounds(cmin, cmax)

    from scipy import optimize
    
#    (xsol , fval) = optimize.fmin(lambda x : Myfun2Minimize ( x, [t_test , ang_test]) ,[mLinit,kvlinit,Iinit,pos0init,vit0init], ftol = 1e-2 , xtol = 1e-2  )    
#    xsol = optimize.fmin(lambda x : Myfun2Minimize ( x, [t_test , ang_test]) ,[mLinit,kvlinit,Iinit,pos0init,vit0init] )    
#    xsol = optimize.minimize(lambda x : Myfun2Minimize ( x, [t_test , ang_test]) ,[mLinit,kvlinit,Iinit,pos0init,vit0init] , bounds=bounds, tol = 1e-2  )
    xsol = optimize.minimize(lambda x : Myfun2Minimize ( x, [t_test , ang_test]) ,[mLinit,kvlinit,Iinit,pos0init,vit0init] , method='Nelder-Mead' , tol = 1e-2, bounds=bounds,options={ 'xatol': 1e-2, 'fatol': 1e-2}  )    
#    xsol = optimize.minimize(lambda x : Myfun2Minimize ( x, [temp_test , ang_test]) ,[mLinit,kvlinit,Iinit,pos0init,vit0init],  method='L-BFGS-B', bounds=bounds,options={ 'xatol': 1e-2, 'ftol': 1e-2}  )
    
    
    print ("totut = ", xsol) 
    print ("totut de x = ", xsol.x) 

    
    
    
    """
    pos0init=0;
    vit0init=(ang_test(10)-ang_test(1))/(t_test(10)-t_test(1)) ; # trouvé en utilisant le graphe

    # utilisation de fminsearchbnd
    options=optimset('fminsearch');
    options.TolFun=1e-2; options.TolX=1e-2; options.TolCon=1e-2; # options de convergence
    [xsol,fval]=fminsearchbnd(@(x)Myfun2Minimize(x,[t_test,ang_test]),[mLinit,kv1init,Iinit,pos0init,vit0init],cmin,cmax,options);
    # si je ne veux pas lancer l'optimisation mais utiliser des valeurs approchées
#     xsol=[mL,kv1init,Iinit,pos0init,vit0init];

    # les coefficients de la regression moindre carreees sont dans xsol
    Ifmincon=xsol(3);
    mLfmincon=xsol(1);
    c=xsol ; 
    """
#    Ifmincon=xsol[2]
#    print ("Ifmincom =" , Ifmincon)
#    mLfmincon=xsol[0]
#    c=xsol 
#    paramR=(c[0:3])
#    pos0R=c[3]
#    vit0R=c[4]
    Ifmincon=xsol.x[2]
    print ("Ifmincom =" , Ifmincon)
    mLfmincon=xsol.x[0]
    c=xsol.x
    paramR=(c[0:3])
    pos0R=c[3]
    vit0R=c[4]
    
    print ("paramR = " , paramR)
    from scipy import integrate
    solEDR = integrate.solve_ivp(lambda t, y:EDCoulomb(t,y,paramR) , (0, temp) , (pos0R, vit0R)  , dense_output = True, max_step = 0.001, rtol = 1e-8,atol = 1e-10)    

    time = np.arange(0.0, len(a), 1)
    
    plt.figure(1)
    plt.figure( figsize=(8, 16))
    #plt.gcf().subplots_adjust(wspace = 0, hspace = 4)
    
    plt.subplot(611)
    plt.ylabel('teta')
    plt.plot(time ,  a,'k' )
    
    plt.subplot(612)
    plt.ylabel('teta')
    plt.plot(time , ang,'k', indexes, pic , 'ro', ind0 , val_ind0 ,'bo' , indTestFin , valTestFin , 'go' )
    
    plt.subplot(613)
    plt.ylabel('teta')
    plt.plot( temp_test ,  ang_test ,'k', temp_test[ind_n], pic_n ,'ro' )
#    plt.plot( temp_test ,  ang_test ,'k' )
    
    plt.subplot(614)
    plt.ylabel('teta')
    plt.plot( temp_test ,  ang_test ,'k' )
    
    plt.subplot(615)
    plt.ylabel('accel')
    plt.plot( solEDR.t ,  solEDR.y[0] ,'k' )
    
    plt.subplot(616)
    plt.ylabel('les deux')
#    plt.plot( solEDR.t ,  solEDR.y[1] ,'k' )
    plt.plot( temp_test/1000 ,  ang_test ,'k' )
    plt.plot( solEDR.t ,  solEDR.y[0] ,'b' )
    plt.show()

    """
    # ---- pour tracer
    # refait le calcul une derniere fois en utilisant les constantes 'c' pour trouver (et tracer) la fonction de regression 
    paramR=c(1:end-2);
    pos0R=c(end-1);
    vit0R=c(end);
    solEDR=ode45(@(t,y)EDCoulomb(t,y,paramR),[0 20],[pos0R vit0R],optionsOde);
    XR=deval(solEDR,t_test)';
    yRegressMC = (XR(:,1)) ; # degres
    hold on
    plot(t_test,yRegressMC,'.c');      # trace la regression en cyan

    mLI(n,:)=[mLfmincon,Ifmincon]; # garde le couple (mL, I) du fmincon
    
    
end

# enleve les lignes pour lesquelles on a dit que la courbe n'était pas OK
mLI(indAEnlever,:)=[];
rendement(indAEnlever,:)=[];


## ========================================================================
#  post-processing, enregistrement
# =========================================================================
# enregistre ces valeurs en dur dans le workspace
switch typeTest

        case 1 # balle
            
        mLI_balle=mLI;
        if exist('mLI_balle.mat','file')  # recupere précédentes valeurs
            mLI_avant=load('mLI_balle.mat');
            print('*******************************************************')
            print('les valeurs caractéristiques du dernier calibrage balle étaient ')
            [mLI_avant.mLI_balle]
        end
            print('*******************************************************')
        print('les nouvelles valeurs du calibrage balle sont:')
        [mLI_balle]
            print('*******************************************************')


        rendement_balle=rendement;
        if exist('rendement_balle.mat','file')  # recupere précédentes valeurs
            rendement_avant=load('rendement_balle.mat');
            print('*******************************************************')
            print('les valeurs de rendement du dernier calibrage balle étaient ')
            [rendement_avant.rendement_balle]
        end
            print('*******************************************************')
        print('les nouvelles valeurs du rendement calibrage balle sont:')
        rendement_balle
            print('*******************************************************')
        # testSauve=input('Voulez-vous les sauvegarder ? (1 si OUI) ');
        testSauve=1;
        
        switch testSauve
            case 1          
                save('rendement_balle','rendement_balle');
                save('mLI_balle','mLI_balle');
            otherwise
        end
        

    case 2 # raquette
        
        mLI_raquette=mLI;
        
        if exist('mLI_raquette.mat','file')  # recupere précédentes valeurs
            mLI_avant=load('mLI_raquette.mat');
            print('*******************************************************')
            print('les valeurs caractéristiques du dernier calibrage raquette étaient ')
            [mLI_avant.mLI_raquette]
        end
        print('*******************************************************')
        print('les nouvelles valeurs du calibrage raquette sont:')
        [mLI_raquette]
         print('*******************************************************')

            
        rendement_raquette=rendement;
        if exist('rendement_raquette.mat','file')  # recupere precedentes valeurs
            rendement_avant=load('rendement_raquette.mat');
            print('*******************************************************')
            print('les rendements du dernier calibrage raquette étaient ')
            print('i''eme ligne pour test calibrage numero i');
            [rendement_avant.rendement_raquette]
        end
        print('*******************************************************')
        print('les nouvelles valeurs de rendement du calibrage raquette sont:')
        print('i''eme ligne pour test calibrage numero i');
        [rendement_raquette]
        print('*******************************************************')
        # testSauve=input('Voulez-vous les sauvegarder ? (1 si OUI) ');
        testSauve=1;
        
        switch testSauve
            case 1
                save('rendement_raquette','rendement_raquette');
                save('mLI_raquette','mLI_raquette');
                # enregistre aussi les données de raquette (masse + posCdG)
                save('donneesGeomRaquette','mRaquette','posCdGRaquette');
            otherwise
        end
        
end


# le calibrage de la balle et de la raquette donnent en sortie le couple
 mLI = [mL , I]
# avec m masse totale
# L position du centre de gravite
# I moment d'inertie

# l'energie à l'impact est ensuite calculée en utilisant ces valeurs de mL et I

# ils donnent aussi les rendements des differents tests dans rendement_raquette.mat et rendement_balle.mat



#Ajouté par Pierrick
# testSauve=input('Voulez-vous sauvegarder toutes les datas de cette calibration ? (oui => 1 / non par default) ');
testSauve=1;
cd(current.data) # Je change mon current directory
        switch testSauve
            case 1
                save(['Calibration_',nomRaquette])
            otherwise # par default, on n'enregistre rien
        end
cd(current.script) # je restaure mon ancien current directory
"""
## ========================================================================
#  post-processing, enregistrement
# =========================================================================
# enregistre ces valeurs en dur dans le workspace
print ("rendement = " , rendement)
print ("mLI =" , mLI)
print ("file tot write")

import os.path
print ("typeTest=",typeTest)
if typeTest == 1 : #balle
    print ("type=test =",typeTest)
    if os.path.isfile('mLI_balle.py'):# recupere précédentes valeurs
        from mLI_balle import mLI_b
        print('*******************************************************')
        print('les valeurs caractéristiques du dernier calibrage balle étaient ', mLI_b)
    print('*******************************************************')
    print('les nouvelles valeurs du calibrage balle sont:' , mLI)
    print('*******************************************************')   
    if os.path.isfile('rendement_balle.py'):
        from rendement_balle import rendement_b
        print('*******************************************************')
        print('les valeurs caractéristiques du dernier rendement balle étaient ', rendement_b)
    print('*******************************************************')
    print('les nouvelles valeurs du rendement balle sont:' , rendement)
    print('*******************************************************')  
    
    testSauve=input('Voulez-vous les sauvegarder ? (1 si OUI) ');
    testSauve=1;
    if testSauve == 1 :
        file_to_write = open("mLI_balle.py","w+")
        file_to_write.write("mLI_b = " + repr(mLI) + "\n")
        file_to_write.flush()
        file_to_write.close()
        file_to_write = open("rendement_balle.py","w+")
#        file_to_write.write("rendement_b = " + repr(rendement) + "\n")
        file_to_write.write("rendement_b = " + (rendement) + "\n")
        file_to_write.close()
if typeTest == 2 : #raquette
    print ("type=test =",typeTest)
    if os.path.isfile('mLI_raquette.py'):# recupere précédentes valeurs
        from mLI_raquette import mLI_r
        print('*******************************************************')
        print('les valeurs caractéristiques du dernier calibrage raquette étaient ', mLI_r)
    print('*******************************************************')
    print('les nouvelles valeurs du calibrage requette sont:' , mLI)
    print('*******************************************************')   
    if os.path.isfile('rendement_raquette.py'):
        from rendement_raquette import rendement_r
        print('*******************************************************')
        print('les valeurs caractéristiques du dernier rendement raquette étaient ', rendement_r)
    print('*******************************************************')
    print('les nouvelles valeurs du rendement raquette sont:' , rendement)
    print('*******************************************************')  
    
    testSauve=input('Voulez-vous les sauvegarder ? (1 si OUI) ');
    testSauve=1;
    if testSauve == 1 :
        file_to_write = open("mLI_raquette.py","w+")
        file_to_write.write("mLI_r = " + repr(mLI) + "\n")
        file_to_write.flush()
        file_to_write.close()
        file_to_write = open("rendement_raquette.py","w+")
        file_to_write.write("rendement_r = " + repr(rendement) + "\n")
        file_to_write.close()
"""   
end

# enleve les lignes pour lesquelles on a dit que la courbe n'était pas OK
mLI(indAEnlever,:)=[];
rendement(indAEnlever,:)=[];





# le calibrage de la balle et de la raquette donnent en sortie le couple
# mLI = [mL , I]
# avec m masse totale
# L position du centre de gravite
# I moment d'inertie

# l'energie à l'impact est ensuite calculée en utilisant ces valeurs de mL et I

# ils donnent aussi les rendements des differents tests dans rendement_raquette.mat et rendement_balle.mat



#Ajouté par Pierrick
# testSauve=input('Voulez-vous sauvegarder toutes les datas de cette calibration ? (oui => 1 / non par default) ');
testSauve=1;
cd(current.data) # Je change mon current directory
        switch testSauve
            case 1
                save(['Calibration_',nomRaquette])
            otherwise # par default, on n'enregistre rien
        end
cd(current.script) # je restaure mon ancien current directory
"""
print ("size = ", a.size)
#time = np.arange(0.0, 64058, 1)
time = np.arange(0.0, len(ang), 1)

plt.figure(1)
plt.figure( figsize=(8, 16))
#plt.gcf().subplots_adjust(wspace = 0, hspace = 4)

plt.subplot(311)
plt.ylabel('teta')
plt.plot(time ,  a,'k' )

plt.subplot(312)
plt.ylabel('teta')
plt.plot(time , ang,'k', indexes, pic , 'ro', ind0 , val_ind0 ,'bo' , indTestFin , valTestFin , 'go' )

plt.subplot(313)
plt.ylabel('teta')
#plt.plot( temp_test ,  ang_test ,'k', ind_n, pic_n ,'ro' )
plt.plot( temp_test ,  ang_test ,'k' )
plt.plot( solEDR.t ,  solEDR.y[0] ,'b' )


plt.show()