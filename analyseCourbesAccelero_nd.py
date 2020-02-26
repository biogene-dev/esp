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

# fait l'analyse des courbes provenant de l'accéléromètre monté sur la raquette
# permet de récupérer  :
# - la mesure d'amplitude pic à pic
# - la mesure de l'énergie totale de l'impact
# - l'énergie pour une bande de fréquence donnée


# ------ mesures provenant de 'Acquisition'
import scipy.io
# open file
from pathlib import Path
import numpy as np
import os.path
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
            val = float(to_test[4])
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
                transAng=(val); # le 0.8 est du à la bande de 10# -> 90# du capteur
#            transAng=(val*360/((2^13)*0.8)); # le 0.8 est du à la bande de 10# -> 90# du capteur
                if transAng < 10000 : a =np.append(a, transAng)
                if int(to_test[0])< 9999999999 : t= np.append(t ,  (float(to_test[0])/1000000) ) 

#            print (len(a))
                
mesureAcc=(a*100/8192)*1.66;


import math
# sphinx_gallery_thumbnail_number = 3
import matplotlib.pyplot as plt

plt.figure( figsize=(8, 6))
#plt.gcf().subplots_adjust(wspace = 0, hspace = 4)
#    plot(t(indTestDebutB(ii):indTestFin(ii)),rad2deg(ang(indTestDebutB(ii):indTestFin(ii))),'.r');
t=t-t[0]
acc0 = np.mean(mesureAcc[0:100])
acc= mesureAcc-acc0
plt.ylabel('teta')
plt.plot(t,acc)
plt.show()

from detect_peaks import detect_peaks
indexes = detect_peaks(acc, mph=10 , mpd=1000)
print ("indexes =", indexes)
#
## PB et CH 27/10/2014
#xe = 1:length(mesureAcc);
#y = mesureAcc(:,2);
## xi = 1:length(mesureAcc)/(length(mesureAcc)/25):length(mesureAcc);
#xi = 1:50000/2000:length(mesureAcc);
#yi = interp1(xe,y,xi,'spline')';
#
## t=mesureAcc(:,1)*1e-6; # en sec
#t=mesureAcc(1:25:length(mesureAcc(:,1)),1)*1e-6; # PB 27/10/2014
## acc=mesureAcc(:,2); # en g
#acc = yi; #PB 27/10/2014
#
#n=length(acc); # nb de points
#feGl=1/(t(2)-t(1));    # freq d'echantillonage en Hz
#TeGl=1/feGl; # periode d'echantillonage en s
#acc0=mean(acc(t<tpsRepos,1));
## soustrait pour mettre à 0 au repos
#acc=acc-acc0;
#
## enleve les outliers (fonction effectuées plusieurs fois)
## for ii=1:6
##     acc=RemoveOutlier(acc,20);
## end
#
#
## trace la courbe resultat
#hf=figure;
#    hold on
#plot(t,acc,'.-') # bleu
#
## =====================================
## cherche le nombre d'impacts
## =====================================
#

## ind=peakfinder(acc,20,80); # cherche les peaks
## 
## indloc=zeros(length(ind),1);
## for ii=2:length(ind)
##     if ind(ii)-ind(ii-1)<tpsEntreMesures/TeGl # alors c'est le meme choc
##         indloc(ii)=1;
##     end
## end
## # enlève les pics qui ne servent à rien car proviennent du même impact
## ind(find(indloc))=[];
#
## si on connait la période d'échantillonage et que celle ci est 2 fois plus
## importante que pour les capteurs angulaires (regarder dans Acquisition),
## alors je peux utiliser directement l'indice trouvé pour chaque impact
## dans analyseCourbesImpactBalle
#
#ind=2*fix(indTestDebutB./25); # PB 27/10/2014 a ajouté './25'
#
#nChoc=length(ind); # nb de chocs
#
###
## =====================================
## zoom sur chaque impact et étude de chaque impact
## =====================================
#
#tAvantChoc=0.02;    # tps en sec pris en compte avant l'impact pour la fft
#tApresChoc=0.4-tAvantChoc;       # tps en sec prix en compte apres l'impact pour la fft # il y avait 0.8s
#
tAvantChoc=200;    # tps en sec pris en compte avant l'impact pour la fft
tApresChoc=200;       # tps en sec prix en compte apres l'impact pour la fft # il y avait 0.8s
max_i = []
min_i = []
amplitude = []
valEnergieTempTout = []
valEnergieFreqTout = []
valEnergieFreqBande = []
for i in range (len(indexes)):
    plt.ylabel('accele')
    plt.plot(t[indexes[i]-100:indexes[i]+200],acc[indexes[i]-100:indexes[i]+200])
    plt.show()
    from numpy.fft import fft
    T= t[indexes[i]+100] - t[indexes[i]-10] 
    fe = 1000
    nfft = np.arange(start=0.0,stop=4,step=1.0/fe)
    
    echantillons = acc[indexes[i]-100:indexes[i]+200]
    tfd = fftn(echantillons)
    N=len(echantillons)
    spectre = np.absolute(tfd)*2/N
    
    freq=np.arange(N)*1.0/T
    
    plt.figure(figsize=(10,4))
    plt.plot(freq,spectre,'r')
    plt.xlabel('f')
    plt.ylabel('A')
    plt.axis([-0.1,fe/2,0,spectre.max()])
    plt.grid()
    plt.show()
    mx = abs(tfd)/len(echantillons); # magnitude de FFT 
#    print ("la magnitude est de =",mx)
#    mx = math.pow(mx,2); # amplitude (magnitude^2)
    amplitude.append( max(echantillons) - min(echantillons) )
    print ("amplitude pic ap ic ", amplitude)
    valEnergieTempTout.append ( np.trapz((np.power(echantillons,2)), t[indexes[i]-100 : indexes[i]+200 ]) )
    print ("energie temmp tout = ", valEnergieTempTout)
    valEnergieFreqTout.append ( np.trapz(spectre) )
    print ("energie freq tout = ", valEnergieFreqTout)
    valEnergieFreqBande.append ( np.trapz(spectre[60:80]) )
    print ("energie freq Bande 60-80  = ", valEnergieFreqBande)
#for ii=1:nChoc
#
#    range=ind(ii)-tAvantChoc/TeGl:ind(ii)+tApresChoc/TeGl;
#    
#    # si le nb de valeurs dans range est impair, je rajoute 1 (nécessaire pour la fft)
#    if isodd(length(range))
#        range(end+1)=range(end)+1;
#    end
#    
#    # plot la fourchette utile en couleur rouge (pour la démarquer)
#    figure(hf);
#    plot(t(range),acc(range),'.r') # rouge
#
#    # ------- FFT (Fast Fourier Transform)
#    x=acc(range);
#    nfft= 2^(nextpow2(length(x))); 
#    fftx = fft(x,nfft);  # fonction FFT de matlab
#    NumUniquePts = ceil((nfft+1)/2); # nb de points unique (la moitié)
#    fftx = fftx(1:NumUniquePts); # garde que la permière moitié (FFT est symmétrique)
#    mx = abs(fftx)/length(x); # magnitude de FFT 
#    mx = mx.^2; # amplitude (magnitude^2)
#    # Comme on a laché la moitié de la FFT, on multiplie par 2 pour garder la même énergie (sauf pour le Nyquist point) 
#    if rem(nfft, 2) # odd nfft n'a pas de Nyquist point
#      mx(2:end) = mx(2:end)*2;
#    else
#      mx(2:end -1) = mx(2:end -1)*2;
#    end
#    fbon = (0:NumUniquePts-1)*feGl/nfft; 
#
#    # ----- Trace spectre de fréquence
#    figure(21)
#    plot(fbon,mx,'.-'); 


    
#
#    valEnergieTempTout(ii)=trapz(t(range),x.^2);                         # énergie du signal dans le domaine temporel (intégrale de x^2)
#    valEnergieFreqTout(ii)=trapz(fbon,mx)*(tApresChoc+tAvantChoc)^2;    # énergie du signal dans le domaine fréquentiel (intégrale de f * TpsTotal^2)
#    # ces 2 énergies doivent être égales
#
#    # ------- filtre le signal avec une "moving average"  => enlève le bruit
##     # filtre sur les m derniers points
##     m=20;
##     yMov=moving(x,m);
#     
#    #CH 28/10/2014
#    Wn = [10 400]/(feGl/2);
#    [b,a] = butter(4,Wn,'bandpass');
#    yMov = filtfilt(b,a,x);
#    
#    #PB 27/10/2014
##     [b,a] = butter(4,(20/(feGl/2)),'high');
##     yMov = filtfilt(b,a,x);
#     
#    #PB 27/10/2014
##     yMov = sgolayfilt(x,3,11); # Ordre 3 et avec une fenetre de 11 points
#
#    # fft sur le signal filtré
#    fftx=fft(yMov,nfft);
#    fftx = fftx(1:NumUniquePts); 
#    mx = abs(fftx)/length(x); 
#    mx = mx.^2; 
#    if rem(nfft, 2) 
#      mx(2:end) = mx(2:end)*2;
#    else
#      mx(2:end -1) = mx(2:end -1)*2;
#    end
#    fMov = (0:NumUniquePts-1)*feGl/nfft; 
#
#
#    # ----- Trace spectre de fréquence : signal filtré
#    figure(21)
#    hold on
#    plot(fMov,mx,'.-c');  # cyan
#    figure(hf);
#    hold on
#    plot(t(range),yMov,'.c')
#
#
#    valEnergieTempFilt(ii)=trapz(t(range),(abs(yMov)).^2); # énergie signal filtré domaine temporel
#    valEnergieFreqFilt(ii)=trapz(fMov,mx)*(tApresChoc+tAvantChoc)^2;; # énergie signal filtré domaine fréquentiel
#
#
#    # ------- filtre pour une bande passante particulière => sélectionne une seule bande passante
#    # fréquences entre fmin et fmax
#    fmax=80; # Hz
#    fmin=60; #Hz
#
#    y=1/feGl*fft(x);
#    y=fftshift(y);
#    Yfilter = [linspace(0,feGl,length(x))]';
#    Yfilter = and((abs(Yfilter-feGl/2)<=fmax),(abs(Yfilter-feGl/2)>=fmin)); # garde que ceux entre les 2 fréquences
#    Yfilter = Yfilter.*y;
#    Yfilter = fftshift(Yfilter);
#    sig_filter=feGl*ifft(Yfilter);
#
#    # ------ trace le signal filtré sur une bande passante (domaine temporel)
##     figure(hf);
##     plot(t(range),real(sig_filter),'.-m')  # magenta
#
#    # ------- nouvelle fft sur le signal filtré
#    fftx = fft(sig_filter,nfft); 
#    fftx = fftx(1:NumUniquePts); 
#    mx = abs(fftx)/length(x); 
#    mx = mx.^2; 
#    if rem(nfft, 2) 
#      mx(2:end) = mx(2:end)*2;
#    else
#      mx(2:end -1) = mx(2:end -1)*2;
#    end
#    fbon = (0:NumUniquePts-1)*feGl/nfft; 
#
#    # ------ trace le signal filtré sur une bande passante (domaine fréquentiel)
#    figure(21)
#    hold on
#    plot(fbon,mx,'.-m');
#
#    valEnergieFreqBande(ii)=trapz(fbon,mx)*(tApresChoc+tAvantChoc)^2; # Energie pour cette bande passante
#    
#    # ------- mesure de l'amplitude pic à pic
#    # utilise les 0.035 premières secondes apres le début (en fait apres 0.0035s du début) de l'impact pour
#    # trouver les max et min
##     [valMax(ii),indMax(ii)]=max(acc(range(round(0.0035/TeGl):round(0.035/TeGl)))); # valeurs max
##     [valMin(ii),indMin(ii)]=min(acc(range(round(0.0035/TeGl):round(0.035/TeGl)))); # valeurs min
#    # amplitude pic à pic calculée sur le signal filtré
#    [valMax(ii),indMax(ii)]=max(yMov);
#    [valMin(ii),indMin(ii)]=min(yMov);
#
#    amplitudePP(ii)=valMax(ii)-valMin(ii);    # amplitude pic à pic
#
#
#    
#    # différence entre les énergies calculées (doivent être ~0)
#    diffLoc=2*(valEnergieFreqTout(ii)-valEnergieTempTout(ii))/(valEnergieFreqTout(ii)+valEnergieTempTout(ii))*100;
#    disp(['Différence entre énergie domaine frequentiel et énergie domaine temporel pour signal ENTIER (#) : ' ...
#        num2str(diffLoc) ]);
#    # erreur si cette différence est trop grande (>2)
#    if abs(diffLoc)>3
#        valEnergieFreqTout(ii)=-1;
#        valEnergieFreqBande(ii)=-1;
#    end
#    
##     disp(['Différence entre énergie domaine frequentiel et énergie domaine temporel pour signal FILTRE (#) : ' ...
##         num2str(2*(valEnergieFreqFilt(ii)-valEnergieTempFilt(ii))/(valEnergieFreqFilt(ii)+valEnergieTempFilt(ii))*100) ]);    
#    
#    clear f N x nfft fftx yMov  # supprime certaines valeurs 
#    
#end
#
## ----- labels figure
#load('nomRaquette.mat');
#
#figure(hf)
#title(['Accélération (g) dans le domaine temporel - Raquette ' nomRaquette]) # titre pour le graphe
#legend('Courbe résultat originale','Domaine intéressant','Signal filtré hautes fréquences','Signal filtré sur bande passante')
#xlabel('temps (secondes)')
#ylabel('accélération (g)')
#figure(21)
#title(['Spectre de puissance (domaine fréquentiel) - Raquette ' nomRaquette]) # titre pour le graphe
#legend('Signal original','Signal filtré hautes fréquences','Signal filtré sur bande passante')
#xlabel('fréquence (Hz)')
#ylabel('puissance')
#
#
#disp(['============> Les valeurs recherchées pour ' nomRaquette ' sont '])
#fprintf('Amplitude pic à pic (g) : #s\n',num2str(amplitudePP,3)) # ;
#fprintf('Energie (J) totale de l''impact : #s\n',num2str(valEnergieFreqTout,3));
#fprintf('Energie (J) pour la bande de fréquence donnée #g - #g Hz : #s\n',fmin,fmax,num2str(valEnergieFreqBande,3))
#
## =====================================
## sauvegarde
## =====================================
## possibilité de sauvegarder en dur les résultats
#testSauve= 1;       # ajouté par Pierrick
## testSauve=input('Voulez-vous sauvegarder les résultats provenant de l''accéléromètre ? (oui => default / non => 0) ');
#switch testSauve
#    case 0
#        # rien
#    otherwise
#        fid=fopen('acceleroResultats.csv','a');
#        maDate=datestr(now);
#        fprintf(fid,'#s\n',maDate);
#        fprintf(fid,'#s\n',nomRaquette);
#        str='';
#        for ii=1:length(amplitudePP)
#            str=[str '#8.4f '];
#        end
#        str=[str '\n'];
#        str2=str; str3=str; # même format d'écriture pour avoir les colonnes alignées
#        fprintf(fid,str,amplitudePP'); # amplitude
#        fprintf(fid,str2,valEnergieFreqTout'); # energie totale
#        fprintf(fid,str3,valEnergieFreqBande'); # energie bande passante
#        
#        fclose(fid);
#
#end
