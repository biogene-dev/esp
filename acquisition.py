# -*- coding: utf-8 -*-
"""
Created on Wed May  8 19:48:22 2019

@author: ash
"""
import initialisation

global nomPortCOM , nomRaquette

if args<2 :
    typeTest=3; # test impact par default si typeTest non spécifié
end
    
# constantes pour connection port COM
# -------------
# initialisiation
brate=12e6; # baudrate
sizeOut=13; # j'envoie 13 valeurs uint8                
sizeIn=1024; #
timeOut=2;
typeAccelero=8; # = 4.4 ou 8 : specificité accélérometre
print([' ***** ------ ****** type d''accéléro :  ' (typeAccelero)]);

# ------------------------------------------------------------
# avant de commencer reset la carte: deconnection
# ------------------------------------------------------------
deconnection;

# ------------------------------------------------------------
# connection carte/PC
# ------------------------------------------------------------

# -------------
c=serial(nomPortCOM,'BaudRate',brate,'OutputBufferSize',sizeOut,'InputBufferSize',sizeIn, ...
        'Timeout',timeOut,'FlowControl','hardware','ByteOrder','littleEndian');   #  

fopen(c);
if strcmp(c.status,'open')
    print(['connecté sur le port ' nomPortCOM ]);
else
    error('problème de connection');
end
    
# ------------------------------------------------------------
## ecrire commande sur la carte + lire mesures
# ------------------------------------------------------------

commandInit=uint8([hex2dec('FF') hex2dec('00') hex2dec('3C') hex2dec('71')]);

# position des capteurs sur la carte
bAcc=uint8(00); # bit accelero
bAng1=uint8(02); # bit capteur angulaire 1
bAng2=uint8(03); # bit capteur angulaire 2
bNull=uint8(255); # ou n'importe quelle valeur >=10

# ordre dans lequel sont lus les données  => doit être de longueur 8
# ceci influence le temps d'echantillonnage
ordreLecture=uint8([bAng1 bAcc bAng2 bAcc bAng1 bAcc bAng2 bAcc]);

# timer en us   => à mettre en format char
tempo=uint8(10); # doit rentrer dans un uint8 (donc le max est 255us)

commandSend=[commandInit ordreLecture tempo]; # creer une ligne de commande à envoyer à la carte

# test commande
if length(commandSend)~=13; error('longueur de la commande pas bonne'); end

# ------------------------------------------------------------
# send data : envoie la commande petit à petit
for i=1:13
    fwrite(c,commandSend(i),'uint8');
    pause(0.01);
end

# ------------------------------------------------------------
# receive data : acquisition proprement dite
i=0;

# initialisation
clear dataRec
clear loc
resultat=zeros(1,tpsAcquisition*1e6/10); # initialisation (taille fonction du temps ~10s pour 1e6)

tInit=tic;
tFin=toc(tInit);

# waitbar init
strWbar='Acquisition des données';
multiWaitbar(strWbar,0,'CanCancel','on');

# boulce de lecture sur un temps donné
# si on veut à l'infini, faire un while 1
while tFin<tpsAcquisition
    
    # recup un vecteur de données de sizeIn
    dataRec=fread(c,sizeIn,'uint8');   

    j=1; # initialisation
    try # essai de lecture sur la carte (ne passe pas si le dataRec n'existe pas)
        
    loc=zeros(1,sizeIn/2); # initialisation
    for k=1:2:sizeIn # prend un bit sur 2 pour avoir le poids fort/ poids faible
        loc(j)=dataRec(k)*256+dataRec(k+1); # transformation de la valeur
        j=j+1;
    end
         
    # save
    resultat(i*sizeIn/2+1:(i+1)*sizeIn/2)=loc;

    
    catch # si il n'est pas arrivé à lire, c'est probablement car il plante sur le dataRec(k)
        str=[lastwarn '  numero ' int2str(i)];
        print(str);
        deconnection # je me deconnecte (toujours !!!)
        break

    end
    
    i=i+1;
    
    tFin=toc(tInit); # nouveau temps à la fin de la boucle (pour voir si on sort)
    
#   pour sortir de la boucle si l'utilisateur clique sur la croix de la waitbar
    abort=multiWaitbar(strWbar,tFin/tpsAcquisition);
    if abort       
        break # sort de la boucle while
    end

end
multiWaitbar('closeall'); # ferme la waitbar
print(['temps total d''acquisition : ' num2str(tFin) 's']); # temps total

resultat=resultat'; # impose d'etre un vecteur colonne

# enleve les 0 de la fin du vecteur dus à une initialisation approximative
resultat(resultat==0)=[];


# ------------------------------------------------------------
# deconnection
# ------------------------------------------------------------
deconnection;


# ------------------------------------------------------------
# Post-traitements des données qui arrivent de la carte
# ------------------------------------------------------------

# ------------------------------------------------------------
# tri dans les mesures

# trouver quelle valeur correspond à quoi dans 'resultat' en fonction de ordreLecture
# ex: Ang1=resultat(1:2:end);

# ---- Angle 1
indAng1=find([ordreLecture,ordreLecture]==bAng1);
indAng1=indAng1(1):(indAng1(2)-indAng1(1)):length(resultat);
Ang1=resultat(indAng1);

# ---- Angle 2
indAng2=find([ordreLecture,ordreLecture]==bAng2);
indAng2=indAng2(1):(indAng2(2)-indAng2(1)):length(resultat);
Ang2=resultat(indAng2);

# ---- Acceleration
indAcc=setdiff([1:length(resultat)],sort([indAng1,indAng2]));
Acc=resultat(indAcc);


# ----------------------------
# conversion binaire, valeur reelle

# fonction de transformation bit -> degre
transAng=@(bit) bit*360/((2^16)*0.8); # le 0.8 est du à la bande de 10# -> 90# du capteur

Ang1=transAng(Ang1);
Ang2=transAng(Ang2);

# fonction de transformation bit -> g
transAcc=@(bit) bit*500*9.5/(2^16*typeAccelero);        # fonction de la spécificité de l'accéléro utilisé : typeAccelero=4.4 ou 8
Acc=transAcc(Acc);


# ----------------------------
# increment de temps

# ---- Angle 1
tAng1=double(tempo)*indAng1; # en us (utilisation de 'double' car temps est un uint)
# ---- Angle 2
tAng2=double(tempo)*indAng2; # en us
# ---- Acceleration
tAcc=double(tempo)*indAcc; # en us



# ------------------------------------------------------------
# Trace toutes les mesures d'un coup
# ------------------------------------------------------------
# figure
# hold on
# plot(tAng1*1e-6,Ang1,'.','color',rgb('MediumBlue'));
# plot(tAng2*1e-6,Ang2,'.','color',rgb('LightCoral'));
# ylabel('ngle en degrés')
# xlabel('Temps en secondes')
# legend('Angle capteur 1','angle capteur 2');
# grid minor
# figure
# plot(tAcc*1e-6,Acc,'.','color',rgb('ForestGreen'))
# ylabel('Accélération en g')
# xlabel('Temps en secondes')
# grid minor


# ------------------------------------------------------------
# Enregistre dans un fichier csv et donne output 'mesures' de la fonction
# ------------------------------------------------------------
switch typeTest
    
    case 1 # calibrage balle
        
        # ---- output
        mesures=[tAng2',Ang2]; # capteur angle 2
        
        # ---- plot
        figure
        plot(tAng2*1e-6,Ang2,'.','color',rgb('LightCoral'));
        ylabel('Angle en degrés')
        xlabel('Temps en secondes')
        grid minor
        title('Calibrage balle')
        
        # ---- save csv  
        fid=fopen([current.data '/calibBalle-' date '.csv'],'a');
        fprintf(fid,'*************************************\n');
        fprintf(fid,'#s\n',datestr(now));
        fprintf(fid,'Angle (deg) pour un increment de temps de #uus\n',tAng2(2)-tAng2(1));
        fprintf(fid,'(non étalonné)\n');
        fprintf(fid,'#7g\n',[Ang2']);
        fclose(fid);
        
    case 2 #  calibrage raquette        
        
        # ---- output
        mesures=[tAng1',Ang1]; # capteur angle 1
        
        # ---- plot
        figure
        plot(tAng1*1e-6,Ang1,'.','color',rgb('MediumBlue'));
        ylabel('Angle en degrés')
        xlabel('Temps en secondes')
        grid minor
        title('Calibrage raquette')
        
        # ---- save        
        fid=fopen([current.data '/calibRaquette-' nomRaquette '-' date '.csv'],'a');
        fprintf(fid,'*************************************\n');
        fprintf(fid,'#s\n',datestr(now));
        fprintf(fid,'Angle (deg) pour un increment de temps de #uus\n',tAng1(2)-tAng1(1));
        fprintf(fid,'(non étalonné)\n');
        fprintf(fid,'#7g\n',[Ang1']);
        fclose(fid);
        
    case 3 # impact 
        
        # ---- output
        mesures{1}=[tAng1',Ang1]; # capteur angle 1 => raquette
        mesures{2}=[tAng2',Ang2]; # capteur angle 2 => balle
        mesures{3}=[tAcc',Acc];   # accéléro
        
        # ---- plot
        figure
        hold on
        plot(tAng1*1e-6,Ang1,'.','color',rgb('MediumBlue'));
        plot(tAng2*1e-6,Ang2,'.','color',rgb('LightCoral'));
        ylabel('Angle en degrés')
        xlabel('Temps en secondes')
        legend('Angle Raquette','Angle Balle');
        grid minor
        figure
        plot(tAcc*1e-6,Acc,'.','color',rgb('ForestGreen'))
        ylabel('Acceleration en g')
        xlabel('Temps en secondes')
        grid minor

        # ---- save        
        fid=fopen([current.data '/impact-' nomRaquette '-' date '.csv'],'a');
        fprintf(fid,'*************************************\n');
        fprintf(fid,'#s\n',datestr(now));
        fprintf(fid,'Angles (deg) pour un increment de temps de #uus\n',tAng1(2)-tAng1(1));         # ! l'increment sur Ang1 et Ang2 doit etre le meme
        fprintf(fid,'(non étalonné)\n');
        fprintf(fid,'#8s #8s\n','balle','raquette');
        fprintf(fid,'#7g #7g\n',[Ang2';Ang1']);
        fprintf(fid,'-------------------------------------\n');
        fprintf(fid,'Accélération (g) pour un increment de temps de #uus\n',tAcc(2)-tAcc(1));
        fprintf(fid,'(non étalonné)\n');
        fprintf(fid,'#8g\n',Acc');
        fclose(fid);
        
    otherwise # rien
        
end



# fonction annexe
function deconnection

# ------------------------------------------------------------
# deconnection
# ------------------------------------------------------------
    
if exist('c','var');
    try
        fclose(c); 
        clear('c'); 
        print(['déconnecté du port ' nomPortCOM])
    catch # rien
    end
end

end

# ------------------------------------------------------------
# end final de la fonction principale
end

