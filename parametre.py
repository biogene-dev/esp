# -*- coding: utf-8 -*-
"""
Created on Tue Mar 19 11:43:39 2019

@author: ash
"""
import math
###########################################################################
#                                                                         #
#       Calcul du centre de gravité de l'ensemble du bras du Thorbot      #
#                                                                         #
###########################################################################


## Définition des masses en [kg]

##################### Pendule raquette ####################################
# masses des grosses pièces
m_corne = 2.074; # masse de la corne sans les 8 vis d'attaches à l'axe de 
                 # rotation, sans la vis de l'aimant et sans les 3 vis entre le link et la corne
m_aimant_avec_vis = 0.378; # masse de l'aimant avec sa vis de fixation
m_link_padel = 0.580; # masse du link pour le padel
m_link_tennis = 0.180; # masse du link pour le tennis
m_attache = 0.828; # masse de l'attache avec la pièce interne, mais sans les 3 vis de fixation de l'attache au link
# m_piece_interne_attache = 0.120;
# m_attache_sans_piece_interne = 0.708;
m_axe_raquette = 4.330;
m_ressort = 0.514; #avec support
m_link_bras_corne = 0.350; # x2



# masses des petites pieces ( les masses sont données pour une vis avec son
# écrou et sa rondelle si il y en a)
m_vis_fixation_corne = 0.016; # x8. Vis de fixation de la corne au bras de rotation
m_vis_corne_link = 0.016; # x3. vis de fixation du link a la corne
m_vis_link_attache = 0.014; #x3. 
m_vis_bras_corne = 0.016; # x8


# masses totales 
m_tot_padel = m_axe_raquette + 2 * m_link_bras_corne + m_ressort + 8 * m_vis_bras_corne 
+ m_corne + m_aimant_avec_vis + 3 * m_vis_corne_link + m_link_padel 
+ 3 * m_vis_link_attache + m_attache;
          
m_tot_tennis = m_axe_raquette + 2 * m_link_bras_corne + m_ressort + 8 * m_vis_bras_corne + m_corne + m_aimant_avec_vis + 3 * m_vis_corne_link + m_link_tennis + 3 * m_vis_link_attache + m_attache;

              
              
######################## Pendule balle ####################################
m_balle = 0.290;
m_tige = 1.212;
m_axe_balle = 4.274;

# masses totales 
m_tot_pendule_balle = m_balle + m_tige + m_axe_balle;






## Définition des positions des centres de gravités par rapport à l'axe de rotation en [m]

##################### Pendule raquette ####################################
# Grosses pièces
x_corne = 0.070;
x_aimant = 0.245;

x_link_tennis = 0.292; #milieu du link
x_link_padel = 0.369; # milieu du link

x_attache_tennis = 0.349; # avec link tennis installé
x_attache_padel = 0.530; # avec link padel installé

x_axe_raquette = 0;
x_ressort = 0; 
x_link_bras_corne = 0; # x2


#Petites pièces
# On néglige les vis de fixation de la corne a l'axe de rotation car les
# moments d'inerties sont compensés de part et d'autre de l'axe
x_vis_corne_link = 0.265; # vis de fixation du link a la corne
x_vis_link_attache_padel = 0.473; # vis de fixation de l'attache au link padel
x_vis_link_attache_tennis = 0.322; # vis de fixation de l'attache au link tennis
x_vis_bras_corne = 0.023;



######################## Pendule balle ####################################
x_balle = 0.385;
x_tige = 0.777;
x_axe_balle = 0; #centre de rotation de l'axe est sur le centre de gravité

## Definition des rayons de l'axe  et du ressort 
r_axe = 0.0325/2;
r_ressort_interne = 0.0325/2;
r_ressort_externe = 0.063/2;
l_link_bras_corne = 0.055; # largeur
L_link_bras_corne = 0.065; # longueur




## Définition des moments d'inerties

##################### Pendule raquette ####################################
# Grosses pièces
I_corne = m_corne * math.pow(x_corne,2);
I_aimant = m_aimant_avec_vis * math.pow(x_aimant,2);
I_link_padel = m_link_padel * math.pow(x_link_padel,2);
I_link_tennis = m_link_tennis * math.pow(x_link_tennis,2);
I_attache_tennis = m_attache * math.pow(x_attache_tennis,2);
I_attache_padel = m_attache * math.pow(x_attache_padel,2);
I_axe_raquette = 0.5 * m_axe_raquette * math.pow(r_axe,2);
I_ressort = 0.5 * m_ressort * (math.pow(r_ressort_externe,2) - math.pow(r_ressort_interne,2));
I_link_axe_corne = m_link_bras_corne *((math.pow(L_link_bras_corne,2)) + math.pow(l_link_bras_corne,2)/12 - math.pow(r_axe,2)/2);

# Petites pièces
I_vis_corne_link = m_vis_corne_link * math.pow(x_vis_corne_link,2);
I_vis_link_attache_padel = m_vis_link_attache * math.pow(x_vis_link_attache_padel,2);
I_vis_link_attache_tennis = m_vis_link_attache * math.pow(x_vis_link_attache_tennis,2);
I_vis_bras_corne = m_vis_bras_corne * math.pow(x_vis_bras_corne,2); # 4 vis au-dessus et 4 vis en-dessous de l'axe de rotation. Donc les moments d'inertie s'annulent


# Moments d'inerties totaux
I_tot_padel = I_corne + I_aimant + I_link_padel + I_attache_padel
+ 3 * I_vis_corne_link + 3 * I_vis_link_attache_padel + I_axe_raquette 
+ I_ressort + 2 * I_link_axe_corne + 8 * I_vis_bras_corne;
              
I_tot_tennis = I_corne + I_aimant + I_link_tennis + I_attache_tennis
+ 3 * I_vis_corne_link + 3 * I_vis_link_attache_tennis + I_axe_raquette
+ I_ressort + 2 * I_link_axe_corne + 8 * I_vis_bras_corne;

               
               
######################## Pendule balle ####################################
I_balle = m_balle * math.pow(x_balle,2);
I_tige = m_tige * math.pow(x_tige,2);
I_axe_balle = 0.5 * m_axe_balle * math.pow(r_axe,2);


# Moments d'inerties totaux
I_tot_pendule_balle = I_balle + I_tige + I_axe_balle ;




## Calcul des centre de gravités totaux

# Formule: m_tot * x_tot^2 = sum( m_i * x_i ^2 )


# Pendule Raquette

# Centre de gravité du bras pour padel
x_tot_padel = math.sqrt(I_tot_padel / m_tot_padel);

# centre de gravité du bras pour tennis
x_tot_tennis = math.sqrt(I_tot_tennis / m_tot_tennis);



#Pendule Balle

# centre de gravité du pendule balle
x_tot_pendule_balle = math.sqrt(I_tot_pendule_balle / m_tot_pendule_balle);

print('------------------------------------------------------------')
print('Centres de gravités:')
print(['Padel: ',str(x_tot_padel),' m ou ',str(100*x_tot_padel),' cm']);
print(['Tennis: ',str(x_tot_tennis),' m ou ',str(100*x_tot_tennis),' cm']);
print(['Pendule balle: ',str(x_tot_pendule_balle),' m ou ',str(100*x_tot_pendule_balle),' cm']);
print('------------------------------------------------------------')

# 
# sqrt((4.163*0.132^2 - 0.182*0.292^2 - 0.349^2*0.886 + 0.886*0.53^2 + 0.584*0.369^2)/4.565)






