##############################################################################################
# PURPOSE
#   Measures the Tidal Strength Estimator (Q) parameter and creates the plot 'Tipical_Q_values_GEMA' (Fig. 4.1 - TCC)
# 
# CREATED BY:
#   Vitor Eduardo Buss Bootz
#
# ADAPTED BY:
#
# CALLING SEQUENCE
#   python Tidal_Strength_Estimator.py     --> In terminal
#
# INPUT PARAMETERS
#   output_file           --> Output filename
#   list_galaxies         --> List of objects
#   starlight_file        --> Starlight results
#   GEMA-1.0.1.fits       --> GEMA sample data
#
# OUTPUT
#   Q_values.csv
#   Plot Tipical_Q_values_GEMA.pdf
#
# REQUIRED SCRIPTS
#
# COMMENTS
#
##############################################################################################

import pandas as pd
import numpy as np
pd.set_option('max_columns', None)
pd.set_option('max_rows', None)

list_galaxies = 'aux_files/galaxy_list.csv' 
starlight_file = 'results_starlight/STARLIGHT_MAIN_RESULTS.csv'
listgal = pd.read_csv(list_galaxies, sep=',')  # Lista de spectros e redshifts
selection = (listgal['onoff'] == 1) & ((listgal['extension'] < 10) | ((listgal['extension'] > 20) & (listgal['extension'] < 100)))
#Desejamos pegar galáxias <10 e maiores que 20 para não pegar as 1extension do sloan e nem as >100. PELAMOR CONSERTAR ESSE PROBLEMA DO INÍCIO, DEIXAR TUDO COMO FLAGS APENAS
listgal = listgal[selection]
listgal.index = range(len(listgal))
logQ_values = pd.DataFrame([])
FINAL_TABLE = pd.DataFrame([])

data = pd.read_csv(starlight_file, sep=',') #Resultados do starlight

for i in listgal.lcgID.drop_duplicates(): #Separando as ids dos grupos para iteração
    Q = pd.DataFrame([])
    R = pd.DataFrame([])
    RQ = pd.DataFrame([])
    
    ### LCG ###
    selection_lcg = (listgal.lcgID == i) & (listgal.flag_lcg == 1)
    ra_lcg = np.radians(float(listgal[selection_lcg].ra)) #Seleciona primeiro o conjunto de galáxias de um mesmo grupo e então seleciona a ra lcg e transforma o valor em float
    dec_lcg = np.radians(float(listgal[selection_lcg].dec))
    extension = float(listgal[selection_lcg].extension)
    Mcor_lcg_starlight = float(10**(data[(data.lcgID == i) & (data.extension == extension)].Mcor_log_Mo)) #Seleciona primeiro o conjunto de galáxias de um mesmo grupo e então seleciona aquela que possui a mesma ra da calculada acima (neste caso, lcg)
    DA_lcg = float(data[(data.lcgID == i) & (data.extension == extension)].DA_Mpc)
    petroR90_r_lcg = float(listgal[selection_lcg].petroR90_r_lcg) * DA_lcg * np.pi * 1000/(180*3600) #Em kpc, pois arcsec.Mpc.(pi_rad/180_deg).(1_deg/3600arcsec).(1000kpc/1Mpc)
    
    #print(i,int(extension),round(float(data[(data.lcgID == i) & (data.extension == extension)].Mcor_log_Mo),2), float(listgal[selection_lcg].logMstar_lcg_sdss))
    
    ### Galaxy i ###
    other_gal = listgal[listgal.flag_lcg == 0]
    for j in range(len(other_gal[other_gal.lcgID == i])):
        group_i = other_gal[other_gal.lcgID == i]
        group_i.index = range(len(group_i))

        ra_j = np.radians(group_i.ra[j])
        dec_j = np.radians(group_i.dec[j])
        extension = group_i.extension[j]
        Mcor_j_starlight = 10**float(data[(data.lcgID == i) & (data.extension == extension)].Mcor_log_Mo)
        DA_i = float(data[(data.lcgID == i) & (data.extension == extension)].DA_Mpc)

        ### DISTANCE CALCULATOR ###
        theta = np.arccos(np.sin(dec_lcg)*np.sin(dec_j)+np.cos(dec_lcg)*np.cos(dec_j)*np.cos(ra_j - ra_lcg)) #In rad
        Rlcg_i = theta * 1000 * (DA_lcg+DA_i) / 2 #DA's mean, result in Kpc
        #Rlcg_i_theta = np.sqrt((DA_lcg*1000)**2 + (DA_i*1000)**2 - 2 * DA_lcg*1000 * DA_i*1000 * np.cos(theta)) # Lei dos cossenos a²=b²+c²-2.b.c.cos(theta)
        ###
        
        Qlcg_i = float((Mcor_j_starlight / Mcor_lcg_starlight) * ( 2*petroR90_r_lcg / Rlcg_i )**3)
        
        Q = Q.append([Qlcg_i])
        R = R.append([Rlcg_i])
        
            
        
        
        #print(i,extension, ra_j, dec_j, 2*petroR90_r_lcg, Rlcg_i)
        #print(i,extension,Qlcg_i)
    
    logQ = float(np.log10(sum(Q.values)))
    
    RQ = pd.concat([R,Q], axis=1)
    RQ.columns = ['Rlcg_i', 'Q']
    RQ = RQ.sort_values(['Rlcg_i'])
    RQ.index = range(len(RQ))
    RQ.Q = np.log10(RQ.Q)
    
    final_row = pd.concat([RQ[:1],pd.DataFrame([RQ.loc[len(RQ)-1].Rlcg_i]),pd.DataFrame([logQ]),pd.DataFrame([2*petroR90_r_lcg])], axis=1)
    final_row.columns = ['R_1st','Q_1st','R_last','Q_group', 'D_LCG']
    
    FINAL_TABLE = FINAL_TABLE.append(final_row)
    
FINAL_TABLE.index = listgal.lcgID.drop_duplicates()
FINAL_TABLE.to_csv('/home/vitorbootz/research/aux_files/Q_values.csv')
FINAL_TABLE

#################################################################################
# PLOT
#################################################################################

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
import matplotlib.lines as mlines
from astropy.io import fits as pf
from astropy.table import Table

hdul = pf.open('aux_files/GEMA-1.0.1.fits', memmap=True)
evt_data = Table(hdul[13].data) # Groups
our_env = evt_data[(evt_data['GroupSize'] > 2) & (evt_data['GroupSize'] < 6)] #Only grupos with 3-5 galaxies
our_env = our_env[our_env['Q_group'] < 9] #Removing outlaiers

our_env['dnn'] = our_env['dnn']*1000 #Mpc to kpc
our_env['dkn'] = our_env['dkn']*1000 #Mpc to kpc

fig, ax = plt.subplots()
cinza1 = ax.plot('dnn', 'Q_nn', 'o', color='#C0C0C0', data=our_env, markersize=2, label='Galáxia mais próxima - GEMA')
cinza2 = ax.plot('dkn', 'Q_group', '^', color='#808080', alpha=0.8,data=our_env, markersize=3, label='Galáxia mais distante - GEMA')
ax.set_xscale('log')

plt.title('Valores típicos de Q para grupos de 3 a 5 galáxias', fontweight='bold')
plt.xlabel('Distância projetada entre a LCG e a galáxia [kpc]', fontweight=False, size=12)
plt.ylabel('log(Qi) entre a LCG e a galáxia mais próxima', fontweight=False, size=9)

outliers = FINAL_TABLE[(FINAL_TABLE.index == 2361)]
useful = FINAL_TABLE[(FINAL_TABLE.index != 3090) & (FINAL_TABLE.index != 2361)]

azul1 = ax.plot('R_1st', 'Q_1st', 'o', color='#66B2FF', data=useful, markersize=5,label='Galáxia mais próxima \nValores: diâmetro da LCG [kpc]')
azul2 = ax.plot('R_last', 'Q_group', '^', color='#0000FF', data=useful, markersize=7, label='Galáxia mais distante')

outlier = ax.plot('R_1st', 'Q_1st', 'o', color='red', data=outliers, markersize=3, label='ID 2361: Par de galáxias')
#ax.plot('R_last', 'Q_group', '^', color='red', data=outliers, markersize=3, label='_nolegend_')

ax.legend(loc='lower left', shadow=True, fontsize='x-small')

text = useful.copy()
useful_text = text[round(text.D_LCG,2) != 6.39] #Só pra não deixar dois números sobrepostos no plot. No caso i=2247
for i in useful_text.index:
    ax.annotate(str(round(useful_text.D_LCG[i],1)), (useful_text.R_1st[i], useful_text.Q_1st[i]+0.2), size=6)    
ax.annotate(str(round(text.D_LCG[2247],1)), (text.R_1st[2247]-10, text.Q_1st[2247]+0.2), size=6)
ax.annotate(str(round(outliers.D_LCG[2361],1)), (outliers.R_1st[2361], outliers.Q_1st[2361]+0.2), size=6)


secax = ax.secondary_yaxis('right')
secax.set_ylabel('Q do grupo', size=12)

plt.grid(alpha=0.5)

plt.savefig('/home/vitorbootz/research/TCC_images/GEMA/Tipical_Q_values_GEMA.pdf', format='pdf', figsize=(6,3), bbox_inches='tight')


