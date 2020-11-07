##############################################################################################
# PURPOSE
#   Measures the oxygen abundances of galaxies using Pilyugin & Grebel (2016)
# 
# CREATED BY:
#   Vitor Eduardo Buss Bootz
#
# ADAPTED BY:
#
# CALLING SEQUENCE
#   python abundancia_oxigenio.py     --> In terminal
#
# INPUT PARAMETERS
#   lcgs_fluxos_corrigidos_sample.csv           --> Corrected flux sample filename
#   lcgs_fluxos_corrigidos_sdss.csv             --> Corrected flux sdss filename
#   lcgs_fluxos_corrigidos_control.csv          --> Corrected flux control filename
#
# OUTPUT
#   abundancias_sample.csv
#   abundancias_control.csv
#   abundancias_sdss.csv
#
# REQUIRED SCRIPTS
#
# COMMENTS
#   The 'flag_sample_sdss' flag allows the measurement of the Oxygen Abundances for different samples.
#   Just choose between 'sample', 'sdss' and 'control'
#
##############################################################################################

#############################
# Abundância Oxigênio
#############################

import numpy as np
import pandas as pd
from matplotlib import pyplot as plt

flag_sample_sdss = 'control'
grupos = pd.read_csv('/home/vitorbootz/research/aux_files/lcgs_fluxos_corrigidos_sample.csv')

if flag_sample_sdss == 'sample':
    sample = pd.read_csv('/home/vitorbootz/research/aux_files/lcgs_fluxos_corrigidos_sample.csv')
    sample = sample[(sample.lcgID != 2361) & (sample.flag_lcg == 1)]
if flag_sample_sdss == 'sdss':
    sample = pd.read_csv('/home/vitorbootz/research/aux_files/lcgs_fluxos_corrigidos_sdss.csv')
    sample = sample[(sample.lcgID != 2361) & (sample.lcgID != 2023) & (sample.ra != 131.36510) & (sample.ra != 164.07920)]
if flag_sample_sdss == 'control':
    sample = pd.read_csv('/home/vitorbootz/research/aux_files/lcgs_fluxos_corrigidos_control.csv')

sample.index = range(len(sample))
    
if flag_sample_sdss == 'sample':
    OIII1 = sample.OIII_5008
if (flag_sample_sdss == 'control') | (flag_sample_sdss == 'sdss'):
    OIII1 = sample.OIII_5007
OIII2 = sample.OIII_4959
Hb = sample.Hb
Ha = sample.Ha
NII1 = sample.NII_6583
NII2 = sample.NII_6548
SII1 = sample.SII_6718
SII2 = sample.SII_6732

a1 = 8.424
a2 = 0.030
a3 = 0.751
a4 = -0.349
a5 = 0.182
a6 = 0.508
a7 = 8.072
a8 = 0.789
a9 = 0.726
a10 = 1.069
a11 = -0.170
a12 = 0.022

N2 = (NII1+NII2)/Hb
S2 = (SII1+SII2)/Hb
R3 = (OIII1+OIII2)/Hb

if flag_sample_sdss == 'control':
    ratios = pd.DataFrame([sample.lcgID, N2, S2, R3, sample.lgm_tot_p50]).T
    ratios.columns = ['lcgID', 'N2', 'S2', 'R3', 'lgm_tot_p50']
if flag_sample_sdss == 'sample':
    ratios = pd.DataFrame([sample.lcgID, N2, S2, R3, sample.Mcor_log_Mo]).T
    ratios.columns = ['lcgID', 'N2', 'S2', 'R3', 'Mcor_log_Mo']
if flag_sample_sdss == 'sdss':
    ratios = pd.DataFrame([sample.lcgID, N2, S2, R3, sample.lgm_tot_p50]).T
    ratios.columns = ['lcgID', 'N2', 'S2', 'R3', 'lgm_tot_p50']
    
ratios['OH'] = -0.99

for i in range(len(ratios)):
    if ratios.N2[i] >= -0.6:
        ratios.OH[i] = a1 + a2*np.log10(ratios.R3[i]/ratios.S2[i]) + a3*np.log10(ratios.N2[i]) + (a4 + a5*np.log10(ratios.R3[i]/ratios.S2[i]) + a6*np.log10(ratios.N2[i]))*np.log10(ratios.S2[i])
    if ratios.N2[i] < -0.6:
        ratios.OH[i] = a7 + a8*np.log10(ratios.R3[i]/ratios.S2[i]) + a9*np.log10(ratios.N2[i]) + (a10 + a11*np.log10(ratios.R3[i]/ratios.S2[i]) + a12*np.log10(ratios.N2[i]))*np.log10(ratios.S2[i])
        
comparacao = pd.DataFrame([])
#for i in range(len(yuri)):
##    for j in range(len(sample)):
#        if (yuri.lcgID[i] == sample.lcgID[j]):
#            comparacao = comparacao.append(pd.DataFrame(yuri.loc[i]).T)

comparacao = comparacao.drop_duplicates('lcgID')
ratios = ratios.drop_duplicates('lcgID')
comparacao.index = range(len(comparacao))
ratios = ratios.sort_values('lcgID')
ratios.index = range(len(ratios))

sample = sample.drop_duplicates('lcgID')
sample.index = range(len(sample))
output = pd.concat([sample, ratios.OH], axis=1) 
if flag_sample_sdss == 'sample':
    output.to_csv('/home/vitorbootz/research/aux_files/abundancias_sample.csv', index=False)
if flag_sample_sdss == 'control':
    output.to_csv('/home/vitorbootz/research/aux_files/abundancias_control.csv', index=False)
if flag_sample_sdss == 'sdss':    
    output.to_csv('/home/vitorbootz/research/aux_files/abundancias_sdss.csv', index=False)
output
