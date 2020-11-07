#-------------------------------------------------------------------------------
#  C00_REDLAW - Computes the interstellar extinction function 
#   A(lambda)/A(V) of Calzetti+2000

def calzetti(wave, rv=4.05):

  #extl = wave; extl[] = 0
    extl = pd.DataFrame([], index = range(len(wave)), columns = [0])
    npts = len(wave)
  # Convert input wavelength to microns
    wave = wave / 1e4
  
    MIN_WAVE = 0.12 # Minimum wave where function is defined.
    for pix in range(0,npts):
        if ((wave[0][pix] >= 0.12) & (wave[0][pix] < 0.63)):
            y = 2.659 * (-2.156 + (1.509 / wave[0][pix]) - (0.198 / wave[0][pix]**2) + (0.011 / wave[0][pix]**3) ) + rv
        
        if ((wave[0][pix] >= 0.63) & (wave[0][pix] < 2.2)):
            y = 2.659 * (-1.857 + 1.040 / wave[0][pix]) + rv
        
        # y = k(lambda) = A(lambda)/E(B-V); R_V = A(V)/E(B-V)
        #  so, A(lambda) / A(V) = y / R_V
        extl[0][pix] = y / rv
        # for(pix in 1:npts)
        
    return(extl)

#############################
# CORREÇÃO POR A(LAMBDA)
# 1) CALCULA Av A PARTIR DE CALZETTI ET AL. (2000)
# 2) CALCULA Alambda/Av utilizando a função acima
# 3) OBTÉM-SE Alambda A PARTIR DE MULTIPLICAÇÃO DE MATRIZES
# 4) APLICA-SE A CORREÇÃO SOBRE OS FLUXOS DAS LINHAS DE CADA ESPECTRO UTILIZANDO F_corr = F_obs*10**(0.4*Alambda)
#############################

import numpy as np
import pandas as pd
from matplotlib import pyplot as plt
pd.set_option('max_columns', None)
pd.set_option('max_rows', None)

flag_sample_sdss = 'sample'
flag_others = 'off'

lam_oiii1 = 5008.239                                      
lam_oiii2 = 4960.295
lam_hbeta = 4862.721
lam_halpha = 6564.61
lam_nii1 = 6585.27
lam_nii2 = 6549.86
lam_oii1 = 3727.092
lam_oii2 = 3729.875
lam_sii1 = 6718.29
lam_sii2 = 6732.68



if flag_sample_sdss == 'control':
    sample = pd.read_csv('/home/vitorbootz/research/aux_files/control_sample_flux_lines.csv')
    sample = sample[(sample.lcgID != 2361) & (sample.lcgID != 2023) & (sample.ra_lcg != 183.40380) & (sample.ra_lcg != 164.07920)]
    sample.columns = ['lcgID', 'ra', 'dec', 'specObjID', 'OII1', 'OII2', 'Hb', 'OIII_4959', 'OIII_5007', 'NII_6583', 'Ha', 'NII_6548', 'SII_6718', 'SII_6732', 'lgm_tot_p50', 'sfr_tot_p50', 'specsfr_tot_p50']
    sample = sample.sort_values('lcgID')
    wave = pd.DataFrame([lam_oii1, lam_oii2, lam_hbeta, lam_oiii2, lam_oiii1, lam_nii1, 
                         lam_halpha, lam_nii2, lam_sii1, lam_sii2])
if flag_sample_sdss == 'sample':
    sample = pd.read_csv('/home/vitorbootz/research/flux_measurements/sample_flux_lines.csv')  
    if flag_others == 'on':
        sample = sample[(sample.lcgID != 2361) & (sample.lcgID != 2023)  & (sample.flag_lcg == 1) &(sample.flag_sdss == 0) & (sample.ra != 183.40380) & (sample.ra != 164.07920)]
    else:
        sample = sample[(sample.lcgID != 2361) & (sample.lcgID != 2023)  & (sample.flag_lcg == 1) & (sample.ra != 183.40380) & (sample.ra != 164.07920)]
    sample = sample.sort_values('lcgID')
    wave = pd.DataFrame([lam_oiii1, lam_oiii2, lam_hbeta, lam_halpha, lam_nii1,
                    lam_nii2, lam_oii1, lam_oii2, lam_sii1, lam_sii2])
if flag_sample_sdss == 'sdss':
    sample = pd.read_csv('/home/vitorbootz/research/flux_measurements/sdss_flux_lines.csv')
    #sample = sample.drop(['specObjID'], axis=1)
    sample.columns = ['lcgID', 'ra', 'dec', 'specObjID', 'OII1', 'OII2', 'Hb', 'OIII_4959', 'OIII_5007', 'NII_6583', 'Ha', 'NII_6548', 'SII_6718', 'SII_6732', 'lgm_tot_p50', 'sfr_tot_p50', 'specsfr_tot_p50']
    sample = sample.sort_values('lcgID')
    wave = pd.DataFrame([lam_oii1, lam_oii2, lam_hbeta, lam_oiii2, lam_oiii1, lam_nii1, 
                         lam_halpha, lam_nii2, lam_sii1, lam_sii2])

sample.index = range(len(sample))

Av = pd.DataFrame([7.98*np.log10(sample.Ha/(2.87*sample.Hb))]).T # Calzetti (2000)

Al_Av = calzetti(wave, 4.05)

matrix = Al_Av @ Av.T # multiplicação matricial N x 1 vs 1 x N

matrix.index = sample.iloc[:, 4:14].columns
matrix.columns = sample.lcgID
matrix = matrix.T

backup = sample.copy()
sample.index = sample.lcgID
sample = sample.drop(['lcgID'], axis=1)
if flag_sample_sdss == 'control':
    sample_corr = pd.DataFrame([sample.index,  sample.ra, sample.dec]).T
if flag_sample_sdss == 'sample':
    sample_corr = pd.DataFrame([sample.index, sample.extension, sample.ra, sample.dec]).T
if flag_sample_sdss == 'sdss':
    sample_corr = pd.DataFrame([sample.index,  sample.ra, sample.dec]).T

m = sample.iloc[:, 3:13].copy() * 10**(0.4*matrix) # o resultado de cada elemento é a multiplicação do elemento ij de cada matriz
m.index = range(len(m))
sample.index = range(len(sample))
sample_corr = pd.concat([backup.lcgID, sample.iloc[:, 0:3] ,m, sample.iloc[:,13:]], ignore_index=True, axis=1)
sample_corr.columns = backup.columns

if flag_sample_sdss == 'control':
    sample_corr.to_csv('/home/vitorbootz/research/aux_files/lcgs_fluxos_corrigidos_control.csv', index=False)
if flag_sample_sdss == 'sample':
    if flag_others == 'on':
        sample_corr.to_csv('/home/vitorbootz/research/aux_files/lcgs_fluxos_corrigidos_sample_lcgs_gemini.csv', index=False)
    else:
        sample_corr.to_csv('/home/vitorbootz/research/aux_files/lcgs_fluxos_corrigidos_sample.csv', index=False)
if flag_sample_sdss == 'sdss':
    sample_corr.to_csv('/home/vitorbootz/research/aux_files/lcgs_fluxos_corrigidos_sdss.csv', index=False)
    
sample_corr = sample_corr.drop_duplicates('lcgID')
sample_corr.index = range(len(sample_corr))
sample = pd.concat([backup.lcgID, sample], axis=1)
sample_corr
