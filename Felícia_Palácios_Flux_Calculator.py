def interp(xin,xout,yin): #Tentar arrumar isso no futuro.. problemas na conversão entre R e Python
    Ninterpol = len(xout)
    yout = np.array([])
    
    for k in range(Ninterpol): 
        t = xin[xin < xout[k]]
        tind = len(t)
        
        if tind <= 0: tind = 1
        if tind >= len(xin): tind = len(xin) - 1
        t1 = xin[tind - 1]
        t2 = xin[tind]
        t3 = xin[tind + 1]
        tx = xout[k]

        A = (tx - t1) / (t3 - t1)
        B = (tx - t2) / (t3 - t2)
        C = (tx - t3) / (t2 - t1)
        D = (tx - t1) / (t3 - t2)
        E = (tx - t2) / (t3 - t1)

        G = (tx - t2) / (t3 - t2)
        H = (tx - t1) / (t2 - t1)
        
        if (t1 != t2) & (t2 != t3):
            yout = np.append(yout, yin[tind+1] * A * B - yin[tind] * D * C + yin[tind-1] * E * C)
                
        if (t1 == t2):
            yout = np.append(yout, yin[tind+1] - yin[tind]) * G + yin[tind] 
        
        if(t2 == t3):
            yout = np.append(yout, yin[tind] - yin[tind-1]) * H + yin[tind-1] 
        
    return(yout)

def interp1(xin,xout,yin):
    Ninterpol = len(xout)
    yout = np.array([])
    
    for k in range(Ninterpol): 
        t = xin[xin < xout[k]]
        tind = len(t) -1 #VERIFICAR SE ESTE -1 NÃO INTERFERE NOS RESULTADOS
        
        if tind <= 0: tind = 1
        if tind >= len(xin): tind = len(xin) - 1
        t1 = xin[tind - 1]
        t2 = xin[tind]
        t3 = xin[tind + 1]
        tx = xout[k]

        A = (tx - t1) / (t3 - t1)
        B = (tx - t2) / (t3 - t2)
        C = (tx - t3) / (t2 - t1)
        D = (tx - t1) / (t3 - t2)
        E = (tx - t2) / (t3 - t1)

        G = (tx - t2) / (t3 - t2)
        H = (tx - t1) / (t2 - t1)
        
        if (t1 != t2) & (t2 != t3):
            yout = np.append(yout, yin[tind+1] * A * B - yin[tind] * D * C + yin[tind-1] * E * C)
                
        if (t1 == t2):
            yout = np.append(yout, yin[tind+1] - yin[tind]) * G + yin[tind] 
        
        if(t2 == t3):
            yout = np.append(yout, yin[tind] - yin[tind-1]) * H + yin[tind-1] 
        
    return(yout)

def interp2(xin,xout,yin):
    Ninterpol = len(xout)
    yout = np.array([])
    
    for k in range(Ninterpol): 
        t = xin[xin < xout[k]]
        tind = len(t) -2 #VERIFICAR SE ESTE -1 NÃO INTERFERE NOS RESULTADOS
        
        if tind <= 0: tind = 1
        if tind >= len(xin): tind = len(xin) - 1
        t1 = xin[tind - 1]
        t2 = xin[tind]
        if k != Ninterpol:
            t3 = xin[tind + 1]
        tx = xout[k]

        A = (tx - t1) / (t3 - t1)
        B = (tx - t2) / (t3 - t2)
        C = (tx - t3) / (t2 - t1)
        D = (tx - t1) / (t3 - t2)
        E = (tx - t2) / (t3 - t1)

        G = (tx - t2) / (t3 - t2)
        H = (tx - t1) / (t2 - t1)
        
        if (t1 != t2) & (t2 != t3):
            yout = np.append(yout, yin[tind+1] * A * B - yin[tind] * D * C + yin[tind-1] * E * C)
                
        if (t1 == t2):
            yout = np.append(yout, yin[tind+1] - yin[tind]) * G + yin[tind] 
        
        if(t2 == t3):
            yout = np.append(yout, yin[tind] - yin[tind-1]) * H + yin[tind-1] 
        
    return(yout)

import numpy as np
import pandas as pd
import astropy.io.fits as fits
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit
import math
import copy
from matplotlib.backends.backend_pdf import PdfPages
import datetime

#declara as função gaussianas e os contínuos que serão ajustados aos dados
def gauss(x, a, b, c, d, cs, e, f, g, h, i, j, k, l, s1, s2):
    step_min_oiii = np.heaviside(lam_final - min_oiii, 0)
    step_max_oiii = np. heaviside(lam_final - max_oiii, 0)
    step_min_hbeta = np.heaviside(lam_final - min_hbeta, 0)
    step_max_hbeta = np.heaviside(lam_final - max_hbeta, 0)
    step_min_halpha = np.heaviside(lam_final - min_nii, 0)
    step_max_halpha = np.heaviside(lam_final - max_nii, 0)
    step_min_oii = np.heaviside(lam_final - min_oii, 0)
    step_max_oii = np.heaviside(lam_final - max_oii, 0)
    step_min_sii = np.heaviside(lam_final - min_sii, 0)
    step_max_sii = np.heaviside(lam_final - max_sii, 0)

    r1 = a*(step_min_oiii - step_max_oiii)
    r2 = b*(step_min_halpha - step_max_halpha)
    r3 = c*(step_min_oii - step_max_oii)
    r4 = d*(step_min_hbeta - step_max_hbeta) 
    r5 = cs*(step_min_sii - step_max_sii)
    
    w_oiii1 = (x - e)/f
    alpha_oiii1 = (1./math.sqrt(2.*math.pi))*np.exp((-w_oiii1**2.)/2.)*g

    cenoiii2 = lam_oiii2 + (lam_oiii2/lam_oiii1)*(e - lam_oiii1)
    sigoiii2 = (lam_oiii2/lam_oiii1)*f
    ampoiii2 = g/3.
    w_oiii2 = (x - cenoiii2)/sigoiii2
    alpha_oiii2 = (1./math.sqrt(2.*math.pi))*np.exp((-w_oiii2**2.)/2.)*ampoiii2

    cenhbeta = lam_hbeta + (lam_hbeta/lam_oiii1)*(e - lam_oiii1)
    sighbeta = (lam_hbeta/lam_oiii1)*f
    w_hbeta = (x - cenhbeta)/sighbeta
    alpha_hbeta = (1./math.sqrt(2.*math.pi))*np.exp((-w_hbeta**2.)/2.)*h

    cenhalpha = lam_halpha + (lam_halpha/lam_oiii1)*(e - lam_oiii1)
    sighalpha = (lam_halpha/lam_oiii1)*f
    w_halpha = (x - cenhalpha)/sighalpha
    alpha_halpha = (1./math.sqrt(2.*math.pi))*np.exp((-w_halpha**2.)/2.)*i

    cennii1 = lam_nii1 + (lam_nii1/lam_oiii1)*(e - lam_oiii1)
    signii1 = (lam_nii1/lam_oiii1)*f
    w_nii1 = (x - cennii1)/signii1
    alpha_nii1 = (1./math.sqrt(2.*math.pi))*np.exp((-w_nii1**2.)/2.)*j

    cennii2 = lam_nii2 + (lam_nii2/lam_oiii1)*(e - lam_oiii1)
    signii2 = (lam_nii2/lam_oiii1)*f
    ampnii2 = j/3.
    w_nii2 = (x - cennii2)/signii2
    alpha_nii2 = (1./math.sqrt(2.*math.pi))*np.exp((-w_nii2**2.)/2.)*ampnii2

    cenoii1 = lam_oii1 + (lam_oii1/lam_oiii1)*(e - lam_oiii1)
    sigoii1 = (lam_oii1/lam_oiii1)*f
    w_oii1 = (x - cenoii1)/sigoii1
    alpha_oii1 = (1./math.sqrt(2.*math.pi))*np.exp((-w_oii1**2.)/2.)*k

    cenoii2 = lam_oii2 + (lam_oii2/lam_oiii1)*(e - lam_oiii1)
    sigoii2 = (lam_oii2/lam_oiii1)*f
    w_oii2 = (x - cenoii2)/sigoii2
    alpha_oii2 = (1./math.sqrt(2.*math.pi))*np.exp((-w_oii2**2.)/2.)*l
    
    censii1 = lam_sii1 + (lam_sii1/lam_oiii1)*(e - lam_oiii1)
    sigsii1 = (lam_sii1/lam_oiii1)*f
    w_sii1 = (x - censii1)/sigsii1
    alpha_sii1 = (1./math.sqrt(2.*math.pi))*np.exp((-w_sii1**2.)/2.)*s1
    
    censii2 = lam_sii2 + (lam_sii2/lam_oiii1)*(e - lam_oiii1)
    sigsii2 = (lam_sii2/lam_oiii1)*f
    w_sii2 = (x - censii2)/sigsii2
    alpha_sii2 = (1./math.sqrt(2.*math.pi))*np.exp((-w_sii2**2.)/2.)*s2
    
    return (r1 + r2 + r3 + r4 + r5 + (alpha_oiii1/f) + (alpha_oiii2/sigoiii2) + (alpha_hbeta/sighbeta)
+ (alpha_halpha/sighalpha) + (alpha_nii1/signii1) + (alpha_nii2/signii2) + (alpha_oii1/sigoii1) 
+ (alpha_oii2/sigoii2) + (alpha_sii1/sigsii1) + (alpha_sii2/sigsii2))
    

#comprimentos de onda de laboratório
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


###
flag_correction = 'off'
###

### FILES & FOLDERS ###
starlight_input_dir = '/home/vitorbootz/research/files_input_starlight/'
starlight_output_dir = '/home/vitorbootz/research/files_output_starlight/'                 
if flag_correction == 'on':  spectra_z_file = '/home/vitorbootz/research/aux_files/to_correct.csv'            
if flag_correction == 'off': spectra_z_file = '/home/vitorbootz/research/aux_files/galaxy_list.csv'
output_dir = '/home/vitorbootz/research/flux_measurements/'

if flag_correction == 'on':  output_file = 'corrections.csv'
if flag_correction == 'off': output_file = 'sample_flux_lines.csv'

### OUTPUT FILE HEADER ###
header = pd.DataFrame(['lcgID','extension','ra','dec','OIII_5008','OIII_4959', 'Hb', 'Ha', 'NII_6583', 'NII_6548', 'OII1', 'OII2', 'SII_6718', 'SII_6732']).T
header.to_csv(output_dir+output_file, index=False, mode='w', header=False)

### SELECTION ###
data = pd.read_csv(spectra_z_file, sep=',')  # Lista de spectros e redshifts
selection = (data['onoff'] == 1) & (data['extension'] < 100) #Seleciona apenas as lcgs do sdss e as do gemini com abertura = 3 arcsec
data = data[selection]
data.index = range(len(data))

for i in range(len(data)):
    
    if data.flag_sdss[i] == 0: 
        file_vis = (starlight_input_dir + 'input_starlight_gemini_LCG' + str(data.lcgID[i]) + '_' + str(data.extension[i]) + '.csv')
        file_syn = (starlight_output_dir + 'output_starlight_gemini_LCG' + str(data.lcgID[i]) + '_' + str(data.extension[i]) + '.csv')

    if data.flag_sdss[i] == 1: 
        file_vis = (starlight_input_dir + 'input_starlight_sdss_LCG' + str(data.lcgID[i]) + '_' + str(data.extension[i]) + '.csv')
        if (data.extension[i] > 10) & (data.extension[i] < 20):
            file_syn = (starlight_output_dir + 'output_starlight_sdss_LCG' + str(data.lcgID[i]) + '_' + str(data.extension[i]-10) + '.csv')
        else:
            file_syn = (starlight_output_dir + 'output_starlight_sdss_LCG' + str(data.lcgID[i]) + '_' + str(data.extension[i]) + '.csv')
        
    
    hdu = pd.read_csv(file_vis, delim_whitespace=True, engine='python', header = None)
    starlight = pd.read_csv(file_syn, skiprows = 306, delim_whitespace=True, engine='python', header = None)

### O range escolhido para o starlight pega o menor e o maior lambda entre os n espectros de cada grupo ###
### Portanto, tem casos de lambdas iniciais e finais maiores e menores nos inputs em relação ao output e vice e versa ###
### Por isso aplico estas correções, para que tanto o input quando o output do starlight fiquem com o mesmo lambda range ###
    if data.flag_sdss[i] == 0:
        if (hdu[0][0] <= starlight[0][0]) & (hdu[0][len(hdu[0])-1] >= starlight[0][len(starlight[0])-1]):
            hdu = hdu[(hdu[0] >= starlight[0][0]) & (hdu[0] <= starlight[0][len(starlight[0])-1])]
            hdu.index = range(len(hdu))
        elif (hdu[0][0] > starlight[0][0]) & (hdu[0][len(hdu[0])-1] <= starlight[0][len(starlight[0])-1]):
            starlight = starlight[(starlight[0] >= hdu[0][0]) & (starlight[0] <= hdu[0][len(hdu[0])-1])]
            starlight.index = range(len(starlight))
        elif (hdu[0][0] >= starlight[0][0]) & (hdu[0][len(hdu[0])-1] >= starlight[0][len(starlight[0])-1]):
            starlight = starlight[starlight[0] >= hdu[0][0]]
            starlight.index = range(len(starlight))
            hdu = hdu[hdu[0] <= starlight[0][len(starlight[0])-1]]
            hdu.index = range(len(hdu))
        elif (hdu[0][0] <= starlight[0][0]) & (hdu[0][len(hdu[0])-1] <= starlight[0][len(starlight[0])-1]):
            hdu = hdu[hdu[0] >= starlight[0][0]]
            hdu.index = range(len(hdu))
            starlight = starlight[starlight[0] <= hdu[0][len(hdu[0])-1]]
            starlight.index = range(len(starlight))
        
        l_obs = hdu[0]
        f_obs = hdu[1]
        sig = hdu[2]
        f_syn = starlight[2]                          
        
#### Os espectros de input do sloan estão bem feios, então aplico uma interpolação aqui ###        
    if data.flag_sdss[i] == 1:

        lamb = hdu[0]
        flux = hdu[1]
        l_starlight = starlight[0]
        f_syn = starlight[2]
        sig   = hdu[2]
        
        try:
            f_obs = interp(lamb, l_starlight, flux)
            sig = interp(lamb, l_starlight, sig)
            l_obs = l_starlight
        except:
            try:
                f_obs = interp1(lamb, l_starlight, flux)
                sig = interp1(lamb, l_starlight, sig)
                l_obs = l_starlight
            except:
                try:
                    f_obs = interp2(lamb, l_starlight, flux)
                    sig = interp2(lamb, l_starlight, sig)
                    l_obs = l_starlight
                except:
                    print('Erro na interpolação: LCG' + str(data.lcgID[i]) + '_' + str(data.extension[i]))
                    continue
################
### Exceções: aqui foi necessário aplicar a correção feita nos arquivos do gemini E interpolar o espectro ###                    
    if (data.lcgID[i] == 2764) & (data.extension[i] == 1):
        hdu = pd.read_csv(file_vis, delim_whitespace=True, engine='python', header = None)
        starlight = pd.read_csv(file_syn, skiprows = 306, delim_whitespace=True, engine='python', header = None)
        starlight = starlight[starlight[0] >= hdu[0][0]]
        starlight.index = range(len(starlight))
        hdu = hdu[hdu[0] <= starlight[0][len(starlight[0])-1]]
        hdu.index = range(len(hdu))
        lamb = hdu[0]
        flux = hdu[1]
        l_starlight = starlight[0]
        f_syn = starlight[2]
        sig   = hdu[2]
        f_obs = interp1(lamb, l_starlight, flux)
        sig = interp1(lamb, l_starlight, sig)
        l_obs = l_starlight
  #  if (data.lcgID[i] == 2023) & (data.extension[i] == 1):
  #      hdu = pd.read_csv(file_vis, delim_whitespace=True, engine='python', header = None)
  #      starlight = pd.read_csv(file_syn, skiprows = 306, delim_whitespace=True, engine='python', header = None)
  #      
  #      interp2023 = pd.read_csv('/home/vitorbootz/research/aux_files/interp_2023.csv')
  #      interp2023.columns = ['n','flux_out']
  #      interp2023 = interp2023[interp2023.n>118].dropna()
  #      
  #      interpsig2023 = pd.read_csv('/home/vitorbootz/research/aux_files/interp_sig_2023.csv')
  #      interpsig2023.columns = ['n','sig_out']
  #      interpsig2023 = interpsig2023[interpsig2023.n>118].dropna()
  #      
  #      starlight = starlight[starlight[0] >= hdu[0][0]]
  #      starlight.index = range(len(starlight))
  #      hdu = hdu[hdu[0] <= starlight[0][len(starlight[0])-1]]
  #      hdu.index = range(len(hdu))
  #      lamb = hdu[0]
  #      flux = hdu[1]
  #      l_starlight = starlight[0]
  #      f_syn = starlight[2]
  #      sig   = hdu[2]
  #      f_obs = interp2023.flux_out
  #      sig = interpsig2023.sig_out
  #      l_obs = l_starlight
################

    cube = (f_obs - f_syn)
    #cube = (f_obs)
    cube.dropna(inplace=True)

    n = len(l_obs)
    l0 = l_obs[0]     #primeiro valor de lambda
    step = l_obs[7] - l_obs[6]          #intervalo delta lambda
    unit = 'Angstrom' #unidade

    #cria um vetor com comprimentos de onda
    if unit == 'Angstrom':
        lam = np.arange(n)*step + l0
    else:
        if unit == 'nm':
            lam = np.arange(n)*step*10 + l0*10

    #centro das linhas observadas
    cen_oiii1 = 5006.843
    cen_oiii2 = lam_oiii2 + (lam_oiii2/lam_oiii1)*(cen_oiii1 - lam_oiii1)
    cen_hbeta = lam_hbeta + (lam_hbeta/lam_oiii1)*(cen_oiii1 - lam_oiii1)
    cen_halpha = lam_halpha + (lam_halpha/lam_oiii1)*(cen_oiii1 - lam_oiii1)
    cen_nii1 = lam_nii1 + (lam_nii1/lam_oiii1)*(cen_oiii1 - lam_oiii1)
    cen_nii2 = lam_nii2 + (lam_nii2/lam_oiii1)*(cen_oiii1 - lam_oiii1)
    cen_oii1 = lam_oii1 + (lam_oii1/lam_oiii1)*(cen_oiii1 - lam_oiii1)
    cen_oii2 = lam_oii2 + (lam_oii2/lam_oiii1)*(cen_oiii1 - lam_oiii1)
    
    cen_sii1 = lam_sii1 + (lam_sii1/lam_oiii1)*(cen_oiii1 - lam_oiii1)
    cen_sii2 = lam_sii2 + (lam_sii2/lam_oiii1)*(cen_oiii1 - lam_oiii1)

    #cria janela do OIII
    min_oiii = cen_oiii2 - 20
    max_oiii = cen_oiii1 + 20

    imin_oiii = (np.abs(lam - min_oiii)).argmin()
    imax_oiii = (np.abs(lam - max_oiii)).argmin()

    lam1 = lam[imin_oiii : imax_oiii]
    cube1 = cube[imin_oiii : imax_oiii]
    sig1 = sig[imin_oiii : imax_oiii]

    #cria janela do H beta
    min_hbeta = cen_hbeta - 20
    max_hbeta = cen_hbeta + 20

    imin_hbeta = (np.abs(lam - min_hbeta)).argmin()
    imax_hbeta = (np.abs(lam - max_hbeta)).argmin()

    lam2 = lam[imin_hbeta : imax_hbeta]
    cube2 = cube[imin_hbeta : imax_hbeta]
    sig2 = sig[imin_hbeta : imax_hbeta]

    #cria janela do NII
    min_nii = cen_nii2 - 20
    max_nii = cen_nii1 + 20

    imin_nii = (np.abs(lam - min_nii)).argmin()
    imax_nii = (np.abs(lam - max_nii)).argmin()

    lam3 = lam[imin_nii : imax_nii]
    cube3 = cube[imin_nii : imax_nii]
    sig3 = sig[imin_nii : imax_nii]

    #cria janela do SII
    min_sii = cen_sii1 - 30
    max_sii = cen_sii2 + 30

    imin_sii = (np.abs(lam - min_sii)).argmin()
    imax_sii = (np.abs(lam - max_sii)).argmin()

    lam5 = lam[imin_sii : imax_sii]
    cube5 = cube[imin_sii : imax_sii]
    sig5 = sig[imin_sii : imax_sii]
    
    #cria janela do OII
    min_oii = cen_oii1 - 20
    max_oii = cen_oii2 + 20

    #abre o arquivo do espectro visual se as linhas de OII cairem fora do infravermelho
    if min_oii < 10000:
        #head_vis = hdu_vis[0].header
        #data_vis = hdu_vis[0].data
        #sig_vis = hdu_vis[1].data/1e-19

        #cube_vis = data_vis/1e-19

        cube_vis = cube
        sig_vis = sig

        n_vis = len(cube_vis)
        l0_vis = l_obs[0]
        step_vis = l_obs[7] - l_obs[6]
        unit_vis = 'Angstrom'

        if unit_vis == 'Angstrom':
            lam_vis = np.arange(n_vis)*step_vis + l0_vis
        else:
            if unit_vis == 'nm':
                lam_vis = np.arange(n_vis)*step_vis*10 + l0_vis*10    

        imin_oii = (np.abs(lam_vis - min_oii)).argmin()
        imax_oii = (np.abs(lam_vis - max_oii)).argmin()

        lam4 = lam_vis[imin_oii : imax_oii]
        cube4 = cube_vis[imin_oii : imax_oii]
        sig4 = sig_vis[imin_oii : imax_oii]

    else:

        imin_oii = (np.abs(lam - min_oii)).argmin()
        imax_oii = (np.abs(lam - max_oii)).argmin()

        lam4 = lam[imin_oii : imax_oii]
        cube4 = cube[imin_oii : imax_oii]
        sig4 = sig[imin_oii : imax_oii]

    #junta as janelas
    lam_fin1 = np.append(lam1,lam2)
    cube_fin1 = np.append(cube1,cube2)
    sig_fin1 = np.append(sig1,sig2)
    lam_fin2 = np.append(lam_fin1, lam3)
    cube_fin2 = np.append(cube_fin1, cube3)
    sig_fin2 = np.append(sig_fin1, sig3)
    lam_fin = np.append(lam_fin2, lam4)
    cube_fin = np.append(cube_fin2, cube4)
    sig = np.append(sig_fin2, sig4)
    
    lam_final = np.append(lam_fin, lam5)
    cube_final = np.append(cube_fin, cube5)
    sig_final = np.append(sig, sig5)

    #chute inicial dos parâmetros 
    if flag_correction == 'off': guess = np.array([10, 10, 10, 10, 10, cen_oiii1, 60*cen_oiii1/(2.99792*1e5), 100, 100, 100, 100, 100, 100, 100, 100])
    if flag_correction == 'on':  guess = np.array([10, 10, 10, 10, 10, cen_oiii1, 60*cen_oiii1/(2.99792*1e5), 100, 100, 100, 100, 100, 100, 100, 100])

    #ajuste das funções aos dados    
    try:
        if (data.lcgID[i] == 3866) & (data.extension[i] == 13):
            p, pv = curve_fit(gauss, lam_final, cube_final, guess, sigma=None, method='lm')
        else:
            p, pv = curve_fit(gauss, lam_final, cube_final, guess, sigma=sig_final, method='lm')
    except:
        print('Erro no curve_fit: LCG'+str(data.lcgID[i])+'_'+str(data.extension[i]))
        continue

    #parâmetros do ajuste
    a1 = p[0] #continuo 1
    a2 = p[1] #continuo 2
    a3 = p[2] #continuo 3
    a4 = p[3] #continuo 4
    a5 = p[4] #continuo 5
    center_oiii1 = p[5]
    sig_oiii1 = copy.deepcopy(p[6])
    amp_oiii1 = p[7]
    amp_hbeta = p[8]
    amp_halpha = p[9]
    amp_nii1 = p[10]
    amp_oii1 = p[11]
    amp_oii2 = p[12]
    amp_sii1 = p[13]
    amp_sii2 = p[14]

    #outros parâmetros
    center_oiii2 = lam_oiii2 + (lam_oiii2/lam_oiii1)*(center_oiii1 - lam_oiii1)
    center_hbeta = lam_hbeta + (lam_hbeta/lam_oiii1)*(center_oiii1 - lam_oiii1)
    center_halpha = lam_halpha + (lam_halpha/lam_oiii1)*(center_oiii1 - lam_oiii1)
    center_nii1 = lam_nii1 + (lam_nii1/lam_oiii1)*(center_oiii1 - lam_oiii1)
    center_nii2 = lam_nii2 + (lam_nii2/lam_oiii1)*(center_oiii1 - lam_oiii1)
    center_oii1 = lam_oii1 + (lam_oii1/lam_oiii1)*(center_oiii1 - lam_oiii1)
    center_oii2 = lam_oii2 + (lam_oii2/lam_oiii1)*(center_oiii1 - lam_oiii1)
    
    center_sii1 = lam_sii1 + (lam_sii1/lam_oiii1)*(center_oiii1 - lam_oiii1)
    center_sii2 = lam_sii2 + (lam_sii2/lam_oiii1)*(center_oiii1 - lam_oiii1)

    sigma_oiii2 = (lam_oiii2/lam_oiii1)*sig_oiii1
    sigma_hbeta = (lam_hbeta/lam_oiii1)*sig_oiii1
    sigma_halpha = (lam_halpha/lam_oiii1)*sig_oiii1
    sigma_nii1 = (lam_nii1/lam_oiii1)*sig_oiii1
    sigma_nii2 = (lam_nii2/lam_oiii1)*sig_oiii1
    sigma_oii1 = (lam_oii1/lam_oiii1)*sig_oiii1
    sigma_oii2 = (lam_oii2/lam_oiii1)*sig_oiii1
    
    sigma_sii1 = (lam_sii1/lam_oiii1)*sig_oiii1
    sigma_sii2 = (lam_sii2/lam_oiii1)*sig_oiii1

    amp_oiii2 = amp_oiii1/3.
    amp_nii2 = amp_nii1/3.


    vel = p[5]*(299792.458/p[4])
    
    #plot dos espectros e ajustes
    fig = plt.figure(figsize=(6,20))
    axis1 = fig.add_subplot(511)
    axis1.plot(lam1,cube1)
    x = lam1
    y = (a1 + amp_oiii1*(1./math.sqrt(2.*math.pi))*np.exp((-((x - center_oiii1)/sig_oiii1)**2.)/2.)/sig_oiii1
         + amp_oiii2*(1./math.sqrt(2.*math.pi))*np.exp((-((x - center_oiii2)/sigma_oiii2)**2.)/2.)/sigma_oiii2)
    axis1.plot(x,y, color = 'red')
    plt.title('LCG' + str(data.lcgID[i]) + '_' + str(data.extension[i])+': OIII Window')

    axis2 = fig.add_subplot(512)
    axis2.plot(lam2, cube2)
    x1 = lam2
    y1 = (a4 + amp_hbeta*(1./math.sqrt(2.*math.pi))*np.exp((-((x1 - center_hbeta)/sigma_hbeta)**2.)/2.)/sigma_hbeta )
    axis2.plot(x1,y1, color = 'red')
    plt.title('LCG' + str(data.lcgID[i]) + '_' + str(data.extension[i])+': H Beta Window')

    axis3 = fig.add_subplot(513)
    axis3.plot(lam3, cube3)
    x2 = lam3
    y2 = (a2 + amp_halpha*(1./math.sqrt(2.*math.pi))*np.exp((-((x2 - center_halpha)/sigma_halpha)**2.)/2.)/sigma_halpha
         + amp_nii1*(1./math.sqrt(2.*math.pi))*np.exp((-((x2 - center_nii1)/sigma_nii1)**2.)/2.)/sigma_nii1
         + amp_nii2*(1./math.sqrt(2.*math.pi))*np.exp((-((x2 - center_nii2)/sigma_nii2)**2.)/2.)/sigma_nii2)
    axis3.plot(x2,y2, color = 'red')
    plt.title('LCG' + str(data.lcgID[i]) + '_' + str(data.extension[i])+': NII Window')
    
    axis4 = fig.add_subplot(514)
    axis4.plot(lam4,cube4)
    x3 = lam4
    y3 = (a3 + amp_oii1*(1./math.sqrt(2.*math.pi))*np.exp((-((x3 - center_oii1)/sigma_oii1)**2.)/2.)/sigma_oii1
         + amp_oii2*(1./math.sqrt(2.*math.pi))*np.exp((-((x3 - center_oii2)/sigma_oii2)**2.)/2.)/sigma_oii2)
    axis4.plot(x3,y3, color = 'red')
    plt.title('LCG' + str(data.lcgID[i]) + '_' + str(data.extension[i])+': OII Window')
    
    axis5 = fig.add_subplot(515)
    axis5.plot(lam5,cube5)
    x4 = lam5
    y4 = (a5 + amp_sii1*(1./math.sqrt(2.*math.pi))*np.exp((-((x4 - center_sii1)/sigma_sii1)**2.)/2.)/sigma_sii1
         + amp_sii2*(1./math.sqrt(2.*math.pi))*np.exp((-((x4 - center_sii2)/sigma_sii2)**2.)/2.)/sigma_sii2)
    axis5.plot(x4,y4, color = 'red')
    
    plt.title('LCG' + str(data.lcgID[i]) + '_' + str(data.extension[i])+': SII Window')
   
    if flag_correction == 'off':
        if data.flag_sdss[i] == 1:
            fig.savefig('/home/vitorbootz/research/flux_measurements/log_images/' + 'sdss_LCG' + str(data.lcgID[i]) + '_' + str(data.extension[i]) + '.png')
        if data.flag_sdss[i] == 0:
            fig.savefig('/home/vitorbootz/research/flux_measurements/log_images/' + 'gemini_LCG' + str(data.lcgID[i]) + '_' + str(data.extension[i]) + '.png')

    if flag_correction == 'on':
        if data.flag_sdss[i] == 1:
            fig.savefig('/home/vitorbootz/research/flux_measurements/to_correct/' + 'sdss_LCG' + str(data.lcgID[i]) + '_' + str(data.extension[i]) + '.png')
        if data.flag_sdss[i] == 0:
            fig.savefig('/home/vitorbootz/research/flux_measurements/to_correct/' + 'gemini_LCG' + str(data.lcgID[i]) + '_' + str(data.extension[i]) + '.png')

    ###########################################################################
    ###########################################################################
    ###########################################################################
    
    ### OUTPUT DATA ###
    
    amp_oiii1 = round(amp_oiii1,5)
    amp_oiii2 = round(amp_oiii2,5)
    amp_hbeta = round(amp_hbeta,5)
    amp_halpha = round(amp_halpha,5)
    amp_nii1 = round(amp_nii1,5)
    amp_oii1 = round(amp_oii1,5)
    amp_oii2 = round(amp_oii2,5)
    amp_sii1 = round(amp_sii1,5)
    amp_sii2 = round(amp_sii2,5)

    fluxes = pd.DataFrame([data.lcgID[i], data.extension[i], data.ra[i], data.dec[i], amp_oiii1, amp_oiii2, amp_hbeta, amp_halpha, amp_nii1, amp_nii2, amp_oii1, amp_oii2, amp_sii1, amp_sii2]).T
    fluxes[0] = int(fluxes[0])
    fluxes[1] = int(fluxes[1])
    print('LCG' + str(data.lcgID[i]) + '_' + str(data.extension[i]))
    if flag_correction == 'on': fluxes.to_csv(output_dir+output_file, index=False, header=None, mode='a')
    if flag_correction == 'off': fluxes.to_csv(output_dir+output_file, index=False, header=None, mode='a')    

