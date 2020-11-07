##################################################################################################
# PURPOSE:
#   Calculate photometric magnitudes for passbands. 
# CALLING SEQUENCE:
#    python calmag.py
# INPUTS:
#    galaxy_list.csv
# PARAMETERS:
#    
# OUTPUT:
#    Magnitude of each galaxy from galaxy_list.csv
# REQUIRED SCRIPTS & DATA:
#   "softening"
#   "vega.out"
#   "vega.out"
#   "Filters_sdss/"
##################################################################################################
import numpy as np
import pandas as pd
from os import path
import time
from astropy.io import fits as pf
from astropy.table import Table

list = dir() #Seriam estas duas linhas o equivalente?
del list

#################################
# DEFINITIONS
#################################
sdss_band = "r"  # g, r, i, z # r have the best response between the SDSS bands

spectra_z_file = '../aux_files/galaxy_list.csv'
sample = pd.read_csv(spectra_z_file, sep=',')
data = pd.read_csv(spectra_z_file, sep=',')
selection = (data['onoff'] == 1) #Selecting only Gemini galaxies
data = data[selection]
data.index = range(len(data))

#################################
# FUNCTIONS
#################################
# --------------------------
# Function for interpolation
# --------------------------
def interp(xin,xout,yin):
    Ninterpol = len(xout)
    yout = np.array([])
    
    for k in range(Ninterpol): 
        t = xin[xin < xout[k]]
        tind = len(t) -2 #VERIFICAR SE ESTE -2 NÃO INTERFERE NOS RESULTADOS
        
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

# ------------------------------------
# Function for computing the magnitude
# ------------------------------------
def getmag(Flux, band, fvega):
    vega_abmag = pd.DataFrame([-0.08, 0.17, 0.40, 0.57]) #References? I founded another table of values
    #vega_abmag = pd.DataFrame([-0.08, 0.16, 0.37, 0.54])
    f0 = fvega
    t = pd.read_csv("softening", delim_whitespace=True)
    t.columns = ("V1","V2","V3","V4","V5")
    
    b = t.V2[t.V1 == band] # "Softening parameter: typical 1-sigma noise of the sky in a PSF aperture in 1" seeing"
    f0b = f0[bands == band] # Vega Flux
    
    mag = (-2.5 / np.log(10)) * (np.arcsinh((Flux/f0b) / (2*b)) + np.log(b)) + float(vega_abmag[bands == band][0])
    #The equation above can be accessed in https://www.sdss.org/dr12/algorithms/magnitudes/
    #The use of arcsinh is an usefull techinque to report a magnitude even in the absence of a formal detection, since negative flux fails in the log way
    
    return(mag)


#################################
# READ FILTER
#################################
t = pd.read_csv("Filters_sdss/"+ sdss_band + ".dat", index_col=False, delim_whitespace=True, skiprows=6, header=None)
t.columns = ("V1","V2","V3","V4","V5")
filter_l = np.array(t.V1) #Lambda filter
filter_curve = np.array(t.V4) #Filter response in sdss_band band (a coefficient)
filter_lmin = filter_l[0] #Minimum lambda with response > 0 in the filter
filter_lmax = filter_l[len(filter_l)-1] #Maximum lambda with response > 0 in the filter

dl = 0.1 #step
filter_lout = np.linspace(filter_lmin, filter_lmax, num = int((filter_lmax-filter_lmin)/dl+1))
filter_lout = np.round(filter_lout, 1)
#By increasing the number of lambda intervals, we decrase the error in the "integration" (in fact, is just a sum of all these small regions).
#The higher the value of filter_lout, more precise the integration will be

filter_curveout = interp(filter_l, filter_lout, filter_curve) #Output of the interpolation of filter_curve of filter_l to filet_lout

#################################
# READ VEGA
#################################
t = pd.read_csv("vega.out", skiprows=4, delim_whitespace=True, index_col=False, header=None)
t.columns = ("V1","V2","V3")

lvega = np.array(t.V2) #Lambda of Vega
fvega = np.array(t.V3) #Flux of Vega
fvegaout = interp(lvega, filter_lout, fvega) #Interpolated flux of Vega in filter_lout lambda interval

#nu1 = 2.99792e8 / ((4702.5 - 928.17/2) * 1e-10)
#nu2 = 2.99792e8 / ((4702.5 + 928.17/2) * 1e-10)
#dnu = abs(nu1 - nu2)

flux_vega = sum(fvegaout * filter_curveout) * dl #Final flux vega

bands = np.array(["g","r","i","z"])
filter_vega = np.array([-99.0,-99.0,-99.0,-99.0]) #Just storage
filter_vega[bands == sdss_band] = flux_vega

mag_vega = getmag(flux_vega, sdss_band, filter_vega) #Magnitude of vega

#####################################################
for i in range(len(data)):
    if data.flag_sdss[i] == 0:
        input_spec = "../files_gemini/csv_" + str(data.lcgID[i]) + "/csv_LCG" + str(data.lcgID[i]) + "_" + str(data.extension[i]) + ".csv"
        t = pd.read_csv(input_spec, delim_whitespace=True, header=None)
    if data.flag_sdss[i] == 1:
        hdul_sdss = pf.open('../files_sdss/' + str(data.lcgID[i]) + '_' + str(data.extension[i]) + '.fits', memmap=True)
        evt_sdss = Table(hdul_sdss[1].data)
        try:
            sdss = pd.DataFrame([10**evt_sdss['loglam'], evt_sdss['flux'], 1/np.sqrt(evt_sdss['ivar'])]).T
        except:
            sdss = pd.DataFrame([10**evt_sdss['LOGLAM'], evt_sdss['FLUX'], 1/np.sqrt(evt_sdss['IVAR'])]).T
        sdss.columns = ['lambda', 'flux', 'error']
        sdss['flux'] = sdss['flux'] / 1e17
        t = sdss[['lambda', 'flux']]

    t.columns = ("V1","V2")
    filter_obsout = interp(t.V1, filter_lout, t.V2) #Final flux of the target in filter_lout lambda interval
    mag_filter = getmag(sum(filter_obsout * filter_curveout) * dl, sdss_band, filter_vega) #Final magnitude of the target in the sdss_band band
    
#################################
# SAVING RESULTS INTO galaxy_list.csv (main file)
#################################

    data.magSlit_r[i] = round(mag_filter,5)
    print("LCG" + str(data.lcgID[i]) + "_" + str(data.extension[i]) + ": " + "magnitude in the " + sdss_band + " band: " + str(float(mag_filter)))
    
for k in range(len(data)):
    for j in range(len(sample)):
        if (sample.ra[j] == data.ra[k]) & (sample.extension[j] == data.extension[k]):
            sample.magSlit_r[j] = data.magSlit_r[k]
            
sample.to_csv('../aux_files/galaxy_list.csv', index=False)
