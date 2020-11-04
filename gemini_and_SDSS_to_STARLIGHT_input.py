#Copyright 2013 Rogerio Riffel (riffel@ufrgs.br)
#
#   Licensed under the Apache License, Version 2.0 (the "License");
#   you may not use this file except in compliance with the License.
#   You may obtain a copy of the License at
#
#       http://www.apache.org/licenses/LICENSE-2.0
#
#   Unless required by applicable law or agreed to in writing, software
#   distributed under the License is distributed on an "AS IS" BASIS,
#   WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
#   See the License for the specific language governing permissions and
#   limitations under the License. 
#
#   Some parts of this script where written by other authors and in some cases 
#   modified by Rogerio Riffel. The original authors are #quoted in each routine.
#
import os, glob
from pylab import *
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.mlab as mlab
import re
from scipy import interpolate
import time
import os
from  scipy import ndimage
from numpy import sin,cos,round,isscalar,array,ndarray,ones_like,pi
from astropy.io.fits import open
from astropy.io import fits as pf
from astropy.table import Table

############################################################################
#
#
#      Convert coodinates from equatorial to galactic coordinates.
#
#     Author: Dorothee Brauer
#     Modified by Eduardo Balbinot
#     Modified by Rogerio Riffel, introducing a correction on the DE Calculation. Not fully tested yet
#
############################################################################

def hms2deg(RA,DE):
    '''
    Convert RA= [hh,mm,ss] and DEC = [DD,mm,ss] to degres.
    Usage: hms2deg(RA,DEC).
    Adapted by: Rogerio Riffel
    '''
    RA0 = array(RA).reshape(-1,3)
    DE0 = array(DE).reshape(-1,3)
    RA = 15.*RA0[:,0] + 15./60.*RA0[:,1] + 15./3600.*RA0[:,2]
    #      DE = DE0[:,0] + 1./60.*DE0[:,1] + 1./3600.*DE0[:,2]
    if DE0[:,0] >=0.:
        DE=((DE0[:,2]/60+DE0[:,1])/60 + DE0[:,0])
    elif DE0[:,0][0] <0.:
        DE=(-1*(DE0[:,2]/60+DE0[:,1])/60 + DE0[:,0])
    return RA,DE

def eq2galCoords(RA,DE,units='observers'):
    deg2rad = pi/180.
    rad2deg = 180./pi

    kpc2km = 3.085678e16
    yr2s   = 31557600.

    # Direction to the North Galactic Pole in Equatorial coordinates
    RA_GP = 15*(12.+49./60.+ 8.13e-4 *(2000.-1950.)) * deg2rad
    DE_GP =    (27.4       - 5.44e-3 *(2000.-1950.)) * deg2rad

    # Galactic longitude of the North Celestial Pole
    l_NCP = (123. - 1.33e-3 *(2000.-1950.)) * deg2rad

    if units == 'observers':
        (RA,DE)=hms2deg(RA,DE)

    if units == 'deg':
      RA = array(RA).reshape(-1)
      DE = array(DE).reshape(-1)

    sdp = sin(DE_GP)
    cdp = sqrt(1. - sdp*sdp)
    sdec= sin(DE*deg2rad)
    cdec= sqrt(1. - sdec*sdec)

    ras = RA*deg2rad - RA_GP
    sgb = sdec*sdp + cdec*cdp*cos(ras)

    b  = arcsin(sgb)*rad2deg

    cgb = sqrt(1. - sgb*sgb)
    sine= cdec*sin(ras) / cgb
    cose= (sdec-sdp*sgb) / (cdp*cgb)

    l = (l_NCP - arctan2(sine,cose))*rad2deg

#    lt0 = argwhere(l < 0.).reshape(-1)
#    l[lt0] = l[lt0]+360.

    return (l,b)

############################################################################
#
#
#      Reddening Laws
#
#      Written/addapted by Rogerio Riffel
#
############################################################################


def ccm(wavelength,rv):
    '''
    CCM -- Compute CCM Extinction Law for a given lambda and a given RV. Default is 3.1.
    
    Usage: ccm(wavelength,rv)
    Implemented by Rogerio Riffel
    '''

    # Convert to inverse microns
    x = 10000. / wavelength

    # Compute a(x) and b(x)
    if (x < 0.3):
        print("Wavelength out of range of extinction function")

    elif (x < 1.1):
        y = x ** 1.61
        a = 0.574 * y
        b = -0.527 * y

    elif (x < 3.3):
        y = x - 1.82
        a = 1 + y * (0.17699 + y * (-0.50447 + y * (-0.02427 +
        y * (0.72085 + y * (0.01979 + y * (-0.77530 + y * 0.32999))))))
        b = y * (1.41338 + y * (2.28305 + y * (1.07233 + y * (-5.38434 +
        y * (-0.62251 + y * (5.30260 + y * (-2.09002)))))))

    elif (x < 5.9):
        y = (x - 4.67) ** 2
        a = 1.752 - 0.316 * x - 0.104 / (y + 0.341)
        b = -3.090 + 1.825 * x + 1.206 / (y + 0.263)

    elif (x < 8.0):
        y = (x - 4.67) ** 2
        a = 1.752 - 0.316 * x - 0.104 / (y + 0.341)
        b = -3.090 + 1.825 * x + 1.206 / (y + 0.263)

        y = x - 5.9
        a = a - 0.04473 * y**2 - 0.009779 * y**3
        b = b + 0.2130 * y**2 + 0.1207 * y**3

    elif (x <= 10.0):
        y = x - 8
        a = -1.072 - 0.628 * y + 0.137 * y**2 - 0.070 * y**3
        b = 13.670 + 4.257 * y - 0.420 * y**2 + 0.374 * y**3

    else:
        print("\n\n >>>>> Lambda=",wavelength,": Wavelength out of range for the dust extintion law <<<<\n\n")

    # Compute A(lambda)/A(V)
    y = a + b / rv
    return (y)
#    print y

def ftzpt(lamb,flx,ebv):
    '''
    Fitzpatrick, Edward L (1999, PASP, 111, 63) extinction law 
    usage: ftzpt(lamb,flx,ebv)
    Implemented by: Daniel Faes
    Addapted by: Rogerio Riffel
    '''
    x = (1.e4)/lamb
    flx_cor = flx*10**(.4*ebv*(1.e-5 + .22707*x + 1.95243*x**2 - 2.67596*x**3 + 2.6507*x**4 - 1.26812*x**5 + 0.27549*x**6 - 0.02212*x**7))
    return(flx_cor)

def seaton(wavelength,ebv):
    '''
    Seaton (1979, MNRAS, 187P, 73) extinction law 
    Implemented by: Natacha Dametto
    Addapted by: Rogerio Riffel
    '''


    X=[1.0, 1.1 ,1.2 ,1.3 ,1.4 ,1.5 ,1.6 ,1.7 ,1.8 ,1.9 ,2.0 ,2.1 ,2.2 ,2.3 ,2.4 ,2.5 ,2.6 ,2.7 ]
    A=[1.36,1.44,1.84,2.04,2.24,2.44,2.66,2.88,3.14,3.36,3.56,3.77,3.96,4.15,4.26,4.40,4.52,4.64]
    x=10000.0/wavelength

    GetPts=interpolate.interp1d(X,A)
    try:
        y=GetPts(x)
    except:
        print("\n\n >>>>> Lambda=",wavelength,": Wavelength out of range for the dust extintion law <<<<\n\n")
    cor= y * ebv
    return (cor)

def calzetti(wavelength,ebv):
    '''
    Calzetti et all. (2000,  ApJ, 533, 682) extinction law.
    Implemented by: Natacha Dametto
    Addapted by: Rogerio Riffel
    '''

    x=wavelength/10000.0
    if x< 0.12:
        print("\n\n >>>>> Lambda=",wavelength,": Wavelength out of range for the dust extintion law <<<<\n\n")
    elif 0.12 <= x <= 0.63:
        k=2.659*(-2.156+(1.509/x)-(0.198/x**2)+(0.011/x**3)) + 4.05
    elif 0.63 <= x <= 2.20:
        k=2.659*(-1.857+(1.040/x)) + 4.05
    else:
        print("\n\n >>>>> Lambda=",wavelength,": Wavelength out of range for the dust extintion law <<<<\n\n")
    cor = k*ebv
    return (cor)

############################################################################
#
#
#      Routine to apply reddening Laws
#
#      Written by Rogerio Riffel
############################################################################

def deredden(spectra,law,ebv,rv=3.1,quiet=False):
    '''
    This routine allows to apply reddening corrections using the following reddening laws:
    
    - ccm - Cardelli, Clayton, and Mathis (1989, ApJ, 345, 245)
    - calzetti - Calzetti et all. (2000,  ApJ, 533, 682)
    - ftzpt - Fitzpatrick, Edward L (1999, PASP, 111, 63)
    - seaton - Seaton (1979, MNRAS, 187P, 73).
    
    usage: deredden(spectra,law,ebv,rv)
    spectra = file name of an ASCII file with 2 or 3 collumns (lambda, flux, eflux) or an array with [lambda, flux and/or eflux]. If there is no eflux, it is assumed as zero, and returns eflux=0 in both cases (array or file name).
    law= one of the above laws 
    ebv = E(B-V)
    rv = 3.1 is the default.
    example: deredden('spectra.dat','ccm',ebv=0.5,rv=3.0)
    
    
    '''
    ebv=ebv
    law=law
    rv=rv
    av=rv*ebv
    if type(spectra) == str:
        filein=spectra
        corrname=spectra
        try:
            (lambdas,flux,eflux)=np.loadtxt(filein, unpack=True)
        except:
            (lambdas,flux)=np.loadtxt(filein, unpack=True)
            eflux=flux*0.0
# Trying to implement the correction for a color, but it does not work in the loops.
#    elif type(spectra) == list:
#        spectra=array(spectra).reshape(-2)
#        try:
#            lambdas=spectra[0]
#            flux=spectra[1]
#            eflux=spectra[2]
#        except:
#            lambdas=spectra[0]
#            flux=spectra[1]
#            eflux=flux*0.0
#        corrname='a list of colors'
    else:
        try:
            lambdas=spectra[:,0]
            flux=spectra[:,1]
            eflux=spectra[:,2]
        except:
            lambdas=spectra[:,0]
            flux=spectra[:,1]
            eflux=flux*0.0
        corrname='an array'
    fluxCor=[]
    efluxCor=[]
    if (law == 'ftzpt'):
        for i in range(0,len(lambdas)):
            fluxCor.append(ftzpt(lambdas[i],flux[i],ebv))
            efluxCor.append(ftzpt(lambdas[i],eflux[i],ebv))
    
    elif(law == 'ccm'):
        for i in range(0,len(lambdas)):
            cor = 10. ** (0.4 * av * ccm(lambdas[i], rv))
            fluxCor.append(flux[i] * cor)
            efluxCor.append(eflux[i] * cor)

    elif(law == 'seaton'):
        for i in range(0,len(lambdas)):
            cor = 10. ** (0.4 *seaton(lambdas[i], ebv))
            fluxCor.append(flux[i] * cor)
            efluxCor.append(eflux[i] * cor)
    elif(law== 'calzetti'):
        for i in range(0,len(lambdas)):
            cor = 10. ** (0.4 *calzetti(lambdas[i], ebv))
#            print lambdas[i], flux[i], flux[i]*cor, calzetti(lambdas[i], ebv)
            fluxCor.append(flux[i] * cor)
            efluxCor.append(eflux[i] * cor)

    if quiet == True:
        print("\n ----------------------------------------------------\n")
        print("Reddening correction sucessfull applied to: ", corrname)
        print("\nUsing: Law=",law,"; E(B-V)=",ebv," & Rv= ",rv)
        print("\n ----------------------------------------------------\n")
    return np.array(lambdas),np.array(fluxCor),np.array(efluxCor)



############################################################################
#
#
#      Routine to get ebv from Schlegel et al. dust maps http://irsa.ipac.caltech.edu/applications/DUST/
#      Copyright 2008 Erik Tollerud
############################################################################

#Copyright 2008 Erik Tollerud
#
#   Licensed under the Apache License, Version 2.0 (the "License");
#   you may not use this file except in compliance with the License.
#   You may obtain a copy of the License at
#
#       http://www.apache.org/licenses/LICENSE-2.0
#
#   Unless required by applicable law or agreed to in writing, software
#   distributed under the License is distributed on an "AS IS" BASIS,
#   WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
#   See the License for the specific language governing permissions and
#   limitations under the License. 
# 
#   This file was modified by Eduardo Balbinot (balbinot@if.ufrgs.br)
#   in 2012 for the purpose of using it in the QR tool developed by LIneA. 
#
#   This file was modified by Rogerio Riffel (riffel@ufrgs.br)
#   in 2013 for the purpose of using it in the ReddeningCorrections tool.


def get_SFD_dust(longitud,lat,dustmap='ebv',pathdustmap='./aux_files/DustMaps',interpolate=True):
    """
    Gets map values from Schlegel, Finkbeiner, and Davis 1998 extinction maps.
    
    `dustmap` can either be a filename (if '%s' appears in the string, it will be
    replaced with 'ngp' or 'sgp'), or one of:
    
    * 'i100' 
        100-micron map in MJy/Sr
    * 'x'
        X-map, temperature-correction factor
    * 't'
        Temperature map in degrees Kelvin for n=2 emissivity
    * 'ebv'
        E(B-V) in magnitudes
    * 'mask'
        Mask values 
        
    For these forms, the files are assumed to lie in the current directory.
    
    Input coordinates are in degrees of galactic latiude and logitude - they can
    be scalars or arrays.
    
    if `interpolate` is an integer, it can be used to specify the order of the
    interpolating polynomial

    """
    
    if type(dustmap) is not str:
        raise ValueError('dustmap is not a string')
    dml=dustmap.lower()
    if dml == 'ebv' or dml == 'eb-v' or dml == 'e(b-v)' :
        dustmapfn=pathdustmap+'/SFD_dust_4096_%s.fits'
    else:
        dustmapfn=dustmap
    
    if isscalar(longitud):
        l=array([longitud])*pi/180
    else:
        l=array(longitud)*pi/180
    if isscalar(lat):
        b=array([lat])*pi/180
    else:
        b=array(lat)*pi/180
        
    if not len(l)==len(b):
        raise ValueError('input coordinate arrays are of different length')
    
    if '%s' not in dustmapfn:
        f=open(dustmapfn)
        try:
            mapds=[f[0].data]
        finally:
            f.close()
        assert mapds[-1].shape[0] == mapds[-1].shape[1],'map dimensions not equal - incorrect map file?'
        
        polename=dustmapfn.split('.')[0].split('_')[-1].lower()
        if polename=='ngp':
            n=[1]
            if sum(b > 0) > 0:
                print('using ngp file when lat < 0 present... put %s wherever "ngp" or "sgp" should go in filename')
        elif polename=='sgp':
            n=[-1]
            if sum(b < 0) > 0:
                print('using sgp file when lat > 0 present... put %s wherever "ngp" or "sgp" should go in filename')
        else:
            raise ValueError("couldn't determine South/North from filename - should have 'sgp' or 'ngp in it somewhere")
        masks = [ones_like(b).astype(bool)]
    else: #need to do things seperately for north and south files
        nmask = b >= 0
        smask = ~nmask
        
        masks = [nmask,smask]
        ns = [1,-1]
        
        mapds=[]
        f=open(dustmapfn%'ngp')
        try:
            mapds.append(f[0].data)
        finally:
            f.close()
        assert mapds[-1].shape[0] == mapds[-1].shape[1],'map dimensions not equal - incorrect map file?'
        f=open(dustmapfn%'sgp')
        try:
            mapds.append(f[0].data)
        finally:
            f.close()
        assert mapds[-1].shape[0] == mapds[-1].shape[1],'map dimensions not equal - incorrect map file?'
    
    retvals=[]
    for n,mapd,m in zip(ns,mapds,masks):
        #project from galactic longitude/latitude to lambert pixels (see SFD98)
        npix=mapd.shape[0]
        
        x=npix/2*cos(l[m])*(1-n*sin(b[m]))**0.5+npix/2-0.5
        y=-npix/2*n*sin(l[m])*(1-n*sin(b[m]))**0.5+npix/2-0.5
        #now remap indecies - numpy arrays have y and x convention switched from SFD98 appendix
        x,y=y,x
        
        if interpolate:
            from scipy.ndimage import map_coordinates
            if type(interpolate) is int:
                retvals.append(map_coordinates(mapd,[x,y],order=interpolate))
            else:
                retvals.append(map_coordinates(mapd,[x,y]))
        else:
            x=round(x).astype(int)
            y=round(y).astype(int)
            retvals.append(mapd[x,y])
            
    if isscalar(longitud) or isscalar(lat):
        for r in retvals:
            if len(r)>0:
                return r[0]
        assert False,'None of the return value arrays were populated - incorrect inputs?'
    else:
        #now recombine the possibly two arrays from above into one that looks like  the original
        retval=ndarray(l.shape)
        for m,val in zip(masks,retvals):
            retval[m] = val
        return retval


def linereddening(line1,line2,obs_ratio,theo_ratio,error_ratio_1,error_ratio_2,law,rv=3.1,verbose='y'):
    '''
    This routine allows to apply reddening corrections using the following reddening laws:
    
    - ccm - Cardelli, Clayton, and Mathis (1989, ApJ, 345, 245)
    - calzetti - Calzetti et all. (2000,  ApJ, 533, 682)
    - ftzpt - Fitzpatrick, Edward L (1999, PASP, 111, 63)
    - seaton - Seaton (1979, MNRAS, 187P, 73).
    
    usage: linereddening(line1,line2,ObsRatio,TeorRatio,law,rv=3.1)
    line1 = Wavelength of the first line
    line2 = Wavelength of the second line
    obs_ratio = line1/line2 observed ratio
    theo_ratio = line1/line2 theoretical ratio
    error_ratio_1= Flux error / Flux for line1
    error_ratio_2= Flux error / Flux for line2
    law= one of the above laws 
    rv = 3.1 is the default.
    verbose= y/n to print or not the input informations
    '''
    if (law == 'ftzpt'):
       f_lamb=(ftzpt(line1,rv) - ftzpt(line2,rv))
    elif(law == 'ccm'):
       f_lamb=(ccm(line1,rv) - ccm(line2,rv))
    elif(law == 'seaton'):
       f_lamb=(seaton(line1,rv) - seaton(line2,rv))
    elif(law== 'calzetti'):
       f_lamb=(calzetti(line1,rv) - calzetti(line2,rv))

    C_ext=np.log10(obs_ratio/theo_ratio)/(-1.0*f_lamb)
     
    sigR=np.sqrt((obs_ratio)**2*((error_ratio_1)**2 + (error_ratio_2)**2)) # error propagation
    Cext_err = abs(1./(f_lamb*(obs_ratio))*np.log10(np.e)*sigR)
    
    if verbose == 'y':
      toprint=('\n*********\nline1='+str(line1)+'\n'+'line2='+str(line2)+'\n'+'obs_ratio='+str(obs_ratio)+'\n' 
+'theo_ratio='+str(theo_ratio)+'\n'+'error_ratio_1='+str(error_ratio_1)+'\n'+'error_ratio_2='+str(error_ratio_2)+'\n' 
+'law='+law +'\n'+'rv='+str(rv)+'\n-----\n'+'C_ext='+str(C_ext)+'\n'+'Cext_err='+str(Cext_err)+'\n*********\n')
      print(toprint)
    return C_ext, Cext_err,f_lamb


###########################################################
#                  Interpolation: spline

# Convertido de R (Marina T.) para Python
###########################################################

def interp(xin,xout,yin):
    Ninterpol = len(xout)
    yout = np.array([])
    
    for k in range(Ninterpol): 
        t = xin[xin < xout[k]]
        tind = len(t)-1 #Verificar este -1. Soluciona o problema, mas pode ter gerado outros escondidos
        
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

##############################################################################################
#                                    sdss2tocorrect_sdss
##############################################################################################
# PURPOSE
#    Convert sdss .fits files to .csv files that are input for extinction and doppler correction codes
#
# CALLING SEQUENCE
#    sdss2tocorrect_sdss('file_name')
#
# INPUT PARAMETERS
#    String 'path/file_name' that is the file name containig a list of groups and extensions of sdss data separated by a space
#    
# OUTPUT
#   tocorrect_sdss_LCG %group% _ %extention% .csv 
#
# REQUIRED SCRIPTS
#   
# COMMENTS
# Developed by Vitor Bootz
##############################################################################################

def sdss2tocorrect_sdss(file_list_name): #Nome do arquivo contendo a lista de objetos no modo "group extension"
    data = pd.read_csv(file_list_name, sep=',')
    selection = (data['flag_sdss'] == 1) & (data['onoff'] == 1)
    data = data[selection]
    data.index = range(len(data))

    for i in range(len(data)):
        group = str(data.lcgID[i])
        extension = str(data.extension[i])
        hdul = pf.open('files_sdss/' + group + '_' + extension + '.fits', memmap=True)
        evt_data = Table(hdul[1].data)
        try:
            csv_file = pd.DataFrame([10**evt_data['loglam'], evt_data['flux'], 1/np.sqrt(evt_data['ivar'])]).T
        except:
            csv_file = pd.DataFrame([10**evt_data['LOGLAM'], evt_data['FLUX'], 1/np.sqrt(evt_data['IVAR'])]).T
        csv_file.columns = ['lambda', 'flux', 'error']

        selection = csv_file['lambda'] == float('inf')
        selection += csv_file['error'] == float('inf')
        dropinf = csv_file[selection]
        for i in range(len(dropinf)):
            item = dropinf.index[i]
            csv_file = csv_file.drop(item)

        csv_file.to_csv('files_tocorrect/tocorrect_sdss_LCG' + group + '_' + extension + '.csv', header=False, index=False, sep=' ')
        print('tocorrect_sdss_LCG' + group + '_' + extension + '.csv')
        
        
##############################################################################################
#                                    csv_gemini2tocorrect_gemini
##############################################################################################
# PURPOSE
#    Convert gemini .csv files to .csv files that are input for extinction and doppler correction codes
#
# CALLING SEQUENCE
#    csv_gemini2tocorrect_gemini('file_name')
#
# INPUT PARAMETERS
#    String 'file_name' that is the file name containig a list of groups and extensions of sdss data separated by a space
#    
# OUTPUT
#   tocorrect_gemini_LCG %group% _ %extention% .csv 
#
# REQUIRED SCRIPTS
#   
# COMMENTS
# Developed by Vitor Bootz
##############################################################################################
def csv_gemini2tocorrect_gemini(file_list_name):
    spectra_z_file = file_list_name
    data = pd.read_csv(spectra_z_file, sep=',')
    selection = (data['flag_sdss'] == 0) & (data['onoff'] == 1)
    data = data[selection]
    data.index = range(len(data))

    for i in range(len(data)):
        group = str(data.lcgID[i])
        extension = str(data.extension[i])
        name_spec = 'files_gemini/csv_' + str(group) + '/csv_LCG' + str(group) + '_' + str(extension) + '.csv'
        name_erro = 'files_gemini/csv_' + str(group) + '/csv_ERROR_LCG' + str(group) + '_' + str(extension) + '.csv'
        spec = pd.read_csv(name_spec, header=None, delim_whitespace=True, names = ['lambda', 'flux'])
        erro = pd.read_csv(name_erro, header=None, delim_whitespace=True, names = ['lambda', 'error'])
        output = pd.DataFrame([spec['lambda'], spec['flux'], erro['error']]).T
        output.to_csv('files_tocorrect/tocorrect_gemini_LCG' + str(group) + '_' + str(extension) + '.csv', index=False, header=False, sep=' ')
        print('tocorrect_gemini_LCG' + str(group) + '_' + str(extension) + '.csv')


##############################################################################################
#                                    radec_sdss_organizer
##############################################################################################
# PURPOSE
# In case of missing "allgalaxies_list" file, these commands will collect and organize all the useful galaxies of "lista_all_galaxies" (result from a sql query in casjobs)   
#
# CALLING SEQUENCE
#    radec_sdss_organizer('file_name')
#
# INPUT PARAMETERS
#    String 'file_name' that is the file name containig the result of casjobs query called "lista_all_galaxies"
#    
# OUTPUT
#   allgalaxies_radec.csv, containing columns [lcgID, ra, dec] with no header
#
# REQUIRED SCRIPTS
#   
# COMMENTS
# Developed by Vitor Bootz
##############################################################################################
def radec_sdss_organizer(allgalaxies_info):
    dados = pd.read_csv('aux_files/allgalaxies_info.csv')
    radec_lcg = dados[['lcgID','ra_lcg','dec_lcg']]
    radec_others = dados[['lcgID','ra','dec']]
    radec_lcg.index = dados.lcgID
    radec_others.index = dados.lcgID
    radec_lcg.columns = [['lcgID','ra', 'dec']]
    radec_others.columns = [['lcgID','ra', 'dec']]
    radec_lcg.drop_duplicates(inplace = True)
    radec = radec_lcg.append(radec_others)
    radec.sort_index(inplace=True)
    radec = radec[['ra','dec']]
    radec.to_csv('allgalaxies_radec.csv', sep=',', header=False)
        
######################################################################################################################
######################################################################################################################

#                                   ISOLATED GROUPS OF EXTREMELY BLUE DWARF GALAXIES

######################################################################################################################
######################################################################################################################

###########################################################################
#   Nome do arquivo contendo a lista de galáxias a ser analisada
#   Chamando funções de transformação de arquivos sdss e gemini para .csv
###########################################################################
spectra_z_file = 'aux_files/galaxy_list.csv'
sdss2tocorrect_sdss(spectra_z_file)
csv_gemini2tocorrect_gemini(spectra_z_file)


##################################################################################
#       Puxando coords dos arquivos de interesse e convertendo de eq para gal
##################################################################################

data = pd.read_csv(spectra_z_file, sep=',')
selection = data['onoff'] == 1
data = data[selection]
data.index = range(len(data))
radec = pd.DataFrame([data.lcgID, data.ra, data.dec]).T
radec.index = range(len(radec))
radec.columns = ['lcgID', 'ra', 'dec']
allcoords = pd.DataFrame(eq2galCoords(radec[['ra']], radec[['dec']], units='deg')).T
allcoords = pd.concat([radec,allcoords], axis = 1)
allcoords.columns = ['lcgID', 'ra','dec','l','b']

###########################################################
#                  Aplicando E(B-V)
###########################################################
ebv_maps = pd.DataFrame(get_SFD_dust(allcoords.l,allcoords.b,dustmap='ebv', pathdustmap='./aux_files/DustMaps',interpolate=True))
ebv_maps.columns = ['ebv']
allcoords = pd.concat([allcoords,ebv_maps], axis=1)

######################################################################################################################
#       Selecionando amostra Gemini e organizando
# Os slits das lcgs das galáxias do sloan e 6193 são == 1
# Para alterar os arquivos de interesse, basta alterar a flag "onoff" no arquivo "galaxy_list.csv"
######################################################################################################################
gemini = pd.DataFrame([data.extension, allcoords.ra, allcoords.dec, allcoords.l, allcoords.b, allcoords.ebv, data.z, data.flag_sdss]).T
gemini.index = allcoords.lcgID.round(0).astype(int)
gemini.columns = ['slit','ra','dec','l','b','ebv','redshift','flag_sdss']
gemini.slit = gemini.slit.round(0).astype(int)
gemini.flag_sdss = gemini.flag_sdss.round(0).astype(int)

###########################################################
#             Aplicando extinção galáctica
###########################################################

all_files = pd.DataFrame([gemini.index, gemini.slit, gemini.flag_sdss]).T
all_files.columns = ['group', 'extension', 'flag_sdss']
all_files
for i in all_files.index:
    group = int(all_files.group[[i]])
    extension = int(all_files.extension[[i]])
    if all_files.flag_sdss[i] == 1:
        name_spec = 'files_tocorrect/tocorrect_sdss_LCG' + str(group) + '_' + str(extension) + '.csv'
    else:
        name_spec = 'files_tocorrect/tocorrect_gemini_LCG' + str(group) + '_' + str(extension) + '.csv'
    group_data = gemini.loc[[group]]
    selection = group_data['slit'] == extension
    ebv = float(group_data[selection].ebv)
    spectra = np.loadtxt(name_spec, delimiter=' ', unpack=True).T
    extincted_spectra = pd.DataFrame(deredden(spectra,'ccm',ebv,rv=3.1,quiet=False)).T
    #extincted_spectra.to_csv('extincted_LCG' + str(group) + '_' + extension + '.csv', index=False, header=False)


###########################################################
#           Convertendo lambda para rest-frame
###########################################################

    #name_spec = 'extincted_LCG' + str(group) + '_' + extension + '.csv'
    #dopcor_spectra = pd.read_csv(name_spec, header=None)
    dopcor_spectra = extincted_spectra
    spectra = dopcor_spectra[0] / (float(group_data[selection].redshift) + 1) # lambda / (z+1)
    dopcor_spectra[0] = spectra
    dopcor_spectra.columns = ['lambda', 'flux', 'error']
    #dopcor_spectra.to_csv('dopcor_extincted_LCG' + str(group) + '_' + extension + '.csv', index=False, header=False)
    #dopcor_spectra


###########################################################
#                   Input da interpolação
###########################################################
    lambda_input = np.array(dopcor_spectra['lambda']) # Lambda original
    flux_input = np.array(dopcor_spectra['flux']) # Fluxo original
    error_input = np.array(dopcor_spectra['error']) # Error original
    newlambda_input = np.array([])
    for j in range(len(dopcor_spectra['lambda'])): 
        newlambda_input = np.append(newlambda_input, int(dopcor_spectra['lambda'][j])) # Lambda de interesse

    
###########################################################
#            Execução da interpolação e output
###########################################################
    newflux = pd.DataFrame(interp(lambda_input, newlambda_input, flux_input))
    newerror = pd.DataFrame(interp(lambda_input, newlambda_input, error_input))
    newlambda = pd.DataFrame(newlambda_input)
    output = pd.concat([newlambda,newflux,newerror], axis=1)
    output.columns = ['lambda', 'flux', 'error']
    if int(all_files.flag_sdss[[i]]) == 0:
        output[['flux']] = output[['flux']] * 1e17
        output[['error']] = output[['error']] * 1e17
    output.drop_duplicates(inplace = True)
    if all_files.flag_sdss[i] == 1:
        output.to_csv('files_input_starlight/input_starlight_sdss_LCG' + str(group) + '_' + str(extension) + '.csv', sep=' ',index=False, header=False)
        print('file_output: ' + 'input_starlight_sdss_LCG' + str(group) + '_' + str(extension) + '.csv')
    else:
        output.to_csv('files_input_starlight/input_starlight_gemini_LCG' + str(group) + '_' + str(extension) + '.csv', sep=' ',index=False, header=False)
        print('file_output: ' + 'input_starlight_gemini_LCG' + str(group) + '_' + str(extension) + '.csv')
