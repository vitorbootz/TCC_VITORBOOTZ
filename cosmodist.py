
##############################################################################################
# PURPOSE
#   Cosmological distance calculator function
#
# CREATED BY:
#   Author: Marina Trevisan (in R)
#
# ADAPTED BY:
#   Vitor Eduardo Buss Bootz (conversion to Python)
#
# CALLING SEQUENCE
#   python cosmodist.py     --> In terminal
#
# INPUT PARAMETERS
#   z --> Redshift
#   H0 --> Hubble constant
#   Omega_m --> Matter density parameter
#   Omega_l --> Cosmological constant
#    
# OUTPUT
#   DH_Mpc, tH_Gyr, DM_Mpc, DA_Mpc, DL_Mpc, VC_Gpc3, TL_Gyr, DM
#
# REQUIRED SCRIPTS
#   
# COMMENTS    
#
##############################################################################################

import numpy as np
import pandas as pd

def cosmodist(z, H0, Omega_m, Omega_l, dz=1e-6): #Isto torna como padrão o dz=1e-6, transformando-a em parâmetro não obrigatório
    c = 299792
    DH_Mpc = c / H0     # Mpc
    tH_Gyr = 978 / H0   # Gyr

    if z > 0:
        zz = np.linspace(start = 0, stop = z, num = int(z/dz)+1)

        E_z = np.sqrt(Omega_m * (1 + zz)**3 + Omega_l)

        DC_Mpc = DH_Mpc * sum(dz / E_z) #Esta divisão só funciona porque zz é array numpy
        DM_Mpc = DC_Mpc       # Omega_z = 0 (flat Universe)

        DA_Mpc = DM_Mpc / (1 + z)
        DL_Mpc = (1 + z) * DM_Mpc

        TL_Gyr = tH_Gyr * sum(dz / ((1 + zz) * E_z))

        DM = 5 * (np.log10(DL_Mpc * 1e6) - 1)

        VC_Gpc3 = 4 * np.pi / 3 * DM_Mpc**3 * 1e-9

        cosmo = pd.DataFrame([DH_Mpc, tH_Gyr, DM_Mpc, DA_Mpc, DL_Mpc, VC_Gpc3, TL_Gyr, DM]).T
    else:
        DM_Mpc = 0
        DA_Mpc = 0
        DL_Mpc = 0
        VC_Gpc3 = 0
        TL_Gyr = 0
        DM = 0
        cosmo = pd.DataFrame([DH_Mpc, tH_Gyr, DM_Mpc, DA_Mpc, DL_Mpc, VC_Gpc3, TL_Gyr, DM]).T
            
    cosmo.columns = ['DH_Mpc', 'tH_Gyr', 'DM_Mpc', 'DA_Mpc', 'DL_Mpc', 'VC_Gpc3', 'TL_Gyr', 'DM']

    return(cosmo)

