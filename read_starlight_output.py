##############################################################################################
# PURPOSE
#   Read STARLIGHT output files
# 
# CREATED BY:
#   Marina Trevisan (in R)
#
# ADAPTED BY:
#   Vitor Eduardo Buss Bootz (Conversion from R to Python and adaptation)
#
# CALLING SEQUENCE
#   python read_starlight_output.py     --> In terminal
#
# INPUT PARAMETERS
#   output_dir            --> Output directory name
#   output_file           --> Output filename
#   starlight_output_dir  --> Directory of input files
#   spectra_z_file        --> List of objects
#   BaseAges              --> Base file
#   DR = 12               --> Spectra file name: DR12 or DR7 format
#   age_format            --> Units of ages in the base file (yr or Gyr)
#   met_format = 'Zsun'   --> Units of metallicitues in the base file (logZ_Zsun or Z)
#   Z_sun = 0.019         --> Solar metallicity
#    
# OUTPUT
#   STARLIGHT_MAIN_RESULTS.csv
#   Starlight_SFH_Mcor.csv
#   Starlight_SFH_Mcor_cumul.csv
#   Starlight_SFH_Mini.csv
#   Starlight_SFH_Mini_cumul.csv
#   Starlight_CEH_Lum.csv
#   Starlight_CEH_Mcor.csv
#   Starlight_CEH_Mini.csv
#
# REQUIRED SCRIPTS
# cosmodist.py
#   
# COMMENTS
#
##############################################################################################

##############################################################################################
# INPUTS
##############################################################################################
import numpy as np
import pandas as pd
from os import path
import time



def read_starlight_sfh(starlight_file):
    met_L_age = pd.DataFrame([np.zeros(n_age)]).T
    met_M_age = pd.DataFrame([np.zeros(n_age)]).T
    alpha_L_age = pd.DataFrame([np.zeros(n_age)]).T
    alpha_M_age = pd.DataFrame([np.zeros(n_age)]).T
    frac_age = pd.DataFrame([np.zeros(n_age)]).T
    mini_frac_age = pd.DataFrame([np.zeros(n_age)]).T
    mcor_frac_age = pd.DataFrame([np.zeros(n_age)]).T

    for aa in range(n_age):
        xx = abs(age/1e9 - ages[aa]/1e9) < 0.001
        frac_age[0][aa] = sum(frac[xx])
        mini_frac_age[0][aa] = sum(mini_frac[xx])
        mcor_frac_age[0][aa] = sum(mcor_frac[xx])
        if frac_age[0][aa] > 0:
            met_L_age[0][aa] = (np.log10(sum(frac[xx] * met[xx]/0.019) / sum(frac[xx])))
            alpha_L_age[0][aa] = (np.log10(sum(frac[xx] * 10**alpha[xx]) / sum(frac[xx])))

            met_M_age[0][aa] = (np.log10(sum(mcor_frac[xx] * met[xx]/0.019) / sum(mcor_frac[xx])))
            alpha_M_age[0][aa] = (np.log10(sum(mcor_frac[xx] * 10**alpha[xx]) / sum(mcor_frac[xx])))

    mcor_frac_age = mcor_frac_age/sum(mcor_frac_age[0])
    mini_frac_age = mini_frac_age/sum(mini_frac_age[0])

### Cumulative ###
    mini_cumul_age = pd.DataFrame([np.zeros(n_age)]).T
    mcor_cumul_age = pd.DataFrame([np.zeros(n_age)]).T
    for aa in range(n_age):
        mini_cumul_age[0][aa] = sum(mini_frac_age[0][ages >= ages[aa]])
        mcor_cumul_age[0][aa] = sum(mcor_frac_age[0][ages >= ages[aa]])

### Write the output files ###
    mini_frac_age = mini_frac_age.T
    mini_cumul_age = mini_cumul_age.T
    mcor_frac_age = mcor_frac_age.T
    mcor_cumul_age = mcor_cumul_age.T

    mini_frac_age.to_csv(output_dir+output_file_Mini, mode='a', sep=' ', index=False, header=False)
    mini_cumul_age.to_csv(output_dir+output_file_Mini_cumul, mode='a', sep=' ', index=False, header=False)
    mcor_frac_age.to_csv(output_dir+output_file_Mcor, mode='a', sep=' ', index=False, header=False)
    mcor_cumul_age.to_csv(output_dir+output_file_Mcor_cumul, mode='a', sep=' ', index=False, header=False)
    

list = dir()
del list 

### FILES & FOLDERS ###
starlight_output_dir = '../files_output_starlight/'                                        # Directory of input files
spectra_z_file = '../aux_files/galaxy_list.csv'                                   # Sufix of input file name
BaseAges = pd.read_csv('../starlight/csv_Base_MILES15_v10.0', sep=',', header=None)     # Base file                                                                  # Spectra file name: DR12 or DR7 format
output_dir = '../results_starlight/'

### OUTPUT FILES ###
output_file = 'STARLIGHT_MAIN_RESULTS.csv'   # Output file

output_file_Mini_SFH       = 'Starlight_SFH_Mini.out'                 #Output File
output_file_Mini_cumul_SFH = 'Starlight_SFH_Mini_cumul.out'           #Output File
output_file_Mcor_SFH       = 'Starlight_SFH_Mcor.out'                 #Output File
output_file_Mcor_cumul_SFH = 'Starlight_SFH_Mcor_cumul.out'           #Output File

output_file_Mini_CEH = 'Starlight_CEH_Mini.out'            # Output file
output_file_Mcor_CEH = 'Starlight_CEH_Mcor.out'            # Output file
output_file_Lum_CEH  = 'Starlight_CEH_Lum.out'             # Output file

### UNITS OF MEASUREMENT ###
age_format = 'yr'         # yr or Gyr
met_format = 'Z'          # logZ_Zsun or Z
Z_sun = 0.019
DR = 12
flag_sdss = 0 # Apenas inicializando uma variável

exec(open("../aux_files/cosmodist.py").read())

### BASE ###
t = pd.DataFrame(BaseAges)
t.columns = ['V1','V2','V3','V4','V5','V6','V7']
ages = np.array(t.V2[np.array(t.V3) == t.V3[1]])
n_age = len(ages)
Nbase = len(t)

######################
### DATA SELECTION ###
######################
data = pd.read_csv(spectra_z_file, sep=',')  # Lista de spectros e redshifts
selection = (data['onoff'] == 1) & (data['extension'] < 100) #Seleciona apenas as lcgs do sdss e as do gemini com abertura = 3 arcsec
data = data[selection]
data.index = range(len(data))

######################
        
### CREATING MAIN OUTPUT FILE ###
Nobj = len(data)
vvv = np.full(Nobj, 'NaN', float) #Cria um vetor onde todos os valores de i=0 até i=Nobj são nan. Vetor a ser preenchido

temp = pd.DataFrame(['lcgID', 'extension', 'ra', 'dec', 'ageL_mean', 'metL_mean', 'aFeL_mean', 'ageM_mean', 'metM_mean', 'aFeM_mean',
          'ageL_gmean', 'metL_gmean', 'aFeL_gmean', 'ageM_gmean', 'metM_gmean', 'aFeM_gmean',
          'chisqrt', 'Mini_log_Mo', 'Mcor_log_Mo', 'magFrac_r', 'v0', 'vd', 'av', 'z', 
          'DM_Mpc', 'DA_Mpc', 'DL_Mpc', 'VC_Gpc3', 'TL_Gyr', 'DM_mag',
          'SN1', 'SN2', 'SN3', 'SN4']).T
temp.to_csv(output_dir+output_file, header=False, index=False, sep=',', mode='w')

### CREATING SFH & CEH OUTPUT FILES ###
temp2 = pd.DataFrame([], columns=[-1,-1,-1, ages[0]/1e9, ages[1]/1e9, ages[2]/1e9, ages[3]/1e9, ages[4]/1e9, ages[5]/1e9, ages[6]/1e9, ages[7]/1e9, ages[8]/1e9, ages[9]/1e9, ages[10]/1e9, ages[11]/1e9, ages[12]/1e9, ages[13]/1e9, ages[14]/1e9])
temp2.to_csv(output_dir+output_file_Mini_SFH, mode='w', sep=' ', index=False)
temp2.to_csv(output_dir+output_file_Mini_cumul_SFH, mode='w', sep=' ', index=False)
temp2.to_csv(output_dir+output_file_Mcor_SFH, mode='w', sep=' ', index=False)
temp2.to_csv(output_dir+output_file_Mcor_cumul_SFH, mode='w', sep=' ', index=False)

temp2.to_csv(output_dir+output_file_Mini_CEH, mode='w', sep=' ', index=False)
temp2.to_csv(output_dir+output_file_Mcor_CEH, mode='w', sep=' ', index=False)
temp2.to_csv(output_dir+output_file_Lum_CEH, mode='w', sep=' ', index=False)

progress = np.around(Nobj * np.linspace(0, 1, 101)).astype(int)

### FILE PICKER ###
for i in range(0,Nobj):
    
    if (data.extension[i] > 10) & (data.extension[i] < 20): #Os nomes do output do starlight são todos <10, mas não quero perder a ordem, mudar output do starlight no futuro
        #for r in range(0,2):
           # if flag_cut == 0:
            #    filename = (starlight_output_dir + '/output_starlight_sdss_LCG' + str(data.lcgID[i]) + '_' + str(data.extension[i]) + '.csv')
        data.extension[i] = data.extension[i] - 10
        flag_sdss = 1
    
    
    if data.flag_sdss[i] == 0: filename = (starlight_output_dir + '/output_starlight_gemini_LCG' + str(data.lcgID[i]) + '_' + str(data.extension[i]) + '.csv')
    if data.flag_sdss[i] == 1: filename = (starlight_output_dir + '/output_starlight_sdss_LCG' + str(data.lcgID[i]) + '_' + str(data.extension[i]) + '.csv')
        
    
    starlight_file = filename
    check = path.exists(starlight_file)
    print(starlight_file)


    if check == True:
        
        ################################################
        # MAIN CODE
        ################################################

        t = pd.read_csv(starlight_file, skiprows = 63, nrows=75, delim_whitespace=True, engine='python', header = None)
        t.columns = ['V1','V2','V3','V4','V5','V6','V7','V8','V9','V10','V11','V12','V13','V14','V15']
        t = t[1:Nbase+1]
        frac = pd.to_numeric(t.V2, errors = 'coerce')
        age = pd.to_numeric(t.V5, errors = 'coerce')
        if(age_format == 'Gyr'):
            age = age * 1e9
        met = pd.to_numeric(t.V6, errors = 'coerce')
        alpha = pd.to_numeric(t.V11, errors = 'coerce')
        mini_frac = pd.to_numeric(t.V3, errors = 'coerce')
        mcor_frac = pd.to_numeric(t.V4, errors = 'coerce')
        frac_star = pd.to_numeric(t.V9, errors = 'coerce')
        ml_ratio = pd.to_numeric(t.V7, errors = 'coerce')

        if met_format == 'logZ_Zsun':
            met = Z_sun * 10**met  
        
        # Main calculator
        # ------------------
        age_L_mean = (sum(frac * age) / sum(frac)) / 1e9
        met_L_mean = np.log10(sum(frac * (met/Z_sun)) / sum(frac))
        alpha_L_mean= np.log10(sum(frac * 10**alpha) / sum(frac))

        age_M_mean = (sum(mcor_frac * age) / sum(mcor_frac)) / 1e9
        met_M_mean = np.log10(sum(mcor_frac * (met/Z_sun)) / sum(mcor_frac))
        alpha_M_mean = np.log10(sum(mcor_frac * 10**alpha) / sum(mcor_frac))

        age_L_gmean = 10**(sum(frac * np.log10(age)) / sum(frac)) / 1e9
        met_L_gmean = sum(frac * np.log10(met/Z_sun)) / sum(frac)
        alpha_L_gmean = sum(frac * alpha) / sum(frac) 

        age_M_gmean = 10**(sum(mcor_frac * np.log10(age)) / sum(mcor_frac)) / 1e9
        met_M_gmean = sum(mcor_frac * np.log10(met/Z_sun)) / sum(mcor_frac)
        alpha_M_gmean = sum(mcor_frac * alpha) / sum(mcor_frac)

        # Get S/N
        # ---------------------------
        t = pd.read_csv(starlight_file, skiprows = 30, nrows=4, header=None, delim_whitespace=True, engine='python')
        
        sn1 = pd.to_numeric(t[0][0], errors = 'coerce')
        sn2 = pd.to_numeric(t[0][1], errors = 'coerce')
        sn3 = pd.to_numeric(t[0][2], errors = 'coerce')
        sn4 = pd.to_numeric(t[0][3], errors = 'coerce')

        # Get vel_disp, Mass
        # ---------------------------
        t = pd.read_csv(starlight_file, skiprows = 49, nrows=12, header=None, sep='                             ', skipinitialspace=True, engine='python', index_col=None)
        
        chisqrt = t[0][0]
        mini = pd.to_numeric(t[0][4], errors = 'coerce')
        mcor = pd.to_numeric(t[0][5], errors = 'coerce')
        v0 = t[0][6]
        vd = t[0][7]
        av = t[0][8]
        vvv[i] = vd 
        
#mcor = sum(mini * mini.frac * 0.01 * frac.star) 
# DH = 4285.714 Mpc = 1.32e+26 m
        
        # Get cosmological distances
        # ---------------------------
        tt = cosmodist(data.z[i], H0 = 70, Omega_m = 0.3, Omega_l = 0.7)

        DH_Mpc = tt.DH_Mpc[0]
        tH_Gyr = tt.tH_Gyr[0]
        DM_Mpc =  tt.DM_Mpc[0]
        DA_Mpc = tt.DA_Mpc[0]
        DL_Mpc = tt.DL_Mpc[0]
        VC_Gpc3 = tt.VC_Gpc3[0]
        TL_Gyr = tt.TL_Gyr[0]
        DM_mag = tt.DM[0]

        dist =  tt.DL_Mpc[0] * 3.086e+24  # Mpc to cm
        
        TL_Gyr = tt.TL_Gyr[0]
        
 
        
        # Get MASS
        # ---------------------------
        Mass_factor = pd.to_numeric(1e-17 * 4 * np.pi * dist**2 / 3.826e33, errors = 'coerce')
        #!!
        if data.flag_sdss[i] == 0: #Gemini
            mass_corr = 10**(-0.4*(data.magSlit_r[i] - data.magPetro_r[i]))
            mini_out = np.log10(mini * Mass_factor / mass_corr)
            mcor_out = np.log10(mcor * Mass_factor / mass_corr)
        if data.flag_sdss[i] == 1: #SDSS
            mass_corr = 10**(-0.4*(data.fiberMag_r[i] - data.magPetro_r[i])) #magFibra_r[i] - magPetro_r[i]
            mini_out = np.log10(mini * Mass_factor / mass_corr)
            mcor_out = np.log10(mcor * Mass_factor / mass_corr)
        
        if flag_sdss == 1:
            data.extension[i] = data.extension[i] + 10
            flag_sdss = 0
            
        ### WRITING MAIN OUTPUT FILES ###
        temp = pd.DataFrame([str(data.lcgID[i]), str(data.extension[i]), str(data.ra[i]), str(data.dec[i]), age_L_mean, met_L_mean, alpha_L_mean, age_M_mean, met_M_mean, alpha_M_mean, 
                  age_L_gmean, met_L_gmean, alpha_L_gmean, age_M_gmean, met_M_gmean, alpha_M_gmean,
                  pd.to_numeric(chisqrt, errors = 'coerce'), mini_out, mcor_out, mass_corr, pd.to_numeric(v0, errors = 'coerce'), pd.to_numeric(vd, errors = 'coerce'), pd.to_numeric(av, errors = 'coerce'), data.z[i], 
                  DM_Mpc, DA_Mpc, DL_Mpc, VC_Gpc3, TL_Gyr, DM_mag,
                  sn1, sn2, sn3, sn4]).T
        temp.to_csv(output_dir+output_file, header=False, index=False, sep=',', mode='a')

        
        ################################################
        # SFH & CEH
        ################################################
        
        # Creating null vectors & matrices
        # ---------------------------------
        met_L_age = pd.DataFrame([np.zeros(n_age)]).T
        met_Mini_age = pd.DataFrame([np.zeros(n_age)]).T
        met_Mcor_age = pd.DataFrame([np.zeros(n_age)]).T
        met_L_age = pd.DataFrame([np.zeros(n_age)]).T
        met_M_age = pd.DataFrame([np.zeros(n_age)]).T
        alpha_L_age = pd.DataFrame([np.zeros(n_age)]).T
        alpha_M_age = pd.DataFrame([np.zeros(n_age)]).T
        frac_age = pd.DataFrame([np.zeros(n_age)]).T
        mini_frac_age = pd.DataFrame([np.zeros(n_age)]).T
        mcor_frac_age = pd.DataFrame([np.zeros(n_age)]).T
        
        # Main calculator
        # ------------------
        for aa in range(n_age):
            xx = abs(age/1e9 - ages[aa]/1e9) < 0.001
            frac_age[0][aa] = sum(frac[xx])
            mini_frac_age[0][aa] = sum(mini_frac[xx])
            mcor_frac_age[0][aa] = sum(mcor_frac[xx])
            if frac_age[0][aa] > 0:
                # sfh
                met_L_age[0][aa] = (np.log10(sum(frac[xx] * met[xx]/0.019) / sum(frac[xx])))
                alpha_L_age[0][aa] = (np.log10(sum(frac[xx] * 10**alpha[xx]) / sum(frac[xx])))
                
                met_M_age[0][aa] = (np.log10(sum(mcor_frac[xx] * met[xx]/0.019) / sum(mcor_frac[xx])))
                alpha_M_age[0][aa] = (np.log10(sum(mcor_frac[xx] * 10**alpha[xx]) / sum(mcor_frac[xx])))
                
                # ceh
                met_L_age[0][aa] = (np.log10(sum(frac[xx] * met[xx]/0.019) / sum(frac[xx])))
                met_Mini_age[0][aa] = (np.log10(sum(mini_frac[xx] * met[xx]/0.019) / sum(mini_frac[xx])))
                met_Mcor_age[0][aa] = (np.log10(sum(mcor_frac[xx] * met[xx]/0.019) / sum(mcor_frac[xx])))
                
        mcor_frac_age = mcor_frac_age/sum(mcor_frac_age[0])
        mini_frac_age = mini_frac_age/sum(mini_frac_age[0])
        
        # Cumulative calculator
        # ---------------------------
        mini_cumul_age = pd.DataFrame([np.zeros(n_age)]).T
        mcor_cumul_age = pd.DataFrame([np.zeros(n_age)]).T
        for aa in range(n_age):
            mini_cumul_age[0][aa] = sum(mini_frac_age[0][ages >= ages[aa]])
            mcor_cumul_age[0][aa] = sum(mcor_frac_age[0][ages >= ages[aa]])

            
        ### WRITING OUTPUT FILES ###
        
        # sfh transpose
        mini_frac_age = mini_frac_age.T
        mini_cumul_age = mini_cumul_age.T
        mcor_frac_age = mcor_frac_age.T
        mcor_cumul_age = mcor_cumul_age.T
        
        # ceh transpose
        met_L_age = met_L_age.T
        met_Mini_age = met_Mini_age.T
        met_Mcor_age = met_Mcor_age.T
        
        # sfh writing
        mini_frac_age.to_csv(output_dir+output_file_Mini_SFH, mode='a', sep=' ', index=False, header=False)
        mini_cumul_age.to_csv(output_dir+output_file_Mini_cumul_SFH, mode='a', sep=' ', index=False, header=False)
        mcor_frac_age.to_csv(output_dir+output_file_Mcor_SFH, mode='a', sep=' ', index=False, header=False)
        mcor_cumul_age.to_csv(output_dir+output_file_Mcor_cumul_SFH, mode='a', sep=' ', index=False, header=False)
        
        # ceh writing
        met_L_age.to_csv(output_dir+output_file_Lum_CEH, mode='a', sep=' ', index=False, header=False)
        met_Mini_age.to_csv(output_dir+output_file_Mini_CEH, mode='a', sep=' ', index=False, header=False)
        met_Mcor_age.to_csv(output_dir+output_file_Mcor_CEH, mode='a', sep=' ', index=False, header=False)
        
        ### PROGRESS ###
        if i in progress:
            print('i=' + str(i+1) + ' ' + 'Nobj=' + str(Nobj) + ' ' + 'i/Nobj=' + str(round(((i+1)/Nobj) * 100,2)) + '%')
            localtime = time.asctime( time.localtime(time.time()) )
            print(localtime)
            
