import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
pd.set_option('max_columns', None)
pd.set_option('max_rows', None)

sdss = pd.read_csv('/home/vitorbootz/research/flux_measurements/sdss_flux_lines.csv')
sample = pd.read_csv('/home/vitorbootz/research/flux_measurements/sample_flux_lines.csv')
main = pd.read_csv('/home/vitorbootz/research/aux_files/galaxy_list.csv')

sdss.sort_values(['lcgID'], inplace=True)
sdss.index = range(len(sdss))
sdss_sample = pd.DataFrame([])
main_sample = pd.DataFrame([])
t = pd.DataFrame([])
for i in range(len(sdss)): sdss.ra[i] = round(sdss.ra[i],3)
for i in range(len(sdss)): sdss.dec[i] = round(sdss.dec[i],3)
    
for i in range(len(sample)): 
    sdss_sample = sdss_sample.append(sdss[(sdss.ra == round(sample.ra[i],3)) & (sdss.dec == round(sample.dec[i],3))], ignore_index=True)

sdss_sample = sdss_sample.drop_duplicates(['ra'])
sdss_sample.index = range(len(sdss_sample))
for i in range(len(sample)): main_sample = main_sample.append(main[(sample.lcgID[i] == main.lcgID) & (sample.extension[i] == main.extension)])
main_sample.index = range(len(main_sample)) 

main_sdss = main_sample[main_sample.flag_sdss == 1]
main_g = main_sample[(main_sample.flag_sdss == 0) & (main_sample.flag_lcg) == 1]
main_sdss.index = range(len(main_sdss)) 
main_g.index = range(len(main_g)) 

for i in range(len(sample)):
    for j in range(len(sdss_sample)):
        if round(sample.ra[i],3) == sdss_sample.ra[j]:
            df = pd.DataFrame([sample.lcgID[i].astype(int), sample.extension[i].astype(int), sample.ra[i], sample.dec[i],
                  sample.OII1[i],sample.OII2[i],sample.Hb[i],sample.OIII_4959[i],sample.OIII_5008[i],sample.Ha[i],sample.NII_6583[i],sample.NII_6548[i], sample.SII_6718[i], sample.SII_6732[i],
                  sdss_sample.oii_3726_flux[j], sdss_sample.oii_3729_flux[j], sdss_sample.h_beta_flux[j],sdss_sample.oiii_4959_flux[j],sdss_sample.oiii_5007_flux[j], sdss_sample.nii_6548_flux[j], sdss_sample.h_alpha_flux[j], sdss_sample.nii_6584_flux[j],sdss_sample.sii_6717_flux[j],sdss_sample.sii_6731_flux[j]]).T
            t = t.append(df)
            
t.columns = ['lcgID', 'extension', 'ra', 'dec', 'OII1','OII2','Hb','OIII_4959','OIII_5007','Ha','NII_6583','NII_6548', 'SII_6718', 'SII_6732','OII1_SDSS','OII2_SDSS','Hb_SDSS','OIII_4959_SDSS','OIII_5007_SDSS','NII_6548_SDSS', 'Ha_SDSS', 'NII_6583_SDSS','SII_6718_SDSS','SII_6732_SDSS']
t.index = range(len(t))

for i in range(len(main)):
    for j in range(len(t)):
        if (main.lcgID[i] == t.lcgID[j]) & (main.extension[i] == t.extension[j]):
            if main.flag_sdss[i] == 0: #Gemini
                # No primeiro bloco, os valores medidos dos fluxos oriundos do Gemini
                mass_corr = 10**(-0.4*(main.magSlit_r[i] - main.magPetro_r[i]))
                t.OII1[j] = t.OII1[j] / mass_corr
                t.OII2[j] = t.OII2[j] / mass_corr
                t.Hb[j] = t.Hb[j] / mass_corr
                t.OIII_4959[j] = t.OIII_4959[j] / mass_corr
                t.OIII_5007[j] = t.OIII_5007[j] / mass_corr
                t.Ha[j] = t.Ha[j] / mass_corr
                t.NII_6583[j] = t.NII_6583[j] / mass_corr
                t.NII_6548[j] = t.NII_6548[j] / mass_corr
                t.SII_6718[j] = t.SII_6718[j] / mass_corr
                t.SII_6732[j] = t.SII_6732[j] / mass_corr
                
                # Neste segundo, os fluxos da tabela do sloan
                mass_corr = 10**(-0.4*(main.fiberMag_r[i] - main.magPetro_r[i]))
                t.OII1_SDSS[j] = t.OII1_SDSS[j] / mass_corr
                t.OII2_SDSS[j] = t.OII2_SDSS[j] / mass_corr
                t.Hb_SDSS[j] = t.Hb_SDSS[j] / mass_corr
                t.OIII_4959_SDSS[j] = t.OIII_4959_SDSS[j] / mass_corr
                t.OIII_5007_SDSS[j] = t.OIII_5007_SDSS[j] / mass_corr
                t.Ha_SDSS[j] = t.Ha_SDSS[j] / mass_corr
                t.NII_6583_SDSS[j] = t.NII_6583_SDSS[j] / mass_corr
                t.NII_6548_SDSS[j] = t.NII_6548_SDSS[j] / mass_corr
                t.SII_6718_SDSS[j] = t.SII_6718_SDSS[j] / mass_corr
                t.SII_6732_SDSS[j] = t.SII_6732_SDSS[j] / mass_corr
                
                
            if main.flag_sdss[i] == 1: #SDSS
                # Neste caso, todos são do sloan, portanto todos estão associados à Fibra apenas
                mass_corr = 10**(-0.4*(main.fiberMag_r[i] - main.magPetro_r[i]))
                t.OII1[j] = t.OII1[j] / mass_corr
                t.OII2[j] = t.OII2[j] / mass_corr
                t.Hb[j] = t.Hb[j] / mass_corr
                t.OIII_4959[j] = t.OIII_4959[j] / mass_corr
                t.OIII_5007[j] = t.OIII_5007[j] / mass_corr
                t.Ha[j] = t.Ha[j] / mass_corr
                t.NII_6583[j] = t.NII_6583[j] / mass_corr
                t.NII_6548[j] = t.NII_6548[j] / mass_corr
                t.SII_6718[j] = t.SII_6718[j] / mass_corr
                t.SII_6732[j] = t.SII_6732[j] / mass_corr
                
                t.OII1_SDSS[j] = t.OII1_SDSS[j] / mass_corr
                t.OII2_SDSS[j] = t.OII2_SDSS[j] / mass_corr
                t.Hb_SDSS[j] = t.Hb_SDSS[j] / mass_corr
                t.OIII_4959_SDSS[j] = t.OIII_4959_SDSS[j] / mass_corr
                t.OIII_5007_SDSS[j] = t.OIII_5007_SDSS[j] / mass_corr
                t.Ha_SDSS[j] = t.Ha_SDSS[j] / mass_corr
                t.NII_6583_SDSS[j] = t.NII_6583_SDSS[j] / mass_corr
                t.NII_6548_SDSS[j] = t.NII_6548_SDSS[j] / mass_corr
                t.SII_6718_SDSS[j] = t.SII_6718_SDSS[j] / mass_corr
                t.SII_6732_SDSS[j] = t.SII_6732_SDSS[j] / mass_corr
                

# Jogados fora: grupo 2361, 2023_1 e 3090_2
useful = t[(t.lcgID != 2361) & (t.extension != 13) & (t.ra != 131.36510) & (t.ra != 183.40380) & (t.ra != 164.07920)]

#####################################
#####################################

mosaic = plt.figure(figsize=(14, 22))
gs = gridspec.GridSpec(nrows=4, ncols=2, figure=mosaic)

mosaic_bpt = plt.figure(figsize=(14, 22))
gs_bpt = gridspec.GridSpec(nrows=1, ncols=1, figure=mosaic_bpt)

mosaic_OH = plt.figure(figsize=(14, 22))
gs_OH = gridspec.GridSpec(nrows=1, ncols=3, figure=mosaic_OH)

#####################################
############ Hb Line ################

fig1 = plt.figure(figsize=(6,4))
axis1 = mosaic.add_subplot(gs[0, 1])
axis1.grid(alpha=0.2, color='grey')

razao_Hb = np.log10(useful.Hb/useful.Hb_SDSS)
median_Hb = razao_Hb.median()
std_Hb = razao_Hb.std()

axis1.plot(np.log10(useful.Hb_SDSS),np.log10(useful.Hb),'x', color='black', ms='9')
axis1.set_xlabel('log(H'+r'$\beta_{SDSS}$' ' / ' '$[10^{-17} erg.s^{-1}.cm^{-2}])$')
axis1.set_ylabel('log(H'+r'$\beta_{Medido}$' ' / ' '$[10^{-17} erg.s^{-1}.cm^{-2}])$')

mini = min(np.log10(useful.Hb_SDSS)-2)
maxi = max(np.log10(useful.Hb)+2)
axis1.set_xlim(mini+1.9,maxi-1.9)
axis1.set_ylim(mini+1.9,maxi-1.9)

axis1.fill_between([mini,maxi], [mini-0.1,maxi-0.1], [mini+0.1,maxi+0.1], 
                   where=([mini-0.1,maxi-0.1] <= [mini+0.1,maxi+0.1]), color='gray', alpha=0.4, interpolate=True, label=r'$\pm 0.1 dex$')
axis1.plot([mini-0.1,maxi-0.1],[mini,maxi], color='black', alpha=0.5, ls='--')
axis1.plot([mini,maxi],[mini,maxi], color='black', alpha=0.5, ls='-')
axis1.plot([mini+0.1,maxi+0.1],[mini,maxi], color='black', alpha=0.5, ls='--')
axis1.legend()

axis1.set_title('Comparação dos fluxos de H'+r'$\beta$')

axis1.text(0.75,0.05,r'log$_{10}\left(\frac{F_{Medido}}{F_{SDSS}}\right) \Rightarrow$', fontsize=15, transform=axis1.transAxes, horizontalalignment='right', verticalalignment='bottom')
axis1.text(0.96,0.12,'Mediana: '+str(round(median_Hb,2)), fontsize=10, transform=axis1.transAxes, horizontalalignment='right', verticalalignment='bottom')
axis1.text(0.95,0.05,'Desvio: '+str(round(std_Hb,2)), fontsize=10, transform=axis1.transAxes, horizontalalignment='right', verticalalignment='bottom')

#fig1.savefig('/home/vitorbootz/research/TCC_images/flux_comparsion/flux_comparsion_Hb.pdf')

##########################################
############ OIII 5007 Line ##############

fig2 = plt.figure(figsize=(6,4))
axis2 = mosaic.add_subplot(gs[1, 1])
axis2.grid(alpha=0.2, color='grey')

razao_OIII_5007 = np.log10(useful.OIII_5007/useful.OIII_5007_SDSS)
median_OIII_5007 = razao_OIII_5007.median()
std_OIII_5007 = razao_OIII_5007.std()

axis2.plot(np.log10(useful.OIII_5007_SDSS),np.log10(useful.OIII_5007),'x', color='black', ms='9')
axis2.set_xlabel('log(OIII'+r'$_{SDSS}^{5007}$' ' / ' '$[10^{-17} erg.s^{-1}.cm^{-2}])$')
axis2.set_ylabel('log(OIII'+r'$_{Medido}^{5007}$' ' / ' '$[10^{-17} erg.s^{-1}.cm^{-2}])$')

mini = min(np.log10(useful.OIII_5007_SDSS)-2)
maxi = max(np.log10(useful.OIII_5007)+2)

axis2.set_xlim(mini+1.9,maxi-1.9)
axis2.set_ylim(mini+1.9,maxi-1.9)

axis2.fill_between([mini,maxi], [mini-0.1,maxi-0.1], [mini+0.1,maxi+0.1], 
                   where=([mini-0.1,maxi-0.1] <= [mini+0.1,maxi+0.1]), color='gray', alpha=0.4, interpolate=True, label=r'$\pm 0.1 dex$')
axis2.plot([mini-0.1,maxi-0.1],[mini,maxi], color='black', alpha=0.5, ls='--')
axis2.plot([mini,maxi],[mini,maxi], color='black', alpha=0.5, ls='-')
axis2.plot([mini+0.1,maxi+0.1],[mini,maxi], color='black', alpha=0.5, ls='--')
axis2.legend()

axis2.set_title('Comparação dos fluxos de [OIII]'+'$\lambda 5007\ \AA$')

axis2.text(0.75,0.05,r'log$_{10}\left(\frac{F_{Medido}}{F_{SDSS}}\right) \Rightarrow$', fontsize=15, transform=axis2.transAxes, horizontalalignment='right', verticalalignment='bottom')
axis2.text(0.96,0.12,'Mediana: '+str(round(median_OIII_5007,2)), fontsize=10, transform=axis2.transAxes, horizontalalignment='right', verticalalignment='bottom')
axis2.text(0.95,0.05,'Desvio: '+str(round(std_OIII_5007,2)), fontsize=10, transform=axis2.transAxes, horizontalalignment='right', verticalalignment='bottom')

#fig2.savefig('/home/vitorbootz/research/TCC_images/flux_comparsion/flux_comparsion_OIII_5007.pdf')

#####################################
############ Ha Line ################

fig3 = plt.figure(figsize=(6,4))
axis3 = mosaic.add_subplot(gs[0, 0])
axis3.grid(alpha=0.2, color='grey')

razao_Ha = np.log10(useful.Ha/useful.Ha_SDSS)
median_Ha = razao_Ha.median()
std_Ha = razao_Ha.std()

axis3.plot(np.log10(useful.Ha_SDSS),np.log10(useful.Ha), 'x', color='black', ms='9')
axis3.set_xlabel('log(H'+r'$\alpha_{SDSS}$' ' / ' '$[10^{-17} erg.s^{-1}.cm^{-2}])$')
axis3.set_ylabel('log(H'+r'$\alpha_{Medido}$' ' / ' '$[10^{-17} erg.s^{-1}.cm^{-2}])$')

mini = min(np.log10(useful.Ha_SDSS)-2)
maxi = max(np.log10(useful.Ha)+2)

axis3.set_xlim(mini+1.9,maxi-1.9)
axis3.set_ylim(mini+1.9,maxi-1.9)

axis3.fill_between([mini,maxi], [mini-0.1,maxi-0.1], [mini+0.1,maxi+0.1], 
                   where=([mini-0.1,maxi-0.1] <= [mini+0.1,maxi+0.1]), color='gray', alpha=0.4, interpolate=True, label=r'$\pm 0.1 dex$')
axis3.plot([mini-0.1,maxi-0.1],[mini,maxi], color='black', alpha=0.5, ls='--')
axis3.plot([mini,maxi],[mini,maxi], color='black', alpha=0.5, ls='-')
axis3.plot([mini+0.1,maxi+0.1],[mini,maxi], color='black', alpha=0.5, ls='--')
axis3.legend()

axis3.set_title('Comparação dos fluxos de H'+r'$\alpha$')

axis3.text(0.75,0.05,r'log$_{10}\left(\frac{F_{Medido}}{F_{SDSS}}\right) \Rightarrow$', fontsize=15, transform=axis3.transAxes, horizontalalignment='right', verticalalignment='bottom')
axis3.text(0.96,0.12,'Mediana: '+str(round(median_Ha,2)), fontsize=10, transform=axis3.transAxes, horizontalalignment='right', verticalalignment='bottom')
axis3.text(0.96,0.05,'Desvio: '+str(round(std_Ha,2)), fontsize=10, transform=axis3.transAxes, horizontalalignment='right', verticalalignment='bottom')

#fig3.savefig('/home/vitorbootz/research/TCC_images/flux_comparsion/flux_comparsion_Ha.pdf')

##########################################
############ NII_6583 Line ###############

fig4 = plt.figure(figsize=(6,4))
axis4 = mosaic.add_subplot(gs[2, 1])
axis4.grid(alpha=0.2, color='grey')

razao_NII_6583 = np.log10(useful.NII_6583/useful.NII_6583_SDSS)
median_NII_6583 = razao_NII_6583.median()
std_NII_6583 = razao_NII_6583.std()

axis4.plot(np.log10(useful.NII_6583_SDSS),np.log10(useful.NII_6583),'x', color='black', ms='9')
axis4.set_xlabel('log(NII'+r'$^{6583}_{SDSS}$' ' / ' '$[10^{-17} erg.s^{-1}.cm^{-2}])$')
axis4.set_ylabel('log(NII'+r'$^{6583}_{Medido}$' ' / ' '$[10^{-17} erg.s^{-1}.cm^{-2}])$')

mini = min(np.log10(useful.NII_6583_SDSS)-2)
maxi = max(np.log10(useful.NII_6583)+2)

axis4.set_xlim(mini+1.9,maxi-1.9)
axis4.set_ylim(mini+1.9,maxi-1.9)

axis4.fill_between([mini,maxi], [mini-0.1,maxi-0.1], [mini+0.1,maxi+0.1], 
                   where=([mini-0.1,maxi-0.1] <= [mini+0.1,maxi+0.1]), color='gray', alpha=0.4, interpolate=True, label=r'$\pm 0.1 dex$')
axis4.plot([mini-0.1,maxi-0.1],[mini,maxi], color='black', alpha=0.5, ls='--')
axis4.plot([mini,maxi],[mini,maxi], color='black', alpha=0.5, ls='-')
axis4.plot([mini+0.1,maxi+0.1],[mini,maxi], color='black', alpha=0.5, ls='--')
axis4.legend()

axis4.set_title('Comparação dos fluxos de [NII]'+'$\lambda 6584\ \AA$')

axis4.text(0.75,0.05,r'log$_{10}\left(\frac{F_{Medido}}{F_{SDSS}}\right) \Rightarrow$', fontsize=15, transform=axis4.transAxes, horizontalalignment='right', verticalalignment='bottom')
axis4.text(0.96,0.12,'Mediana: '+str(round(median_NII_6583,2)), fontsize=10, transform=axis4.transAxes, horizontalalignment='right', verticalalignment='bottom')
axis4.text(0.95,0.05,'Desvio: '+str(round(std_NII_6583,2)), fontsize=10, transform=axis4.transAxes, horizontalalignment='right', verticalalignment='bottom')
#fig4.savefig('/home/vitorbootz/research/TCC_images/flux_comparsion/flux_comparsion_NII_6583.pdf')

#########################################
############ NII/Ha ratio ###############

fig5 = plt.figure(figsize=(6,4))
axis5 = mosaic_bpt.add_subplot(gs[0, 0])
axis5.grid(alpha=0.2, color='grey')

razao_NII_Ha = np.log10(useful.NII_6583/useful.Ha) - np.log10(useful.NII_6583_SDSS/useful.Ha_SDSS)
median_NII_Ha = razao_NII_Ha.median()
std_NII_Ha = razao_NII_Ha.std()

axis5.plot(np.log10(useful.NII_6583_SDSS/useful.Ha_SDSS),np.log10(useful.NII_6583/useful.Ha),'x', color='black', ms='9')
axis5.set_xlabel('log([NII/H'+r'$\alpha]_{SDSS}$' ' / ' '$[10^{-17} erg.s^{-1}.cm^{-2}])$')
axis5.set_ylabel('log([NII/H'+r'$\alpha]_{Medido}$' ' / ' '$[10^{-17} erg.s^{-1}.cm^{-2}])$')

mini = min(np.log10(useful.NII_6583_SDSS/useful.Ha_SDSS)-2)
maxi = max(np.log10(useful.NII_6583/useful.Ha)+2)

axis5.set_xlim(mini+1.9,maxi-1.9)
axis5.set_ylim(mini+1.9,maxi-1.9)

axis5.fill_between([mini,maxi], [mini-0.1,maxi-0.1], [mini+0.1,maxi+0.1], 
                   where=([mini-0.1,maxi-0.1] <= [mini+0.1,maxi+0.1]), color='gray', alpha=0.4, interpolate=True, label=r'$\pm 0.1 dex$')
axis5.plot([mini-0.1,maxi-0.1],[mini,maxi], color='black', alpha=0.5, ls='--')
axis5.plot([mini,maxi],[mini,maxi], color='black', alpha=0.5, ls='-')
axis5.plot([mini+0.1,maxi+0.1],[mini,maxi], color='black', alpha=0.5, ls='--')
axis5.legend()

axis5.set_title('Comparação da razão [NII/Ha]')

axis5.text(0.75,0.05,r'log$_{10}\left(\frac{F^{NII/Ha}_{Medido}}{F^{NII/Ha}_{SDSS}}\right) \Rightarrow$', fontsize=15, transform=axis5.transAxes, horizontalalignment='right', verticalalignment='bottom')
axis5.text(0.96,0.13,'Mediana: '+str(round(median_NII_Ha,2)), fontsize=10, transform=axis5.transAxes, horizontalalignment='right', verticalalignment='bottom')
axis5.text(0.95,0.07,'Desvio: '+str(round(std_NII_Ha,2)), fontsize=10, transform=axis5.transAxes, horizontalalignment='right', verticalalignment='bottom')
#fig5.savefig('/home/vitorbootz/research/TCC_images/flux_comparsion/flux_comparsion_NII_Ha.pdf')

#########################################
############ OIII/Hb ratio ##############


fig6 = plt.figure(figsize=(6,4))
axis6 = mosaic_bpt.add_subplot(gs[0, 1])
axis6.grid(alpha=0.2, color='grey')

razao_OIII_Hb = np.log10(useful.OIII_5007/useful.Hb) - np.log10(useful.OIII_5007_SDSS/useful.Hb_SDSS)
median_OIII_Hb = razao_OIII_Hb.median()
std_OIII_Hb = razao_OIII_Hb.std()

axis6.plot(np.log10(useful.OIII_5007_SDSS/useful.Hb_SDSS),np.log10(useful.OIII_5007/useful.Hb),'x', color='black', ms='9')
axis6.set_xlabel('log([OIII/H'+r'$\beta]_{SDSS}$' ' / ' '$[10^{-17} erg.s^{-1}.cm^{-2}])$')
axis6.set_ylabel('log([OIII/H'+r'$\beta]_{Medido}$' ' / ' '$[10^{-17} erg.s^{-1}.cm^{-2}])$')

mini = min(np.log10(useful.OIII_5007_SDSS/useful.Hb_SDSS)-2)
maxi = max(np.log10(useful.OIII_5007/useful.Hb)+2)

axis6.set_xlim(mini+1.9,maxi-1.9)
axis6.set_ylim(mini+1.9,maxi-1.9)

axis6.fill_between([mini,maxi], [mini-0.1,maxi-0.1], [mini+0.1,maxi+0.1], 
                   where=([mini-0.1,maxi-0.1] <= [mini+0.1,maxi+0.1]), color='gray', alpha=0.4, interpolate=True, label=r'$\pm 0.1 dex$')
axis6.plot([mini-0.1,maxi-0.1],[mini,maxi], color='black', alpha=0.5, ls='--')
axis6.plot([mini,maxi],[mini,maxi], color='black', alpha=0.5, ls='-')
axis6.plot([mini+0.1,maxi+0.1],[mini,maxi], color='black', alpha=0.5, ls='--')
axis6.legend()

axis6.set_title('Comparação da razão [OIII/Hb]')

axis6.text(0.75,0.05,r'log$_{10}\left(\frac{F^{OIII/H\beta}_{Medido}}{F^{OIII/H\beta}_{SDSS}}\right) \Rightarrow$', fontsize=15, transform=axis6.transAxes, horizontalalignment='right', verticalalignment='bottom')
axis6.text(0.96,0.13,'Mediana: '+str(round(median_OIII_Hb,2)), fontsize=10, transform=axis6.transAxes, horizontalalignment='right', verticalalignment='bottom')
axis6.text(0.95,0.07,'Desvio: '+str(round(std_OIII_Hb,2)), fontsize=10, transform=axis6.transAxes, horizontalalignment='right', verticalalignment='bottom')
#fig6.savefig('/home/vitorbootz/research/TCC_images/flux_comparsion/flux_comparsion_OIII_Hb.pdf')

#####################################
############ SII_6718 Line ###############


fig7 = plt.figure(figsize=(6,4))
axis7 = mosaic.add_subplot(gs[3, 0])
axis7.grid(alpha=0.2, color='grey')

razao_SII_6583 = np.log10(useful.SII_6718/useful.SII_6718_SDSS)
median_SII_6583 = razao_SII_6583.median()
std_SII_6583 = razao_SII_6583.std()

axis7.plot(np.log10(useful.SII_6718_SDSS),np.log10(useful.SII_6718),'x', color='black', ms='9')
axis7.set_xlabel('log(SII'+r'$^{6718}_{SDSS}$' ' / ' '$[10^{-17} erg.s^{-1}.cm^{-2}])$')
axis7.set_ylabel('log(SII'+r'$^{6718}_{Medido}$' ' / ' '$[10^{-17} erg.s^{-1}.cm^{-2}])$')

mini = min(np.log10(useful.SII_6718_SDSS)-2)
maxi = max(np.log10(useful.SII_6718)+2)

axis7.set_xlim(mini+1.9,maxi-1.9)
axis7.set_ylim(mini+1.9,maxi-1.9)

axis7.fill_between([mini,maxi], [mini-0.1,maxi-0.1], [mini+0.1,maxi+0.1], 
                   where=([mini-0.1,maxi-0.1] <= [mini+0.1,maxi+0.1]), color='gray', alpha=0.4, interpolate=True, label=r'$\pm 0.1 dex$')
axis7.plot([mini-0.1,maxi-0.1],[mini,maxi], color='black', alpha=0.5, ls='--')
axis7.plot([mini,maxi],[mini,maxi], color='black', alpha=0.5, ls='-')
axis7.plot([mini+0.1,maxi+0.1],[mini,maxi], color='black', alpha=0.5, ls='--')
axis7.legend()

axis7.set_title('Comparação dos fluxos de [SII]'+'$\lambda 6718\ \AA$')

axis7.text(0.75,0.08,r'log$_{10}\left(\frac{F_{Medido}}{F_{SDSS}}\right) \Rightarrow$', fontsize=15, transform=axis7.transAxes, horizontalalignment='right', verticalalignment='bottom')
axis7.text(0.96,0.13,'Mediana: '+str(round(median_SII_6583,2)), fontsize=10, transform=axis7.transAxes, horizontalalignment='right', verticalalignment='bottom')
axis7.text(0.95,0.07,'Desvio: '+str(round(std_SII_6583,2)), fontsize=10, transform=axis7.transAxes, horizontalalignment='right', verticalalignment='bottom')
#fig7.savefig('/home/vitorbootz/research/TCC_images/flux_comparsion/flux_comparsion_SII_6718.pdf')

#####################################
############ SII_6732 Line ###############

fig8 = plt.figure(figsize=(6,4))
axis8 = mosaic.add_subplot(gs[3, 1])
axis8.grid(alpha=0.2, color='grey')

razao_SII_6732 = np.log10(useful.SII_6732/useful.SII_6732_SDSS)
median_SII_6732 = razao_SII_6732.median()
std_SII_6732 = razao_SII_6732.std()

axis8.plot(np.log10(useful.SII_6732_SDSS),np.log10(useful.SII_6732),'x', color='black', ms='9')
axis8.set_xlabel('log(SII'+r'$^{6732}_{SDSS}$' ' / ' '$[10^{-17} erg.s^{-1}.cm^{-2}])$')
axis8.set_ylabel('log(SII'+r'$^{6732}_{Medido}$' ' / ' '$[10^{-17} erg.s^{-1}.cm^{-2}])$')

mini = min(np.log10(useful.SII_6732_SDSS)-2)
maxi = max(np.log10(useful.SII_6732)+2)

axis8.set_xlim(mini+1.9,maxi-1.9)
axis8.set_ylim(mini+1.9,maxi-1.9)

axis8.fill_between([mini,maxi], [mini-0.1,maxi-0.1], [mini+0.1,maxi+0.1], 
                   where=([mini-0.1,maxi-0.1] <= [mini+0.1,maxi+0.1]), color='gray', alpha=0.4, interpolate=True, label=r'$\pm 0.1 dex$')
axis8.plot([mini-0.1,maxi-0.1],[mini,maxi], color='black', alpha=0.5, ls='--')
axis8.plot([mini,maxi],[mini,maxi], color='black', alpha=0.5, ls='-')
axis8.plot([mini+0.1,maxi+0.1],[mini,maxi], color='black', alpha=0.5, ls='--')
axis8.legend()

axis8.set_title('Comparação dos fluxos de [SII]'+'$\lambda 6732\ \AA$')

axis8.text(0.75,0.06,r'log$_{10}\left(\frac{F_{Medido}}{F_{SDSS}}\right) \Rightarrow$', fontsize=15, transform=axis8.transAxes, horizontalalignment='right', verticalalignment='bottom')
axis8.text(0.96,0.13,'Mediana: '+str(round(median_SII_6732,2)), fontsize=10, transform=axis8.transAxes, horizontalalignment='right', verticalalignment='bottom')
axis8.text(0.94,0.07,'Desvio: '+str(round(std_SII_6732,2)), fontsize=10, transform=axis8.transAxes, horizontalalignment='right', verticalalignment='bottom')
#fig8.savefig('/home/vitorbootz/research/TCC_images/flux_comparsion/flux_comparsion_SII_6732.pdf')

##########################################
############ NII_6548 Line ###############

fig9 = plt.figure(figsize=(6,4))
axis9 = mosaic.add_subplot(gs[2, 0])
axis9.grid(alpha=0.2, color='grey')

razao_NII_6548 = np.log10(useful.NII_6548/useful.NII_6548_SDSS)
median_NII_6548 = razao_NII_6548.median()
std_NII_6548 = razao_NII_6548.std()

axis9.plot(np.log10(useful.NII_6548_SDSS),np.log10(useful.NII_6548),'x', color='black', ms='9')
axis9.set_xlabel('log(NII'+r'$^{6548}_{SDSS}$' ' / ' '$[10^{-17} erg.s^{-1}.cm^{-2}])$')
axis9.set_ylabel('log(NII'+r'$^{6548}_{Medido}$' ' / ' '$[10^{-17} erg.s^{-1}.cm^{-2}])$')

mini = min(np.log10(useful.NII_6548_SDSS)-2)
maxi = max(np.log10(useful.NII_6548)+2)

axis9.set_xlim(mini+1.9,maxi-1.9)
axis9.set_ylim(mini+1.9,maxi-1.9)

axis9.fill_between([mini,maxi], [mini-0.1,maxi-0.1], [mini+0.1,maxi+0.1], 
                   where=([mini-0.1,maxi-0.1] <= [mini+0.1,maxi+0.1]), color='gray', alpha=0.4, interpolate=True, label=r'$\pm 0.1 dex$')
axis9.plot([mini-0.1,maxi-0.1],[mini,maxi], color='black', alpha=0.5, ls='--')
axis9.plot([mini,maxi],[mini,maxi], color='black', alpha=0.5, ls='-')
axis9.plot([mini+0.1,maxi+0.1],[mini,maxi], color='black', alpha=0.5, ls='--')
axis9.legend()

axis9.set_title('Comparação dos fluxos de [NII]'+'$\lambda 6548\ \AA$')

axis9.text(0.75,0.06,r'log$_{10}\left(\frac{F_{Medido}}{F_{SDSS}}\right) \Rightarrow$', fontsize=15, transform=axis9.transAxes, horizontalalignment='right', verticalalignment='bottom')
axis9.text(0.96,0.13,'Mediana: '+str(round(median_NII_6548,2)), fontsize=10, transform=axis9.transAxes, horizontalalignment='right', verticalalignment='bottom')
axis9.text(0.95,0.07,'Desvio: '+str(round(std_NII_6548,2)), fontsize=10, transform=axis9.transAxes, horizontalalignment='right', verticalalignment='bottom')
#fig9.savefig('/home/vitorbootz/research/TCC_images/flux_comparsion/flux_comparsion_NII_6548.pdf')

##########################################
############ OIII 4959 Line ##############

fig10 = plt.figure(figsize=(6,4))
axis10 = mosaic.add_subplot(gs[1, 0])
axis10.grid(alpha=0.2, color='grey')

razao_OIII_4959 = np.log10(useful.OIII_4959/useful.OIII_4959_SDSS)
median_OIII_4959 = razao_OIII_4959.median()
std_OIII_4959 = razao_OIII_4959.std()

axis10.plot(np.log10(useful.OIII_4959_SDSS),np.log10(useful.OIII_4959),'x', color='black', ms='9')
axis10.set_xlabel('log(OIII'+r'$_{SDSS}^{4959}$' ' / ' '$[10^{-17} erg.s^{-1}.cm^{-2}])$')
axis10.set_ylabel('log(OIII'+r'$_{Medido}^{4959}$' ' / ' '$[10^{-17} erg.s^{-1}.cm^{-2}])$')

mini = min(np.log10(useful.OIII_4959_SDSS)-2)
maxi = max(np.log10(useful.OIII_4959)+2)

axis10.set_xlim(mini+1.9,maxi-1.9)
axis10.set_ylim(mini+1.9,maxi-1.9)

axis10.fill_between([mini,maxi], [mini-0.1,maxi-0.1], [mini+0.1,maxi+0.1], 
                   where=([mini-0.1,maxi-0.1] <= [mini+0.1,maxi+0.1]), color='gray', alpha=0.4, interpolate=True, label=r'$\pm 0.1 dex$')
axis10.plot([mini-0.1,maxi-0.1],[mini,maxi], color='black', alpha=0.5, ls='--')
axis10.plot([mini,maxi],[mini,maxi], color='black', alpha=0.5, ls='-')
axis10.plot([mini+0.1,maxi+0.1],[mini,maxi], color='black', alpha=0.5, ls='--')
axis10.legend()

axis10.set_title('Comparação dos fluxos de [OIII]'+'$\lambda 4959\ \AA$')

axis10.text(0.75,0.05,r'log$_{10}\left(\frac{F_{Medido}}{F_{SDSS}}\right) \Rightarrow$', fontsize=15, transform=axis10.transAxes, horizontalalignment='right', verticalalignment='bottom')
axis10.text(0.96,0.12,'Mediana: '+str(round(median_OIII_4959,2)), fontsize=10, transform=axis10.transAxes, horizontalalignment='right', verticalalignment='bottom')
axis10.text(0.95,0.05,'Desvio: '+str(round(std_OIII_4959,2)), fontsize=10, transform=axis10.transAxes, horizontalalignment='right', verticalalignment='bottom')

#fig10.savefig('/home/vitorbootz/research/TCC_images/flux_comparsion/flux_comparsion_OIII_4959.pdf')

#########################################
############ (OIII1+OIII2)/Hb ratio ##############


fig11 = plt.figure(figsize=(6,4))
axis11 = mosaic_OH.add_subplot(gs[0, 0])
axis11.grid(alpha=0.2, color='grey')

razao_OIII_Hb = np.log10((useful.OIII_5007+useful.OIII_4959)/useful.Hb) - np.log10((useful.OIII_5007_SDSS+useful.OIII_4959_SDSS)/useful.Hb_SDSS)
median_OIII_Hb = razao_OIII_Hb.median()
std_OIII_Hb = razao_OIII_Hb.std()

axis11.plot(np.log10((useful.OIII_5007_SDSS+useful.OIII_4959_SDSS)/useful.Hb_SDSS),np.log10((useful.OIII_5007+useful.OIII_4959)/useful.Hb),'x', color='black', ms='9')
axis11.set_xlabel(r'log([R3]$_{SDSS}$' + ' / ' + r'$[10^{-17} erg.s^{-1}.cm^{-2}])$')
axis11.set_ylabel(r'log([R3]$_{Medido}$' + ' / ' + r'$[10^{-17} erg.s^{-1}.cm^{-2}])$')

mini = min(np.log10((useful.OIII_5007_SDSS+useful.OIII_4959_SDSS)/useful.Hb_SDSS)-2)
maxi = max(np.log10((useful.OIII_5007+useful.OIII_4959)/useful.Hb)+2)

axis11.set_xlim(mini+1.9,maxi-1.9)
axis11.set_ylim(mini+1.9,maxi-1.9)

axis11.fill_between([mini,maxi], [mini-0.1,maxi-0.1], [mini+0.1,maxi+0.1], 
                   where=([mini-0.1,maxi-0.1] <= [mini+0.1,maxi+0.1]), color='gray', alpha=0.4, interpolate=True, label=r'$\pm 0.1 dex$')
axis11.plot([mini-0.1,maxi-0.1],[mini,maxi], color='black', alpha=0.5, ls='--')
axis11.plot([mini,maxi],[mini,maxi], color='black', alpha=0.5, ls='-')
axis11.plot([mini+0.1,maxi+0.1],[mini,maxi], color='black', alpha=0.5, ls='--')
axis11.legend()

axis11.set_title('Comparação da razão R3 = ([OIII'+r'$_{5007}$]+[OIII'+r'$_{4959}]$)/H'+r'$\beta$')

axis11.text(0.75,0.05,r'log$_{10}\left(\frac{F^{R3}_{Medido}}{F^{R3}_{SDSS}}\right) \Rightarrow$', fontsize=15, transform=axis11.transAxes, horizontalalignment='right', verticalalignment='bottom')
axis11.text(0.96,0.13,'Mediana: '+str(round(median_OIII_Hb,2)), fontsize=10, transform=axis11.transAxes, horizontalalignment='right', verticalalignment='bottom')
axis11.text(0.95,0.07,'Desvio: '+str(round(std_OIII_Hb,2)), fontsize=10, transform=axis11.transAxes, horizontalalignment='right', verticalalignment='bottom')
#fig11.savefig('/home/vitorbootz/research/TCC_images/flux_comparsion/flux_comparsion_R3.pdf')

#########################################
############ (NII1+NII2)/Hb ratio ##############


fig12 = plt.figure(figsize=(6,4))
axis12 = mosaic_OH.add_subplot(gs[0, 1])
axis12.grid(alpha=0.2, color='grey')

razao_NII_Hb = np.log10((useful.NII_6548+useful.NII_6583)/useful.Hb) - np.log10((useful.NII_6548_SDSS+useful.NII_6583_SDSS)/useful.Hb_SDSS)
median_NII_Hb = razao_NII_Hb.median()
std_NII_Hb = razao_NII_Hb.std()

axis12.plot(np.log10((useful.NII_6548_SDSS+useful.NII_6583_SDSS)/useful.Hb_SDSS),np.log10((useful.NII_6548+useful.NII_6583)/useful.Hb),'x', color='black', ms='9')
axis12.set_xlabel(r'log([N2]$_{SDSS}$' + ' / ' + r'$[10^{-17} erg.s^{-1}.cm^{-2}])$')
axis12.set_ylabel(r'log([N2]$_{Medido}$' + ' / ' + r'$[10^{-17} erg.s^{-1}.cm^{-2}])$')

mini = min(np.log10((useful.NII_6548_SDSS+useful.NII_6583_SDSS)/useful.Hb_SDSS)-2)
maxi = max(np.log10((useful.NII_6548+useful.NII_6583)/useful.Hb)+2)

axis12.set_xlim(mini+1.9,maxi-1.9)
axis12.set_ylim(mini+1.9,maxi-1.9)

axis12.fill_between([mini,maxi], [mini-0.1,maxi-0.1], [mini+0.1,maxi+0.1], 
                   where=([mini-0.1,maxi-0.1] <= [mini+0.1,maxi+0.1]), color='gray', alpha=0.4, interpolate=True, label=r'$\pm 0.1 dex$')
axis12.plot([mini-0.1,maxi-0.1],[mini,maxi], color='black', alpha=0.5, ls='--')
axis12.plot([mini,maxi],[mini,maxi], color='black', alpha=0.5, ls='-')
axis12.plot([mini+0.1,maxi+0.1],[mini,maxi], color='black', alpha=0.5, ls='--')
axis12.legend()

axis12.set_title('Comparação da razão N2 = ([NII'+r'$_{6548}$]+[NII'+r'$_{6583}]$)/H'+r'$\beta$')

axis12.text(0.75,0.05,r'log$_{10}\left(\frac{F^{N2}_{Medido}}{F^{N2}_{SDSS}}\right) \Rightarrow$', fontsize=15, transform=axis12.transAxes, horizontalalignment='right', verticalalignment='bottom')
axis12.text(0.96,0.13,'Mediana: '+str(round(median_NII_Hb,2)), fontsize=10, transform=axis12.transAxes, horizontalalignment='right', verticalalignment='bottom')
axis12.text(0.95,0.07,'Desvio: '+str(round(std_NII_Hb,2)), fontsize=10, transform=axis12.transAxes, horizontalalignment='right', verticalalignment='bottom')
#fig12.savefig('/home/vitorbootz/research/TCC_images/flux_comparsion/flux_comparsion_N2.pdf')

#########################################
############ (SII1+SII2)/Hb ratio ##############


fig13 = plt.figure(figsize=(6,4))
axis13 = mosaic_OH.add_subplot(gs[1, 0])
axis13.grid(alpha=0.2, color='grey')

razao_SII_Hb = np.log10((useful.SII_6718+useful.SII_6732)/useful.Hb) - np.log10((useful.SII_6718_SDSS+useful.SII_6732_SDSS)/useful.Hb_SDSS)
median_SII_Hb = razao_SII_Hb.median()
std_SII_Hb = razao_SII_Hb.std()

axis13.plot(np.log10((useful.SII_6718_SDSS+useful.SII_6732_SDSS)/useful.Hb_SDSS),np.log10((useful.SII_6718+useful.SII_6732)/useful.Hb),'x', color='black', ms='9')
axis13.set_xlabel(r'log([S2]$_{SDSS}$' + ' / ' + r'$[10^{-17} erg.s^{-1}.cm^{-2}])$')
axis13.set_ylabel(r'log([S2]$_{Medido}$' + ' / ' + r'$[10^{-17} erg.s^{-1}.cm^{-2}])$')

mini = min(np.log10((useful.SII_6718_SDSS+useful.SII_6732_SDSS)/useful.Hb_SDSS)-2)
maxi = max(np.log10((useful.SII_6718+useful.SII_6732)/useful.Hb)+2)

axis13.set_xlim(mini+1.9,maxi-1.9)
axis13.set_ylim(mini+1.9,maxi-1.9)

axis13.fill_between([mini,maxi], [mini-0.1,maxi-0.1], [mini+0.1,maxi+0.1], 
                   where=([mini-0.1,maxi-0.1] <= [mini+0.1,maxi+0.1]), color='gray', alpha=0.4, interpolate=True, label=r'$\pm 0.1 dex$')
axis13.plot([mini-0.1,maxi-0.1],[mini,maxi], color='black', alpha=0.5, ls='--')
axis13.plot([mini,maxi],[mini,maxi], color='black', alpha=0.5, ls='-')
axis13.plot([mini+0.1,maxi+0.1],[mini,maxi], color='black', alpha=0.5, ls='--')
axis13.legend()

axis13.set_title('Comparação da razão S2 = ([SII'+r'$_{6718}$]+[SII'+r'$_{6732}]$)/H'+r'$\beta$')

axis13.text(0.75,0.05,r'log$_{10}\left(\frac{F^{S2}_{Medido}}{F^{S2}_{SDSS}}\right) \Rightarrow$', fontsize=15, transform=axis13.transAxes, horizontalalignment='right', verticalalignment='bottom')
axis13.text(0.96,0.13,'Mediana: '+str(round(median_SII_Hb,2)), fontsize=10, transform=axis13.transAxes, horizontalalignment='right', verticalalignment='bottom')
axis13.text(0.95,0.07,'Desvio: '+str(round(std_SII_Hb,2)), fontsize=10, transform=axis13.transAxes, horizontalalignment='right', verticalalignment='bottom')
#fig13.savefig('/home/vitorbootz/research/TCC_images/flux_comparsion/flux_comparsion_S2.pdf')

####################################################

mosaic.savefig('/home/vitorbootz/research/TCC_images/flux_comparsion/mosaico.jpg', dpi=300)
mosaic_bpt.savefig('/home/vitorbootz/research/TCC_images/flux_comparsion/mosaico_bpt.jpg', dpi=300)
mosaic_OH.savefig('/home/vitorbootz/research/TCC_images/flux_comparsion/mosaico_OH.jpg', dpi=300)

plt.show()


