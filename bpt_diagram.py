#############################
# BPT DIAGRAM
#############################

"""
SDSS Line-ratio Diagrams
------------------------
This shows how to plot line-ratio diagrams for the SDSS spectra.  These
diagrams are often called BPT plots [1]_, Osterbrock diagrams [2]_,
or Kewley diagrams [3]_. The location of the dividing line is taken from
from Kewley et al 2001.

References
~~~~~~~~~~
.. [1] Baldwin, J. A.; Phillips, M. M.; Terlevich, R. (1981)
       http://adsabs.harvard.edu/abs/1981PASP...93....5B
.. [2] Osterbrock, D. E.; De Robertis, M. M. (1985)
       http://adsabs.harvard.edu/abs/1985PASP...97.1129O
.. [3] Kewley, L. J. `et al.` (2001)
       http://adsabs.harvard.edu/abs/2001ApJ...556..121K
"""
# Author: Jake VanderPlas <vanderplas@astro.washington.edu>
# License: BSD
#   The figure is an example from astroML: see http://astroML.github.com
import numpy as np
import pandas as pd
import random
from matplotlib import pyplot as plt
import astropy.io.fits as fits
from astroML.datasets import fetch_sdss_corrected_spectra
from astroML.datasets.tools.sdss_fits import log_OIII_Hb_NII

data = fetch_sdss_corrected_spectra()

i = np.where((data['lineindex_cln'] == 4) | (data['lineindex_cln'] == 5))

FITS = fits.open('/home/vitorbootz/research/TCC_images/bpt/gal_line_dr7_v5_2.fit')  # Read the compress data
sample = pd.read_csv('/home/vitorbootz/research/flux_measurements/sample_flux_lines.csv')
control = pd.read_csv('/home/vitorbootz/research/aux_files/control_sample_flux_lines.csv')

# Jogados fora: grupo 2361, 2023_1 e 3090_2
useful = sample[(sample.lcgID != 2361) & (sample.lcgID != 1864) & (sample.extension != 13) & (sample.ra != 131.36510) & (sample.ra != 183.40380) & (sample.ra != 164.07920)]
control = control[control.lcgID != 1864]

useful_lcgs = useful[useful.flag_lcg == 1]
useful_others = useful[useful.flag_lcg == 0]

sample_NII_Ha_lcgs = np.log10(useful_lcgs.NII_6583/useful_lcgs.Ha)
sample_OIII_Hb_lcgs = np.log10(useful_lcgs.OIII_5008/useful_lcgs.Hb)
sample_NII_Ha_others = np.log10(useful_others.NII_6583/useful_others.Ha)
sample_OIII_Hb_others = np.log10(useful_others.OIII_5008/useful_others.Hb)
control_NII_Ha = np.log10(control.nii_6584_flux/control.h_alpha_flux)
control_OIII_Hb = np.log10(control.oiii_5007_flux/control.h_beta_flux)


head      = FITS[1].data
infheader = head.columns

Hbeta  = FITS[1].data['H_BETA_FLUX']          # Reading the column with Hbeta line 
OIII   = FITS[1].data['OIII_5007_FLUX']       # Reading the column with OIII line
Halpha = FITS[1].data['H_ALPHA_FLUX']         # Reading the column with Halpha line
NII    = FITS[1].data['NII_6584_FLUX']   

NII_Ha = np.log10(NII/Halpha)
OIII_Hb = np.log10(OIII/Hbeta)
f = plt.figure(1)
ax = f.add_subplot(1,1,1)

plt.scatter(NII_Ha, OIII_Hb,
            c='grey', s=0.1, lw=1, alpha=0.05)
plt.scatter(control_NII_Ha, control_OIII_Hb,
            c='#008D33', s=100, lw=0.3, marker='d', label = 'LCGs isoladas')
plt.scatter(sample_NII_Ha_lcgs, sample_OIII_Hb_lcgs,
            c='#00BCD9', s=100, lw=2, marker='x', label = 'LCGs em grupos')
plt.scatter(sample_NII_Ha_others, sample_OIII_Hb_others,
            c='#EB5E5C', s=100, lw=1, marker='+', label = 'Galáxias vizinhas')


# Kewley+01 ------------------------------------------
X = np.linspace(-2,0.3)
Y = (0.61/( X  - 0.47  )) + 1.19

# Schawinski+07 --------------------------------------
X3 = np.linspace(-0.180,1.5)
Y3 = 1.05*X3 + 0.45

# Kauffmann+03 ---------------------------------------
Xk = np.linspace(-2,0.)
Yk = 0.61/(Xk -0.05) + 1.3

# Regions --------------------------------------------
ax.plot(X,   Y, '-' , color='black', alpha=0.8, lw=2, label='Kewley+01'    ) # Kewley+01
ax.plot(X3, Y3, '-.', color='black', alpha=0.8, lw=2, label='Schawinski+07') # Schawinski+07
ax.plot(Xk, Yk, '--', color='black', alpha=0.8, lw=2, label='Kauffmann+03' ) # Kauffmann+03

ax.legend(loc='lower left', shadow=True, fontsize='medium')

plt.xlim(-2.0, 1.0)
plt.ylim(-1.2, 1.5)

plt.xlabel(r'$\mathrm{log_{10}([NII]/H\alpha)}$', fontsize='large')
plt.ylabel(r'$\mathrm{log_{10}([OIII]/H\beta)}$', fontsize='large')

plt.savefig('/home/vitorbootz/research/TCC_images/bpt/bpt_diagram.png', format='png', dpi=400, bbox_inches='tight')

plt.show()
