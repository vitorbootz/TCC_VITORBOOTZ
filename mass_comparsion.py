import pandas as pd
sample = pd.read_csv('/home/vitorbootz/research/aux_files/galaxy_list.csv')
sdss = pd.read_csv('/home/vitorbootz/research/flux_measurements/sdss_flux_lines.csv')
main = pd.read_csv('/home/vitorbootz/research/results_starlight/STARLIGHT_MAIN_RESULTS.csv', sep=',')
cut = pd.read_csv('/home/vitorbootz/research/results_starlight/CUT_STARLIGHT_MAIN_RESULTS.csv', sep=',')

pd.set_option('max_columns', None)
pd.set_option('max_rows', None)

sample['lgm_tot_p50'] = -99.9
for i in range(len(sample)):
    for j in range(len(sdss)):
        if round(sample.ra[i], 2) == round(sdss.ra[j], 2):
            sample.lgm_tot_p50[i] = sdss.lgm_tot_p50[j]
            
selection_gem = (sample['onoff'] == 1) & (sample['flag_lcg'] == 1) & (sample['extension'] != 1) & (sample['flag_sdss'] == 0) & (sample['extension'] < 100)
selection_sdss = (sample['onoff'] == 1) & (sample['flag_lcg'] == 1) & (sample['extension'] != 1) & (sample['flag_sdss'] == 1)

sample_gem = sample[selection_gem]
sample_sdss = sample[selection_sdss]
sample_gem.index = range(len(sample_gem))
sample_sdss.index = range(len(sample_sdss))

mass_star_gem = pd.DataFrame([])
mass_star_sdss = pd.DataFrame([])
mass_star_cut = pd.DataFrame([])
lgm_tot_p50 = pd.DataFrame([])

for i in range(len(sample_gem)):
    selection_gem = (main['lcgID'] == sample_gem.lcgID[i]) & (main['extension'] == sample_gem.extension[i])
    selection_sdss = (main['lcgID'] == sample_sdss.lcgID[i]) & (main['extension'] == sample_sdss.extension[i])
    selection_cut =  (cut['lcgID'] == sample_gem.lcgID[i]) & (cut['extension'] == sample_gem.extension[i]+1000)
    mass_star_gem = mass_star_gem.append(main[selection_gem])
    mass_star_sdss = mass_star_sdss.append(main[selection_sdss])
    mass_star_cut = mass_star_cut.append(cut[selection_cut])

mass_star_gem.index = range(len(mass_star_gem))
mass_star_sdss.index = range(len(mass_star_sdss))
mass_star_cut.index = range(len(mass_star_cut))

mass_comparsion = sample_sdss[['lcgID', 'logMstar_lcg_sdss', 'lgm_tot_p50']].drop_duplicates()
mass_comparsion.index = range(len(mass_comparsion))

mass_comparsion = pd.concat([mass_comparsion,round(mass_star_gem.Mcor_log_Mo, 2), round(mass_star_sdss.Mcor_log_Mo, 2), round(mass_star_cut.Mcor_log_Mo, 2)], axis=1)
mass_comparsion.columns = ['lcgID', 'lgm_Izotov', 'lgm_tot_p50_SDSS', 'lgm_starlight_gem', 'lgm_starlight_sdss', 'lgm_starlight_sdss_cut']

############
### PLOT ###
############

import matplotlib
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.backends.backend_pdf import PdfPages
pp = PdfPages('/home/vitorbootz/research/TCC_images/comparacao_massas/mass_comparsion.pdf')

labels = mass_comparsion.lcgID

x = np.arange(len(labels))  # the label locations
width = 0.15  # the width of the bars

fig, ax = plt.subplots()
rects1 = ax.bar(x - 2*width, mass_comparsion.lgm_tot_p50_SDSS, width, label='lgm_tot_p50 SDSS')
rects2 = ax.bar(x - 1*width, mass_comparsion.lgm_Izotov, width, label='Izotov et al.')
rects3 = ax.bar(x, mass_comparsion.lgm_starlight_gem, width, label='STARLIGHT Gemini')
rects4 = ax.bar(x + 1*width, mass_comparsion.lgm_starlight_sdss, width, label='STARLIGHT SDSS')
rects5 = ax.bar(x + 2*width, mass_comparsion.lgm_starlight_sdss_cut, width, label='STARLIGHT SDSS-Gemini')



# Add some text for labels, title and custom x-axis tick labels, etc.
ax.set_ylabel('log(M$\star/$M$\odot$)', size=12)
ax.set_xlabel('LCG ID', size=12)
ax.set_xticks(x)
ax.set_xticklabels(labels)
ax.set_ylim(7,10.5)
ax.legend(loc='lower right', shadow=True, fontsize='x-small')

def autolabel(rects):
    """Attach a text label above each bar in *rects*, displaying its height."""
    for rect in rects:
        height = rect.get_height()


autolabel(rects1)
autolabel(rects2)

plt.grid(alpha=0.3, color='black')
ax.set_axisbelow(True)
plt.yticks(np.arange(7, 11, 0.5))
ax.tick_params(axis='both', which='major', labelsize=12)


plt.savefig(pp, format='pdf', figsize=(6,3))
pp.close()


