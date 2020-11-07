import numpy as np
import pandas as pd
import scipy
import seaborn as sns
import matplotlib.pyplot as plt
from scipy import stats

sample = pd.read_csv('/home/vitorbootz/research/aux_files/abundancias_sample.csv')
control = pd.read_csv('/home/vitorbootz/research/aux_files/abundancias_control.csv')
sdss = pd.read_csv('/home/vitorbootz/research/aux_files/abundancias_sdss.csv')

fig, ax = plt.subplots()
ax.grid(axis='y', alpha=0.8)

control.OH.plot.hist(density=True, ax=ax, histtype='step', color='#EBA592', lw=3, fill=True, alpha=0.9, bins=5, edgecolor='#EBA592', label='LCGs isoladas')
sample.OH.plot.hist(density=True, ax=ax, histtype='step', color='#429E83', linewidth=3, hatch='/', alpha=0.9, bins=5, label='LCGs em grupos')

ax1 = sns.kdeplot(control.OH, color='#EB5E5C', label='_nolegend_')
ax2 = sns.kdeplot(sample.OH, color='#429E83', shade=False, ls='--', label='_nolegend_') 
stats, pvalue = stats.ks_2samp(control.OH, sample.OH, mode = "asymp")

m1 = ax1.axvline(control.OH.median(), c='#EB5E5C', ls='--', lw=1)
m2 = ax2.axvline(sample.OH.median(), c='#429E83', ls='--', lw=1)

ax.text(0.65,0.9, 'KS '+r'$\it{p}$'+'-value = ' + str(round(pvalue,2)), fontsize=12, transform=ax.transAxes)

ax.set_xlabel('12+log(O/H)', fontsize=15)
ax.set_ylabel('Densidade', fontsize=15)
ax.legend(loc='upper left')

fig.savefig('/home/vitorbootz/research/TCC_images/abundancia_oxigenio/hist_abundancia_LCGs_sample_control.pdf', format='pdf', bbox_inches='tight')


sample = pd.read_csv('/home/vitorbootz/research/aux_files/abundancias_sample.csv')
control = pd.read_csv('/home/vitorbootz/research/aux_files/abundancias_control.csv')
sdss = pd.read_csv('/home/vitorbootz/research/aux_files/abundancias_sdss.csv')
sample_gemini = pd.read_csv('/home/vitorbootz/research/aux_files/abundancias_sample_lcgs_gemini.csv')

sample_toplot = pd.DataFrame([])
sdss_toplot = pd.DataFrame([])
for i in range(len(sample_gemini)):
    for j in range(len(sample)):
        if (sample_gemini.lcgID[i] == sample.lcgID[j]):
            sample_toplot = sample_toplot.append(sample[sample.lcgID == sample_gemini.lcgID[i]])
    for j in range(len(sdss)):
        if (sample_gemini.lcgID[i] == sdss.lcgID[j]):
            sdss_toplot = sdss_toplot.append(sdss[sdss.lcgID == sample_gemini.lcgID[i]])
    
sdss_toplot.index = range(len(sdss_toplot))
sample_toplot.index = range(len(sample_toplot))

labels = sample_toplot.lcgID

x = np.arange(len(labels))  # the label locations
width = 0.15  # the width of the bars

fig, ax = plt.subplots()
rects1 = ax.bar(x - 1*width, sample_gemini.OH, width, label='Fluxos dos espectros do Gemini')
rects2 = ax.bar(x, sample_toplot.OH, width, label='Fluxos dos espectros do SDSS')
rects3 = ax.bar(x + 1*width, sdss_toplot.OH, width, label='Fluxos obtidos diretamente do SDSS')

ax.set_ylabel('12+log(O/H)', size=12)
ax.set_xlabel('LCG ID', size=12)
ax.set_xticks(x)
ax.set_xticklabels(labels)
ax.set_ylim(7,8.5)
ax.legend(loc='lower right', shadow=True, fontsize='medium')

def autolabel(rects):
    """Attach a text label above each bar in *rects*, displaying its height."""
    for rect in rects:
        height = rect.get_height()


autolabel(rects1)
autolabel(rects2)

plt.grid(alpha=0.3, color='black')
ax.set_axisbelow(True)
plt.yticks(np.arange(7, 8.5, 0.1))
ax.tick_params(axis='both', which='major', labelsize=12)
