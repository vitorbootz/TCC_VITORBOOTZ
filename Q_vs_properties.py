import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
import seaborn as sns
from scipy.stats import spearmanr
from scipy.stats import kendalltau
from scipy.stats import pearsonr
pd.set_option('max_columns', None)
pd.set_option('max_rows', None)

prop = pd.read_csv('/home/vitorbootz/research/aux_files/properties_sample.csv')
tidal = pd.read_csv('/home/vitorbootz/research/aux_files/Q_values.csv')
oxig = pd.read_csv('/home/vitorbootz/research/aux_files/abundancias_sample.csv')
tidal = tidal[(tidal.lcgID != 2361) & (tidal.lcgID != 1864)]
oxig = oxig[(oxig.lcgID != 2361) & (tidal.lcgID != 1864)]
tidal.index = range(len(tidal))
oxig.index = range(len(oxig))
prop.index = range(len(prop))

#####################################
############ Concentração ################

fig2 = plt.figure(figsize=(6,4))
axis2 = fig2.add_subplot()
axis2.grid(alpha=0.2, color='grey')

razao_Q_OH = np.log10(abs(tidal.Q_group/prop.log_R90_R50))
median_Q_OH = razao_Q_OH.median()
std_Q_OH = razao_Q_OH.std()


#coef = np.polyfit(tidal.Q_group,prop.log_R90_R50,1)
#poly1d_fn = np.poly1d(coef)
#axis2.plot(tidal.Q_group, poly1d_fn(tidal.Q_group), '--k')

minix = min(tidal.Q_group)
maxix = max(tidal.Q_group)
miniy = min(prop.log_R90_R50)
maxiy = max(prop.log_R90_R50)
plt.xlim(minix-0.2,maxix+0.2)
plt.ylim(miniy-0.07,maxiy+0.02)

coef = np.polyfit(tidal.Q_group,prop.log_R90_R50,1)
m,b = np.poly1d(coef)
spearman, p = spearmanr(tidal.Q_group,prop.log_R90_R50)
kendall, tau = kendalltau(tidal.Q_group,prop.log_R90_R50)
pearson, pp = pearsonr(tidal.Q_group,prop.log_R90_R50)

sns.regplot(tidal.Q_group, prop.log_R90_R50, color='k', marker='x', ci=95, scatter_kws={"s": 1})

axis2.plot(tidal.Q_group,prop.log_R90_R50,'x', color='black', ms='9')
axis2.set_xlabel('Q', fontsize = 15)
axis2.set_ylabel(r'log$_{10}\left[\frac{petroR90}{petroR50}\right]$', fontsize = 15)

axis2.annotate('Coef. Ang. = '+str(round(m,3))+'\nI.C. = 95%', xy=(260, 40), xycoords='axes points',
            size=12, ha='center', va='top',
            bbox=dict(boxstyle='round', fc='w'))

axis2.annotate(r'Pearson $\rho$: '+ str(round(pearson,3)) + '\np-value: '+ str(round(pp,3)), xy=(72, 40), xycoords='axes points',
             size=12, ha='center', va='top',
             bbox=dict(boxstyle='round', fc='w'))


axis2.set_title('Q vs Concentração', fontsize=15)

plt.show()
fig2.savefig('/home/vitorbootz/research/TCC_images/Q_vs_properties/Q_concentracao.pdf', bbox_inches='tight')

#####################################
############ sSFR ################
fig3 = plt.figure(figsize=(6,4))
axis3 = fig3.add_subplot()
axis3.grid(alpha=0.2, color='grey')

razao_Q_OH = np.log10(abs(tidal.Q_group/prop.specsfr_tot_p50))
median_Q_OH = razao_Q_OH.median()
std_Q_OH = razao_Q_OH.std()


minix = min(tidal.Q_group)
maxix = max(tidal.Q_group)
miniy = min(prop.specsfr_tot_p50)
maxiy = max(prop.specsfr_tot_p50)
plt.xlim(minix-0.2,maxix+0.2)
plt.ylim(miniy-1,maxiy+0.3)

coef = np.polyfit(tidal.Q_group,prop.specsfr_tot_p50,1)
m,b = np.poly1d(coef)
spearman, p = spearmanr(tidal.Q_group,prop.specsfr_tot_p50)
kendall, tau = kendalltau(tidal.Q_group,prop.specsfr_tot_p50)
pearson, pp = pearsonr(tidal.Q_group,prop.specsfr_tot_p50)

sns.regplot(tidal.Q_group, prop.specsfr_tot_p50, color='k', marker='x', ci=95, scatter_kws={"s": 1})
axis3.plot(tidal.Q_group,prop.specsfr_tot_p50,'x', color='black', ms='9')
axis3.set_xlabel('Q', fontsize = 15)
axis3.set_ylabel('sSFR', fontsize = 15)

axis3.annotate('Coef. Ang. = '+str(round(m,3))+'\nI.C. = 95%', xy=(260, 40), xycoords='axes points',
            size=12, ha='center', va='top',
            bbox=dict(boxstyle='round', fc='w'))

axis3.annotate(r'Pearson $\rho$: '+ str(round(pearson,3)) + '\np-value: '+ str(round(pp,3)), xy=(72, 41), xycoords='axes points',
             size=12, ha='center', va='top',
             bbox=dict(boxstyle='round', fc='w'))

axis3.set_title('Q vs Taxa de formação estelar específica', fontsize=15)


fig3.savefig('/home/vitorbootz/research/TCC_images/Q_vs_properties/Q_sSFR.pdf', bbox_inches='tight')


#####################################
############ Abundância ################
tidal = tidal[(tidal.lcgID != 1864) & (tidal.lcgID != 2361) & (tidal.lcgID != 2023)]
tidal.index = range(len(tidal))
oxig = pd.read_csv('/home/vitorbootz/research/aux_files/abundancias_sample.csv')
oxig = oxig[(oxig.lcgID != 1864) & (oxig.lcgID != 2361) & (oxig.lcgID != 2023)]
oxig.index = range(len(oxig))

fig1 = plt.figure(figsize=(6,4))
axis1 = fig1.add_subplot()
axis1.grid(alpha=0.2, color='grey')

razao_Q_OH = np.log10(abs(tidal.Q_group/oxig.OH))
median_Q_OH = razao_Q_OH.median()
std_Q_OH = razao_Q_OH.std()

minix = min(tidal.Q_group)
maxix = max(tidal.Q_group)
miniy = min(oxig.OH)
maxiy = max(oxig.OH)
plt.xlim(minix-0.2,maxix+0.2)
plt.ylim(miniy-0.1,maxiy+0.02)

coef = np.polyfit(tidal.Q_group,oxig.OH,1)
m,b = np.poly1d(coef)
spearman, p = spearmanr(tidal.Q_group,oxig.OH)
kendall, tau = kendalltau(tidal.Q_group,oxig.OH)
pearson, pp = pearsonr(tidal.Q_group,oxig.OH)

sns.regplot(tidal.Q_group, oxig.OH, color='k', marker='x', ci=95, scatter_kws={"s": 1})
axis1.plot(tidal.Q_group,oxig.OH,'x', color='black', ms='9')
axis1.set_xlabel('Q', fontsize = 15)
axis1.set_ylabel('12 + log (O/H)', fontsize = 15)

axis1.set_title('Q vs Abundância de oxigênio', fontsize=15)

axis1.annotate('Coef. Ang. = '+str(round(m,3))+'\nI.C. = 95%', xy=(260, 40), xycoords='axes points',
            size=12, ha='center', va='top',
            bbox=dict(boxstyle='round', fc='w'))

axis1.annotate(r'Pearson $\rho$: '+ str(round(pearson,3)) + '\np-value: '+ str(round(pp,3)), xy=(72, 40), xycoords='axes points',
             size=12, ha='center', va='top',
             bbox=dict(boxstyle='round', fc='w'))


fig1.savefig('/home/vitorbootz/research/TCC_images/Q_vs_properties/Q_OH.pdf', bbox_inches='tight')


