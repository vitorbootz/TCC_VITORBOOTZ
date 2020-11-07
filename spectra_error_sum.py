import numpy as np
from astropy.io import fits

def calc_erro(arquivo_1, arquivo_2, arquivo_3): #Define a funcao que calcula a soma quadratica dos erros dos 3 arquivos
    hdul_1 = fits.open(arquivo_1) #Abre o arquivo Header Data Unit List, formado por um header e um data
    hdul_2 = fits.open(arquivo_2)
    hdul_3 = fits.open(arquivo_3)
    error_1 = hdul_1[0].data #Selecionamos a parte data e transormamos em um array numpy
    error_2 = hdul_2[0].data
    error_3 = hdul_3[0].data
    erro_lcg = np.sqrt(error_1**2 + error_2**2 + error_3**2) / 3 #Calculamos o erro
    return(erro_lcg) #Retornamos o resultado
    
grupo = int(input("ID do grupo: "))
extensoes = int(input("Quantidade de extensões: "))
for i in range(2,extensoes+1):
    arquivo_1 = str('error_LCG' + str(grupo) + '_1_' + str(i) + '.fits')
    arquivo_2 = str('error_LCG' + str(grupo) + '_2_' + str(i) + '.fits')
    arquivo_3 = str('error_LCG' + str(grupo) + '_3_' + str(i) + '.fits')
    erro_lcg = calc_erro(arquivo_1, arquivo_2, arquivo_3) #Chama a funcao
    hdu = fits.PrimaryHDU(erro_lcg) #Transforma de array numpy para HDUL novamente
    output = str('ERROR_LCG' + str(grupo) + '_' + str(i) + '.fits')
    hdu.writeto(output)
