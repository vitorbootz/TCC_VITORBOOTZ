import numpy as np
import pandas as pd

def red_calc(data):
    vel_mean = sum(data['vel']*(1/pow(data['error'],2))) / sum(1/pow(data['error'],2)) #Calculo da média ponderada pelo erro
    redshift = vel_mean / 299792.458 #Convertendo de velocidade radial (km/s) para redshift
    return redshift
    
group = int(input("ID group: "))
extensions = int(input("Number of extensions: "))
duplicates = int(input("Quantidade de extensões duplicadas: ")) #Contador
print("############## REDSHIFTS {} ##############".format(group))
for i in range(2,extensions+1): #Iteração que cobre todas as galáxias do conjunto de dados
    name = str(group) + '_' + str(i) + '.csv' #Nome: NumeroDoGrupo_Extensão.csv
    file = pd.read_csv(name, header=None, delim_whitespace=True, names = ['vel', 'error'])
    if sum(file['error']) != 0: #Ignoramos as extensões completamente nulas
        selection = (file['error'] != 0) #Selecionamos apenas as linhas não nulas
        data = file[selection] #Criamos uma nova variável com linhas não nulas, exclusivamente
        redshift = red_calc(data) #Chamamos a função red_calc
        print("Slit " + str(i) + ": {}".format(round(redshift, 5)))
    else:
        print("Slit " + str(i) + ": -99")
    if duplicates > 0: #Esta nova iteração verifica quais são as extensões duplicadas
        try: #A função tentará o conjunto de comandos abaixo, se ocorrer erro o código pula direto para except
            name = str(group) + '_' + str(i) + '_2.csv'
            file = pd.read_csv(name, header=None, delim_whitespace=True, names = ['vel', 'error'])
            if sum(file['error']) != 0:
                selection = (file['error'] != 0)
                data = file[selection]
                redshift = red_calc(data)
                print("Slit " + str(i) + "_2: {}".format(round(redshift, 5)))
            else:
                print("Slit " + str(i) + ": -99")
            duplicates -= 1 #Se esta sequência de comandos for bem sucedida, então encontramos uma extensão duplicada
        except FileNotFoundError as error: #Se este erro específico ocorrer, exijo que a interação prossiga sem erros
            continue
