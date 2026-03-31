import numpy as np

def simular_salto_kmc(energia_barreira, temperatura):
    kb = 8.617e-5  # Constante de Boltzmann em eV/K
    nu = 1e13      # Frequência de tentativa (Hz)
    
    # 1. Calcular a taxa de transição (Arrhenius)
    taxa = nu * np.exp(-energia_barreira / (kb * temperatura))
    
    # 2. Gerar o tempo de residência (estocástico)
    u = np.random.rand()
    tempo_espera = -np.log(u) / taxa
    
    return tempo_espera

# Simulando uma sequência de saltos na sua HEA
tempos = []
for i in range(1000):
    # 'eb' viria do seu mapa de energias gerado com ASE + AFLOW
    eb = mapa_de_energias_da_hea[i] 
    t = simular_salto_kmc(eb, temperatura=300)
    tempos.append(t)

# Agora você converte esses tempos em bits (0 ou 1)
mediana = np.median(tempos)
bits = [1 if t > mediana else 0 for t in tempos]
