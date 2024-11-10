import numpy as np
from numba import prange, jit, njit
import timeit

def divide_em_blocos(imagem, num_blocos):
    blocos = []
    bloco_tamanho = imagem.shape[0] // num_blocos
    for i in range(0, imagem.shape[0], bloco_tamanho):
        bloco = imagem[i:i + bloco_tamanho, :]
        blocos.append(bloco)
    return blocos

@njit#(parallel=True)
def calcular_He(mascara, Lista, lenVet):
    prob = np.zeros(len(Lista))
    for i in range(len(Lista)):
        prob[i] = np.sum(mascara == Lista[i]) / lenVet
    He = 0.0
    for p in prob:
        if p > 0:
            He += -1.0 * p * np.log2(p)
    return He, prob

@njit#(parallel=True)
def convolucaoNumba1(ImArray, rows, cols, janela):
    arrayHe = np.empty((rows, cols), dtype=np.float64)
    arrayHemax = np.empty((rows, cols), dtype=np.float64)
    arraySDL = np.empty((rows, cols), dtype=np.float64)
    arrayLMC = np.empty((rows, cols), dtype=np.float64)
    for row in prange(rows):
        for col in prange(cols):
            Lx = max(0, col - janela)
            Ux = min(cols, col + janela + 1)
            Ly = max(0, row - janela)
            Uy = min(rows, row + janela + 1)
            mascara = ImArray[Ly:Uy, Lx:Ux].flatten()
            lenVet = mascara.size
            Lista = list(set(mascara))
            if len(Lista) == 1 and Lista.count(0) == 1:
                arrayHe[row, col] = 0
                arrayHemax[row, col] = 0
                arrayLMC[row, col] = 0
                arraySDL[row, col] = 0
            else:
                prob = np.zeros(len(Lista))
                for i in range(len(Lista)):
                    prob[i] = np.sum(mascara == Lista[i]) / lenVet
                #He = np.sum([(-1.0 * p * np.log2(p))  for p in prob if p>0])
                He = np.sum(-1.0 * prob * np.log2(prob))
                #He, prob = calcular_He(mascara, Lista, lenVet)
                N = len(Lista) * 1.0
                if N == 1:
                    C = 0
                else:
                    Hmax = np.log2(N)
                    C = He / Hmax
                SDL = (1 - C) * C
                D = np.sum((prob - (1 / N)) ** 2) if N > 1 else 0
                LMC = D * C
                arrayHe[row, col] = He
                arrayHemax[row, col] = C
                arrayLMC[row, col] = LMC
                arraySDL[row, col] = SDL
    return arrayHe, arrayHemax, arraySDL, arrayLMC

@njit#(parallel=True)
def convolucaoNumba2(ImArray, rows, cols, janela):
    arrayHe = np.empty((rows, cols), dtype=np.float64)
    arrayHemax = np.empty((rows, cols), dtype=np.float64)
    arraySDL = np.empty((rows, cols), dtype=np.float64)
    arrayLMC = np.empty((rows, cols), dtype=np.float64)
    for row in prange(rows):
        for col in prange(cols):
            Lx = max(0, col - janela)
            Ux = min(cols, col + janela + 1)
            Ly = max(0, row - janela)
            Uy = min(rows, row + janela + 1)
            mascara = ImArray[Ly:Uy, Lx:Ux].flatten()
            lenVet = mascara.size
            Lista = list(set(mascara))
            if len(Lista) == 1 and Lista.count(0) == 1:
                arrayHe[row, col] = 0
                arrayHemax[row, col] = 0
                arrayLMC[row, col] = 0
                arraySDL[row, col] = 0
            else:
                #prob = np.zeros(len(Lista))
                #for i in range(len(Lista)):
                #    prob[i] = np.sum(mascara == Lista[i]) / lenVet
                #He = np.sum([(-1.0 * p * np.log2(p))  for p in prob if p>0])
                #He = np.sum(-1.0 * prob * np.log2(prob))
                He, prob = calcular_He(mascara, Lista, lenVet)
                N = len(Lista) * 1.0
                if N == 1:
                    C = 0
                else:
                    Hmax = np.log2(N)
                    C = He / Hmax
                SDL = (1 - C) * C
                D = np.sum((prob - (1 / N)) ** 2) if N > 1 else 0
                LMC = D * C
                arrayHe[row, col] = He
                arrayHemax[row, col] = C
                arrayLMC[row, col] = LMC
                arraySDL[row, col] = SDL
    return arrayHe, arrayHemax, arraySDL, arrayLMC

@njit#(parallel=True)
def convolucaoNumba3(ImArray, rows, cols, janela):
    arrayHe = np.empty((rows, cols), dtype=np.float64)
    arrayHemax = np.empty((rows, cols), dtype=np.float64)
    arraySDL = np.empty((rows, cols), dtype=np.float64)
    arrayLMC = np.empty((rows, cols), dtype=np.float64)
    for row in prange(rows):
        for col in prange(cols):
            Lx = max(0, col - janela)
            Ux = min(cols, col + janela + 1)
            Ly = max(0, row - janela)
            Uy = min(rows, row + janela + 1)
            mascara = ImArray[Ly:Uy, Lx:Ux].flatten()
            lenVet = mascara.size
            Lista = list(set(mascara))
            if len(Lista) == 1 and Lista.count(0) == 1:
                arrayHe[row, col] = 0
                arrayHemax[row, col] = 0
                arrayLMC[row, col] = 0
                arraySDL[row, col] = 0
            else:
                prob = np.zeros(len(Lista))
                for i in range(len(Lista)):
                    prob[i] = np.sum(mascara == Lista[i]) / lenVet
                He = np.sum([(-1.0 * p * np.log2(p))  for p in prob if p>0])
                #He = np.sum(-1.0 * prob * np.log2(prob))
                #He, prob = calcular_He(mascara, Lista, lenVet)
                N = len(Lista) * 1.0
                if N == 1:
                    C = 0
                else:
                    Hmax = np.log2(N)
                    C = He / Hmax
                SDL = (1 - C) * C
                D = np.sum((prob - (1 / N)) ** 2) if N > 1 else 0
                LMC = D * C
                arrayHe[row, col] = He
                arrayHemax[row, col] = C
                arrayLMC[row, col] = LMC
                arraySDL[row, col] = SDL
    return arrayHe, arrayHemax, arraySDL, arrayLMC

@njit#(parallel=True)
def convolucaoNumba4(ImArray, rows, cols, janela):
    arrayHe = np.empty((rows, cols), dtype=np.float64)
    arrayHemax = np.empty((rows, cols), dtype=np.float64)
    arraySDL = np.empty((rows, cols), dtype=np.float64)
    arrayLMC = np.empty((rows, cols), dtype=np.float64)
    for row in prange(rows):
        for col in prange(cols):
            Lx = max(0, col - janela)
            Ux = min(cols, col + janela + 1)
            Ly = max(0, row - janela)
            Uy = min(rows, row + janela + 1)
            mascara = ImArray[Ly:Uy, Lx:Ux].flatten()
            lenVet = mascara.size
            Lista = list(set(mascara))
            if len(Lista) == 1 and Lista.count(0) == 1:
                arrayHe[row, col] = 0
                arrayHemax[row, col] = 0
                arrayLMC[row, col] = 0
                arraySDL[row, col] = 0
            else:
                prob = np.zeros(len(Lista))
                for i in range(len(Lista)):
                    prob[i] = np.sum(mascara == Lista[i]) / lenVet
                He = 0.0
                for p in prob:
                    if p > 0:
                        He += -1.0 * p * np.log2(p)
                N = len(Lista) * 1.0
                if N == 1:
                    C = 0
                else:
                    Hmax = np.log2(N)
                    C = He / Hmax
                SDL = (1 - C) * C
                D = np.sum((prob - (1 / N)) ** 2) if N > 1 else 0
                LMC = D * C
                arrayHe[row, col] = He
                arrayHemax[row, col] = C
                arrayLMC[row, col] = LMC
                arraySDL[row, col] = SDL
    return arrayHe, arrayHemax, arraySDL, arrayLMC

def minha_funcao(parametros, algoritmo='algoritmo1'):
    if algoritmo == 'algoritmo1':
        convolucaoNumba1(*parametros)
    elif algoritmo == 'algoritmo2':
        convolucaoNumba2(*parametros)
    elif algoritmo == 'algoritmo3':
        convolucaoNumba3(*parametros)
    elif algoritmo == 'algoritmo4':
        convolucaoNumba4(*parametros)


linhas=[100, 500, 1000]
janelas=[1,2,3]
for jan in janelas:
    for linha in linhas:
        # Defina uma função que executa o algoritmo específico
        def teste_algoritmo(algoritmo):
            # Aqui você passa os parâmetros necessários para sua função
            arrayR = np.random.randint(0, 255, (linha, linha))
            janela = jan
            rows, cols = arrayR.shape
            parametros=(arrayR,rows,cols,janela)
            minha_funcao(parametros, algoritmo=algoritmo)

        # Meça o tempo de execução para cada algoritmo
        tempo_algoritmo1 = timeit.timeit(lambda: teste_algoritmo('algoritmo1'), number=100)
        tempo_algoritmo2 = timeit.timeit(lambda: teste_algoritmo('algoritmo2'), number=100)
        #tempo_algoritmo3 = timeit.timeit(lambda: teste_algoritmo('algoritmo3'), number=100)
        tempo_algoritmo4 = timeit.timeit(lambda: teste_algoritmo('algoritmo4'), number=100)

        # Exiba os resultados
        print("linhas: "+str(linha)+" Janela: "+ str(1+2*(jan))+"x"+str(1+2*(jan)))
        print(f"Tempo do Algoritmo 1: {tempo_algoritmo1:.4f} segundos")
        print(f"Tempo do Algoritmo 2: {tempo_algoritmo2:.4f} segundos")
        #print(f"Tempo do Algoritmo 3: {tempo_algoritmo3:.4f} segundos")
        print(f"Tempo do Algoritmo 4: {tempo_algoritmo4:.4f} segundos")
        print()


