import numpy as np
import os
from numba import njit, prange, jit

@njit
def calcular_He(mascara, Lista, lenVet):
    prob = np.zeros(len(Lista))
    for i in range(len(Lista)):
        prob[i] = np.sum(mascara == Lista[i]) / lenVet
    He = 0.0
    for p in prob:
        if p > 0:
            He += -1.0 * p * np.log2(p)
    return He, prob

@jit
def convolucaoNumba (E, ImArray, rows, cols, janela, opcao):
    for row in prange(rows):
        for col in prange(cols):
            Lx=max(0,col-janela)
            Ux=min(cols,col+janela+1)
            Ly=max(0,row-janela)
            Uy=min(rows,row+janela+1)
            mascara=ImArray[Ly:Uy,Lx:Ux].flatten()
            He=0.0
            lenVet=mascara.size
            Lista=list(set(mascara))
            if len(Lista)==1 and Lista.count(0)==1:
                E[row,col]=0
            else:
                He, prob = calcular_He(mascara, Lista, lenVet)
                if opcao==0:
                    E[row,col]=He
                N=len(Lista)*1.0
                if N == 1:
                    C=0
                else:
                    Hmax=np.log2(N)
                    C=He/Hmax
                if opcao==1:
                    E[row,col]=C
                if opcao==2:
                    SDL=(1-C)*C
                    E[row,col]=SDL
                if opcao==3:
                    D = 0.0
                    D = np.sum((prob - (1 / N)) ** 2) if N > 1 else 0
                    LMC=D*C
                    E[row,col]=LMC
    return E

@jit(nopython=True, parallel=True)
def convolucaoCube(listArray, rows, cols, imagens, kernel):
    arrayHe = np.empty((rows, cols), dtype=np.float64)
    arrayHemax = np.empty((rows, cols), dtype=np.float64)
    arraySDL = np.empty((rows, cols), dtype=np.float64)
    arrayLMC = np.empty((rows, cols), dtype=np.float64)
    for row in prange(rows):
        for col in prange(cols):
            janela = kernel
            if kernel == 0:
                # Usa uma lista NumPy em vez de uma lista Python
                mascara = np.empty(imagens, dtype=np.float64)
                for imagem in range(imagens):
                    mascara[imagem] = listArray[imagem][row][col]
            else:
                Lx = max(0, col - janela)
                Ux = min(cols, col + janela + 1)
                Ly = max(0, row - janela)
                Uy = min(rows, row + janela + 1)
                # Usa uma lista NumPy para mascara
                mascara = np.empty(0, dtype=np.float64)
                for imagem in range(imagens):
                    imArray = listArray[imagem]
                    # Usa np.concatenate para unir as partes da máscara
                    mascara = np.concatenate((mascara, imArray[Ly:Uy, Lx:Ux].flatten()))
            He = 0.0
            lenVet = len(mascara)
            # Lista única com NumPy
            Lista = np.unique(mascara)
            if len(Lista) == 1 and Lista[0] == 0:
                arrayHe[row, col] = 0
                arrayHemax[row, col] = 0
                arraySDL[row, col] = 0
                arrayLMC[row, col] = 0
            else:
                # Cálculo de probabilidades com NumPy
                prob = np.array([(mascara == i).sum() / lenVet for i in Lista])
                # Cálculo manual de He sem usar list comprehension e np.sum
                He = 0.0
                for p in prob:
                    if p > 0:
                        He += -1.0 * p * np.log2(p)  # Acumulando o valor de He
                arrayHe[row, col] = He
                N = len(Lista) * 1.0
                if N == 1:
                    C = 0
                else:
                    Hmax = np.log2(N)
                    C = He / Hmax
                arrayHemax[row, col] = C
                SDL = (1 - C) * C
                arraySDL[row, col] = SDL
                # Cálculo de D com NumPy
                D = np.sum((prob - (1 / N)) ** 2) if N > 0 else 0
                LMC = D * C
                arrayLMC[row, col] = LMC
    return (arrayHe, arrayHemax, arraySDL, arrayLMC)