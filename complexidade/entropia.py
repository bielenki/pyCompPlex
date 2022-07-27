import numpy as np
import os
from numba import njit, prange, jit
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
                prob=[(mascara[mascara==i]).size/(lenVet*1.0) for i in Lista]
                He = np.sum([(-1.0*p*np.log2(p)) for p in prob if p>0])
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
                    D = np.sum([((p-(1/N))**2) for p in prob])
                    LMC=D*C
                    E[row,col]=LMC
    return E
@jit
def convolucaoCube (listArray, rows, cols, imagens, kernel):
    arrayHe = np.empty((rows,cols), dtype=float)
    arrayHemax = np.empty((rows,cols), dtype=float)
    arraySDL = np.empty((rows,cols), dtype=float)
    arrayLMC = np.empty((rows,cols), dtype=float)
    for row in prange(rows):
        for col in prange(cols):
            janela = kernel
            if kernel==1:
                mascara=[]
                for imagem in range(imagens):
                    mascara.append(listArray[imagem][row][col])
            else:
                Lx=max(0,col-janela)
                Ux=min(cols,col+janela+1)
                Ly=max(0,row-janela)
                Uy=min(rows,row+janela+1)
                mascara=[]
                mascara1 = []
                for imagem in range(imagens):
                    imArray = listArray[imagem]
                    mascara1 = list( imArray[Ly:Uy,Lx:Ux].flatten())
                    mascara = mascara + mascara1
                    mascara1=[]
            He=0.0
            lenVet=len(mascara)
            Lista=list(set(mascara))
            if len(Lista)==1 and Lista.count(0)==1:
                arrayHe[row,col]=0
                arrayHemax[row,col]=0
                arraySDL[row,col]=0
                arrayLMC[row,col]=0
            else:
                prob=[(mascara.count(i))/(lenVet*1.0) for i in Lista]
                He = np.sum([(-1.0*p*np.log2(p)) for p in prob if p>0])
                arrayHe[row,col]=He
                N=len(Lista)*1.0
                if N == 1:
                    C=0
                else:
                    Hmax=np.log2(N)
                    C=He/Hmax
                arrayHemax[row,col]=C
                SDL=(1-C)*C
                arraySDL[row,col]=SDL
                D = 0.0
                D = np.sum([((p-(1/N))**2) for p in prob])
                LMC=D*C
                arrayLMC[row,col]=LMC
    return (arrayHe, arrayHemax, arraySDL, arrayLMC)