import numpy as np
from numba import prange, jit
@jit
def convolucaoNumba (ImArray, rows, cols, janela, opcao):
    for row in prange(rows):
        for col in prange(cols):
            Lx=max(0,col-janela)
            Ux=min(cols,col+janela+1)
            Ly=max(0,row-janela)
            Uy=min(rows,row+janela+1)
            mascara=ImArray[Ly:Uy,Lx:Ux].flatten()
            lenVet=mascara.size
            Lista=list(set(mascara))
            if len(Lista)==1 and Lista.count(0)==1:
                He=0
                C=0
                SDL=0
                LMC=0
            else:
                prob = np.zeros(len(Lista))
                for i in range(len(Lista)):
                    prob[i] = np.sum(mascara == Lista[i]) / lenVet
                He = np.sum(-1.0*prob*np.log2(prob)) #for p in prob if p>0])
                N=len(Lista)*1.0
                if N == 1:
                    C=0
                else:
                    Hmax=np.log2(N)
                    C=He/Hmax
                SDL=(1-C)*C
                D = np.sum((prob - (1 / N)) ** 2) if N > 1 else 0
                LMC=D*C

    return He,C, SDL, LMC
