import argparse
import numpy as np
from numba import prange, njit
from osgeo import gdal
import time
from concurrent.futures import ProcessPoolExecutor
import warnings
from numba.core.errors import NumbaWarning, NumbaPendingDeprecationWarning
import multiprocessing as mp
warnings.filterwarnings("ignore", category=NumbaWarning)
warnings.filterwarnings("ignore", category=NumbaPendingDeprecationWarning)

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

@njit
def convolucaoNumba(E, bloco, janela, opcao, isFirst, isLast):
    rows, cols = bloco.shape
    if isFirst:
        rowInicio = 0
        rowFinal = rows - janela
    elif isLast:
        rowInicio = janela
        rowFinal = rows
    else:
        rowInicio = janela
        rowFinal = rows - janela
    for row in prange(rowInicio, rowFinal):  # Ajuste nas linhas
        for col in prange(cols):  # Sem ajuste nas colunas, pois elas não têm bordas
            Lx = max(0, col - janela)
            Ux = min(cols, col + janela + 1)
            Ly = max(0, row - janela)
            Uy = min(rows, row + janela + 1)
            mascara = bloco[Ly:Uy, Lx:Ux].ravel()
            He = 0.0
            lenVet = mascara.size
            Lista = list(set(mascara))
            if len(Lista) == 1 and Lista.count(0) == 1:
                if isFirst:
                    E[row, col] = 0
                else:
                    E[row - janela, col] = 0  # Ajuste no índice apenas para a linha
            else:
                He, prob = calcular_He(mascara, Lista, lenVet)
                if opcao == 0:
                    if isFirst:
                        E[row, col] = He
                    else:
                        E[row - janela, col] = He
                N = len(Lista) * 1.0
                if N == 1:
                    C = 0
                else:
                    Hmax = np.log2(N)
                    C = He / Hmax
                if opcao == 1:
                    if isFirst:
                        E[row, col] = C
                    else:
                        E[row - janela, col] = C
                SDL = (1 - C) * C
                if opcao == 2:
                    if isFirst:
                        E[row, col] = SDL
                    else:
                        E[row - janela, col] = SDL
                D = np.sum((prob - (1 / N)) ** 2) if N > 1 else 0
                LMC = D * C
                if opcao == 3:
                    if isFirst:
                        E[row, col] = LMC
                    else:
                        E[row - janela, col] = LMC
    return E

def divide_em_blocos(imagem, num_blocos, janela):
    blocos = []
    bloco_tamanho = imagem.shape[0] // num_blocos
    for i in range(0, imagem.shape[0], bloco_tamanho):
        bloco = imagem[max(0, i - janela):min(imagem.shape[0], i + bloco_tamanho + janela), :]
        blocos.append(bloco)
    return blocos

def processar_imagem_por_blocos(NameImagem, Janela, metrica, outFN, num_blocos=4):
    Im = gdal.Open(NameImagem)
    cols = Im.RasterXSize
    rows = Im.RasterYSize
    NrBandas = Im.RasterCount
    opcao = metrica
    with ProcessPoolExecutor() as executor:
        for band in range(1, NrBandas + 1):
            print("Processando banda " + str(band) + "...")
            banda_img = Im.GetRasterBand(band)
            ImArray = banda_img.ReadAsArray().astype(np.float32)
            blocos = divide_em_blocos(ImArray, num_blocos, Janela)
            resultados = []
            for idx, bloco in enumerate(blocos):
                # Verifica se o bloco é o primeiro ou o último para ajustar o tamanho de E
                isFirst = (idx == 0)
                isLast = (idx == len(blocos) - 1)
                if idx == 0:  # Primeiro bloco (sem borda superior)
                    E = np.empty((bloco.shape[0] - Janela, bloco.shape[1]))
                elif idx == len(blocos) - 1:  # Último bloco (sem borda inferior)
                    E = np.empty((bloco.shape[0] - Janela, bloco.shape[1]))
                else:  # Blocos intermediários (com borda completa)
                    E = np.empty((bloco.shape[0] - 2 * Janela, bloco.shape[1]))
                # Submete o bloco para processamento
                resultados.append(executor.submit(convolucaoNumba, E, bloco, Janela, opcao, isFirst, isLast))
            # Combinar blocos processados de volta à imagem final
            imagem_processada = np.vstack([f.result() for f in resultados])
            driver = gdal.GetDriverByName('GTiff')
            outDS = driver.Create(outFN.replace(".tif", "") + "_B" + str(band) + ".tif", cols, rows, 1,
                                  gdal.GDT_Float32)
            outDS.SetGeoTransform(Im.GetGeoTransform())
            outDS.SetProjection(Im.GetProjection())
            outBand = outDS.GetRasterBand(1)
            outBand.WriteArray(imagem_processada)
            del outDS
    del Im

if __name__ == '__main__':
    try:
        TI = time.time()
        print('Inicio do programa')
        parser = argparse.ArgumentParser()
        parser.add_argument('--texto1', type=str, required=True, help="Primeiro texto")
        parser.add_argument('--numero1', type=int, required=True, help="Primeiro número inteiro")
        parser.add_argument('--numero2', type=int, required=True, help="Segundo número inteiro")
        parser.add_argument('--texto2', type=str, required=True, help="Segundo texto")
        args = parser.parse_args()
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            processar_imagem_por_blocos(args.texto1, args.numero1, args.numero2, args.texto2, mp.cpu_count())
        TF = time.time()
        print("Tempo de processamento: " + str(TF - TI))
    except Exception as e:
        print("Ocorreu um erro:", e)
        input("Pressione Enter para sair...")