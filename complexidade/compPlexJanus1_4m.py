import argparse
import numpy as np
from numba import prange, njit
from osgeo import gdal
import multiprocessing as mp
from concurrent.futures import ProcessPoolExecutor
import time
import warnings
from numba.core.errors import NumbaWarning, NumbaPendingDeprecationWarning
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
def convolucaoNumba(bloco, janela, isFirst, isLast):
    rows, cols = bloco.shape
    # Ajuste do tamanho dos arrays com base na posição do bloco (primeiro, último ou intermediário)
    if isFirst:
        array_shape = (rows - janela, cols)
        rowInicio = 0
        rowFinal = rows - janela
    elif isLast:
        array_shape = (rows - janela, cols)
        rowInicio = janela
        rowFinal = rows
    else:
        array_shape = (rows - 2 * janela, cols)
        rowInicio = janela
        rowFinal = rows - janela
    # Criação dos arrays de métricas internos
    arrayHe = np.empty(array_shape, dtype=np.float32)
    arrayHmax = np.empty(array_shape, dtype=np.float32)
    arraySDL = np.empty(array_shape, dtype=np.float32)
    arrayLMC = np.empty(array_shape, dtype=np.float32)
    for row in prange(rowInicio, rowFinal):
        for col in prange(cols):
            Lx = max(0, col - janela)
            Ux = min(cols, col + janela + 1)
            Ly = max(0, row - janela)
            Uy = min(rows, row + janela + 1)
            mascara = bloco[Ly:Uy, Lx:Ux].ravel()
            lenVet = mascara.size
            Lista = list(set(mascara))
            if len(Lista) == 1 and Lista.count(0) == 1:
                if isFirst:
                    arrayHe[row , col] = 0
                    arrayHmax[row , col] = 0
                    arraySDL[row , col] = 0
                    arrayLMC[row , col] = 0
                else:
                    arrayHe[row - janela, col] = 0
                    arrayHmax[row - janela, col] = 0
                    arraySDL[row - janela, col] = 0
                    arrayLMC[row - janela, col] = 0
            else:
                He, prob = calcular_He(mascara, Lista, lenVet)
                if isFirst:
                    arrayHe[row , col] = He
                else:
                    arrayHe[row - janela, col] = He
                N = len(Lista) * 1.0
                Hmax = np.log2(N) if N > 1 else 0
                C = He / Hmax if Hmax > 0 else 0
                if isFirst:
                    arrayHmax[row, col] = C
                else:
                    arrayHmax[row - janela, col] = C
                SDL = (1 - C) * C
                if isFirst:
                    arraySDL[row, col] = SDL
                else:
                    arraySDL[row - janela, col] = SDL
                D = np.sum((prob - (1 / N)) ** 2) if N > 1 else 0
                LMC = D * C
                if isFirst:
                    arrayLMC[row, col] = LMC
                else:
                    arrayLMC[row - janela, col] = LMC
    return arrayHe, arrayHmax, arraySDL, arrayLMC


def divide_em_blocos(imagem, num_blocos, janela):
    blocos = []
    bloco_tamanho = imagem.shape[0] // num_blocos
    for i in range(0, imagem.shape[0], bloco_tamanho):
        bloco = imagem[max(0, i - janela):min(imagem.shape[0], i + bloco_tamanho + janela), :]
        blocos.append(bloco)
    return blocos

def processar_imagem_por_blocos(NameImagem, Janela, outFN, num_blocos=4):
    Im = gdal.Open(NameImagem)
    cols = Im.RasterXSize
    rows = Im.RasterYSize
    NrBandas = Im.RasterCount
    with ProcessPoolExecutor() as executor:
        for band in range(1, NrBandas + 1):
            print("Processando banda " + str(band) + "...")
            banda_img = Im.GetRasterBand(band)
            ImArray = banda_img.ReadAsArray().astype(np.float64)
            blocos = divide_em_blocos(ImArray, num_blocos, Janela)
            resultados = []
            for idx, bloco in enumerate(blocos):
                # Verificar se o bloco é o primeiro ou o último
                isFirst = (idx == 0)
                isLast = (idx == len(blocos) - 1)
                resultados.append(executor.submit(convolucaoNumba, bloco, Janela, isFirst, isLast))
                # Inicializar listas para combinar os resultados de cada métrica
            combined_He = []
            combined_Hmax = []
            combined_SDL = []
            combined_LMC = []
            # Processar os resultados dos blocos
            for resultado in resultados:
                arrayHe, arrayHmax, arraySDL, arrayLMC = resultado.result()
                combined_He.append(arrayHe)
                combined_Hmax.append(arrayHmax)
                combined_SDL.append(arraySDL)
                combined_LMC.append(arrayLMC)
            # Gravar cada métrica em arquivos separados
            for metrica, arrays in zip(['He', 'Hmax', 'SDL', 'LMC'],
                                           [combined_He, combined_Hmax, combined_SDL, combined_LMC]):
                imagem_final = np.vstack(arrays)
                driver = gdal.GetDriverByName('GTiff')
                outDS = driver.Create(outFN.replace(".tif", f"_{metrica}_B{band}.tif"), cols, rows, 1, gdal.GDT_Float32)
                outDS.SetGeoTransform(Im.GetGeoTransform())
                outDS.SetProjection(Im.GetProjection())
                outBand = outDS.GetRasterBand(1)
                outBand.WriteArray(imagem_final)
                del outDS
    del Im

if __name__ == '__main__':
    try:
        TI = time.time()
        print('inicio subprograma')
        parser = argparse.ArgumentParser()
        parser.add_argument('--texto1', type=str, required=True, help="Primeiro texto")
        parser.add_argument('--numero1', type=int, required=True, help="Primeiro Número inteiro")
        parser.add_argument('--texto2', type=str, required=True, help="Segundo texto")
        args = parser.parse_args()
        processar_imagem_por_blocos(args.texto1, args.numero1, args.texto2, mp.cpu_count())
        TF = time.time()
        print("Tempo de processamento: " + str(TF - TI))
    except Exception as e:
        print("Ocorreu um erro:", e)
        input("Pressione Enter para sair...")