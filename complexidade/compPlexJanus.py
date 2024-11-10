import os, sys
import numpy as np
from numba import njit, prange
import multiprocessing as mp
import time
from osgeo import gdal
import tkinter as tk
from tkinter import ttk
import argparse
import warnings
from numba.core.errors import NumbaWarning, NumbaPendingDeprecationWarning
warnings.filterwarnings("ignore", category=NumbaWarning)
warnings.filterwarnings("ignore", category=NumbaPendingDeprecationWarning)

def get_resource_path(relative_path):
    """Obtém o caminho absoluto para os arquivos, quer esteja em modo de desenvolvimento ou empacotado no executável."""
    try:
        base_path = sys._MEIPASS
    except Exception:
        base_path = os.path.abspath(".")
    return os.path.join(base_path, relative_path)

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

@njit #(parallel=True)
def convolNumba(ImArray, rows, cols, janela):
    arrayHe = np.empty((rows, cols), dtype=np.float64)
    arrayHemax = np.empty((rows, cols), dtype=np.float64)
    arraySDL = np.empty((rows, cols), dtype=np.float64)
    arrayLMC = np.empty((rows, cols), dtype=np.float64)
    for row in prange(rows):
        for col in prange(cols):
            Lx = max(0, col - janela)
            Ux = min(cols, col + janela +1)
            Ly = max(0, row - janela)
            Uy = min(rows, row + janela +1)
            mascara = ImArray[Ly:Uy, Lx:Ux].ravel()
            lenVet = mascara.size
            Lista = list(set(mascara))
            if len(Lista) == 1 and Lista.count(0) == 1:
                arrayHe[row, col] = 0
                arrayHemax[row, col] = 0
                arrayLMC[row, col] = 0
                arraySDL[row, col] = 0
            else:
                He, prob = calcular_He(mascara, Lista, lenVet)
                arrayHe[row, col] = He
                N = len(Lista) * 1.0
                Hmax = np.log2(N) if N > 1 else 0
                C = He / Hmax if Hmax != 0 else 0
                arrayHemax[row, col] = C
                SDL = (1 - C) * C
                arraySDL[row, col] = SDL
                D = np.sum((prob - (1 / N)) ** 2) if N > 1 else 0
                LMC = D * C
                arrayLMC[row, col] = LMC
    return arrayHe, arrayHemax, arraySDL, arrayLMC

def salvar_geotiff_com_georef(nome_arquivo_saida, array, geotransform, projecao, tipo_dado):
    driver = gdal.GetDriverByName('GTiff')
    linhas, colunas = array.shape
    dataset_saida = driver.Create(nome_arquivo_saida, colunas, linhas, 1, tipo_dado)
    if dataset_saida is None:
        raise ValueError(f"Erro ao criar o arquivo GeoTIFF: {nome_arquivo_saida}")
    dataset_saida.SetGeoTransform(geotransform)
    dataset_saida.SetProjection(projecao)
    banda = dataset_saida.GetRasterBand(1)
    banda.WriteArray(array)
    dataset_saida.FlushCache()
    dataset_saida = None

def processar_imagem_geotiff(nome_arquivo, Saida, janela, dirSaida, queue):
    inicio_imagem = time.time()
    dataset = gdal.Open(nome_arquivo)
    if dataset is None:
        print(f'Erro ao abrir o arquivo {nome_arquivo}')
        return
    num_bandas = dataset.RasterCount
    geotransform = dataset.GetGeoTransform()
    projecao = dataset.GetProjection()
    for banda_idx in range(1, num_bandas + 1):
        banda = dataset.GetRasterBand(banda_idx)
        array = banda.ReadAsArray()
        linhas, colunas = array.shape
        he, hmax, sdl, lmc = convolNumba(array, linhas, colunas, janela)
        sufixo_banda = f"B{banda_idx}"
        salvar_geotiff_com_georef(
            os.path.join(dirSaida, "He", f"He_{os.path.basename(nome_arquivo).replace('.tif', '')}_{sufixo_banda}.tif"),
            he, geotransform, projecao, gdal.GDT_Float32)
        salvar_geotiff_com_georef(
            os.path.join(dirSaida, "Hmax", f"Hmax_{os.path.basename(nome_arquivo).replace('.tif', '')}_{sufixo_banda}.tif"),
            hmax, geotransform, projecao, gdal.GDT_Float32)
        salvar_geotiff_com_georef(
            os.path.join(dirSaida, "sdl", f"SDL_{os.path.basename(nome_arquivo).replace('.tif', '')}_{sufixo_banda}.tif"),
            sdl, geotransform, projecao, gdal.GDT_Float32)
        salvar_geotiff_com_georef(
            os.path.join(dirSaida, "lmc", f"LMC_{os.path.basename(nome_arquivo).replace('.tif', '')}_{sufixo_banda}.tif"),
            lmc, geotransform, projecao, gdal.GDT_Float32)
        fim_imagem = time.time()
        print(f'arquivo {nome_arquivo} processado em '+str(fim_imagem-inicio_imagem) + "seg")
    queue.put(1)

def atualizar_progresso(progress_bar, total, queue, root, p):
    """Função para verificar a fila e atualizar a barra de progresso"""
    processed = progress_bar['value'] // (100 / total)
    while not queue.empty():  # Enquanto houver elementos na fila, eles serão processados
        queue.get()  # Retira um item da fila
        processed += 1
        progress_bar['value'] = (processed / total) * 100  # Atualiza a barra de progresso
        progress_bar.update_idletasks()
    if not p.is_alive():  # Verifica se o processo foi finalizado
        progress_bar['value'] = 100  # Força a barra de progresso a 100% ao finalizar
        progress_bar.update_idletasks()
        root.after(500, root.destroy)  # Fecha a janela após 500ms
    else:
        root.after(100, lambda: atualizar_progresso(progress_bar, total, queue, root,
                                                    p))  # Continua verificando a cada 100ms

def iniciar_processamento_com_progresso(args, total_imagens):
    root = tk.Tk()
    root.title("Progresso do processamento")
    progress_bar = ttk.Progressbar(root, orient='horizontal', mode='determinate', length=400)
    progress_bar.pack(pady=20)
    queue = mp.Manager().Queue()
    p = mp.Process(target=worker_com_progresso, args=(args, queue))
    p.start()
    root.after(100, lambda: atualizar_progresso(progress_bar, total_imagens, queue, root, p))
    root.mainloop()  # O mainloop deve terminar após root.destroy() ser chamado
    p.join()  # Espera o processo terminar após o fechamento da janela

def worker_com_progresso(args, queue):
    num_processos = mp.cpu_count()
    with mp.Pool(processes=num_processos) as pool:
        pool.starmap(processar_imagem_geotiff, [(arg[0], arg[1], arg[2], arg[3], queue) for arg in args])

def inicioJanus(diretorio, diretorioSaida, janela):
    # Função auxiliar para criar diretórios de saída, se não existirem
    def criar_diretorios_saida(diretorioSaida, pastas=["he", "hmax", "sdl", "lmc"]):
        for pasta in pastas:
            caminho_completo = os.path.join(diretorioSaida, pasta)
            if not os.path.exists(caminho_completo):
                os.makedirs(caminho_completo)
                print(f"Diretório criado: {caminho_completo}")
            else:
                print(f"Diretório já existe: {caminho_completo}")
    criar_diretorios_saida(diretorioSaida)
    arrays_banda_r = []
    nomes_arquivos = []
    for arquivo in os.listdir(diretorio):
        if arquivo.endswith('.tiff') or arquivo.endswith('.tif'):
            caminho_arquivo = os.path.join(diretorio, arquivo)
            print(caminho_arquivo)
            arrays_banda_r.append(caminho_arquivo)
            nomes_arquivos.append(arquivo)
    if not arrays_banda_r:
        print("Nenhum arquivo TIFF encontrado no diretório de entrada.")
        return
    argumentos = [(arrays_banda_r[i], nomes_arquivos[i], janela, diretorioSaida) for i in range(len(arrays_banda_r))]
    input("Pressione Enter para continuar...")
    print(argumentos)
    input("Pressione Enter para continuar...")
    iniciar_processamento_com_progresso(argumentos, len(nomes_arquivos))

if __name__ == '__main__':
    print('inicio subprpograma')
    parser = argparse.ArgumentParser()
    parser.add_argument('--texto1', type=str, required=True, help="Primeiro texto")
    parser.add_argument('--texto2', type=str, required=True, help="Segundo texto")
    parser.add_argument('--numero', type=int, required=True, help="Número inteiro")
    args = parser.parse_args()
    print('antes da função')
    inicioJanus(args.texto1, args.texto2, args.numero)
    print('depois da função')
    input("Pressione Enter para continuar...")