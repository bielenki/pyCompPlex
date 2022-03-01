**CompPlex Tools: Ferramenta para cálculo de medidas de complexidade em imagens
de sensoriamento remoto acoplada em SIG livre**

**1 Introdução**

O conceito de entropia da informação é especialmente interessante quando se
estudam paisagens a partir de imagens de sensores remotos, como satélites e
veículos aéreos não tripulados (VANTs). A variabilidade nos valores de pixels em
uma imagem de sensor remoto representa a diversidade de informações presentes em
uma paisagem e suas unidades e pode servir, por exemplo, para estimar a mudança
na quantidade de informações no sistema causada pela fragmentação.

Os scripts CompPlex Tools foram criados para permitir o cálculo da entropia da
informação (He), variabilidade (He/Hmax) e medidas de López-Ruiz, Mancini e
Calbet (LMC) e Shiner, Davison e Landsberg (SDL). O CompPlex HeROI possibilita o
cálculo dessas medidas para diferentes regiões de interesse (ROIs) selecionadas
em uma imagem de satélite da área de estudo, seguida da comparação da
complexidade de seus padrões, além de possibilitar a geração de assinaturas de
complexidade para cada ROI. O CompPlex Janus possibilita espacializar os
resultados dessas quatro medidas em mapas de complexidade da paisagem enquanto
que o CompPlex Chronos possibilita uma análise multitemporal das métricas, pixel
a pixel, em imagens de diferentes datas.

**2 Cálculos das medidas de complexidade**

Aplicando-se a teoria informacional de Shannon (1949) aos dados de reflectância
em uma banda de uma imagem de SR, e estes sendo representados por sua
discretização em números digitais (DN) à medida que a ocorrência de um
determinado grupo de valores de DN se torna mais provável do que outros valores,
a entropia da imagem decresce. O valor máximo da entropia da imagem neste caso
somente seria atingido quando a ocorrência dos valores de DN na imagem é
equiprovável, não havendo tendência de concentração de probabilidades de um
determinado valor.

Sendo N o número de estados de DN (quantidade de valores de DN sem repetições)
em uma amostra de pixel selecionados na imagem temos que a entropia máxima
admitida para esta amostra é:

Dividindo-se o total de valores de um determinado estado de DN pelo total de
pixel da amostra temos a probabilidade P(DN) de ocorrência deste valor dentro da
amostra. A entropia de Shannon para a amostra é então calculada como:

Outras medidas de complexidade utilizadas baseiam-se na percepção do
desequilíbrio entre os estados de informação.

Segundo apresentam Lopez-Ruiz, Mancini e Calbet (2010) o desequilíbrio D pode
ser mensurado segundo uma distância entre o estado atual do sistema e a condição
de equilíbrio que é calculada como:

Sendo:

Temos que:

Shiner, Davison, and Landsberg Shiner, Davison, and Landsberg Shiner, Davison e Landsberg (1999) já haviam proposto calcular o desequilíbrio
D’ como o complemento do equilíbrio, sendo:

E assim, temos que:

**3 A Ferramenta CompPlex Tools**

A ferramenta CompPlex foi desenvolvida com a linguagem de script Python na forma
de um Plugin para o software livre de SIG QGIS 3. Foram utilizadas, e assim são
requisitos para a sua utilização, as bibliotecas Python GDAL, NUMPY, NUMBA,
RASTERIO e PANDAS. A inclusão da ferramenta ao QGIS se dá por meio de um plugin
criado pela ferramenta Plugin Builder, em que todos os arquivos python
necessários são armazenados num diretório, que deve ser copiado para o
diretórios de plugins do QGIS para o perfil padrão.

Este plugin insere à interface do QGIS uma nova toolbar de onde pode-se acessar
as três ferramentas para cálculo das métricas de complexidade.

![](media/85af18dbaa0b2068b64093d0bcc4ff56.png)

Barra de ferramentas CompPlex Tools

A ferramenta possui três funções principais. Uma para o cálculo das métricas a
partir de regiões de interesse (HeROI), outra para cálculo das métricas para a
imagem completa a partir de Kernels (Janus) e outra para a análise multitemporal
das métricas, pixel a pixel, em imagens de diferentes datas (Chronos).

**3.1 CompPlex HeROI**

O objetivo da ferramenta é apresentar os resultados dos cálculos de complexidade
em uma imagem de sensoriamento remoto para determinadas regiões de interesse
(ROIs). As ROIs são feições do tipo polígono que delimitam as áreas de interesse
para os cálculos. Os resultados são armazenados num arquivo de texto CSV e
incluídos ao projeto QGIS na forma de uma tabela.

Os cálculos são realizados para todas as bandas de uma imagem, no caso desta ser
multibanda, e cada ROI é identificada por meio de um identificador escolhido da
tabela de atributos do plano de informação ROI Layer.

A utilização da ferramenta se dá por meio de uma caixa de diálogo que é acessada
pela ToolBar CompPlex Tools na área de ferramentas do QGIS .

Nesta caixa de diálogo tem-se acesso a três caixas de seleção (CS). As duas
primeiras CS permitem selecionar layers presentes no projeto. A primeira CS
possui um filtro para layers do tipo imagem e a segunda CS um filtro para layers
vetoriais do tipo polígonos. A terceira CS está vinculada a segunda e permite
selecionar um dos campos da tabela de atributos da layer selecionada na segunda
CS.

Um botão abre uma caixa de diálogo para escolher o diretório e nome de arquivo,
com um filtro para arquivos texto do tipo CSV, para salvar o resultado do
processamento.

Depois da seleção dos parâmetros a ferramenta é executada clicando-se no botão
OK, os cálculos são realizados e salvos no arquivo selecionado e uma tabela com
o nome Resultados é adicionada ao projeto ficando visível no painel de layers,
podendo ser aberta com o comando Open Attribute Table.

![Interface gráfica do usuário, Aplicativo Descrição gerada
automaticamente](media/cd3e1d52fb2a7b9e06741be1e2bf7c69.png)

Caixa de diálogo para seleção de parâmetros da ferramenta HeROI

O algoritmo possui duas funções principais. A primeira delas seleciona os pixels
que se sobrepõem a uma feição (polígono) e armazena seus valores num array e,
além disso, já calcula para esta feição os valores de estatísticas descritivas
(count, min, max, mean e std) para os pixels selecionados. A segunda usa o array
de saída da primeira função para calcular os valores de complexidade (He, Hmax,
He/Hmax, N, SDL e LMC) referentes àquela feição.

Para realizar o cálculo para todas as bandas e todas as feições são utilizados
dois loopings encadeados, o primeiro para o número de bandas da imagem e o
segundo para o número de feições da ROI Layer.

Os resultados dos cálculos de cada feição vão sendo armazenados numa estrutura
de dados do tipo tabela da biblioteca Pandas do Python e ao final do looping
esta tabela é convertida e salva para um arquivo texto CSV.

**3.2 CompPlex Janus**

No caso da ferramenta para o cálculo das métricas de complexidade para toda uma
imagem (Janus) e que o resultado seja também uma imagem com os valores da
métrica selecionada o algoritmo funciona como os algoritmos tradicionais de
filtragem. O usuário seleciona o tamanho de uma janela móvel para realizar os
cálculos por meio de uma convolução.

Neste caso a convolução tem o papel de fazer uma filtragem para extração de
informações de interesse na imagem aos quais são aplicados a função da métrica
selecionada. Mais especificamente, o uso desse filtro é feito através de
matrizes denominadas máscaras ou kernels - como são mais conhecidos na prática.

Durante a aplicação da convolução em uma imagem, o kernel vai se deslocando ao
longo da imagem, como uma janela móvel, que vai selecionando os valores de DN
aos quais se aplica a função da métrica selecionada e o resultado deste cálculo
vai formando uma nova imagem com o seu valor ocupando a posição central do
kernel.

Na figura podemos ver a seleção de DNs por meio do Kernel e o valor da função de
complexidade formando uma imagem de saída.

![Entendendo as convoluções](media/7ce94391c3443af7559a3eebdca803d0.jpeg)

Princípio da janela móvel para a convolução

A utilização da ferramenta se dá por meio de uma caixa de diálogo que é acessada
pela ToolBar CompPlex Tools na área de ferramentas do QGIS .

Nesta caixa de diálogo tem-se acesso a três caixas de seleção (CS). A primeira
CS permite selecionar rasters presentes no projeto, esta CS possui um filtro
para layers do tipo imagem. A segunda CS possui tamanhos pré determinados de
janelas móveis que podem ser selecionados para os cálculos. A terceira CS
permite selecionar qual a métrica de complexidade será utilizada para calcular a
nova imagem.

Um botão abre uma caixa de diálogo para escolher o diretório e nome de arquivo,
com um filtro para arquivos texto do tipo tiff, para salvar o resultado do
processamento.

Depois da seleção dos parâmetros a ferramenta é executada clicando-se no botão
OK, os cálculos são realizados e a imagem é salva no arquivo selecionado.

![Interface gráfica do usuário, Aplicativo Descrição gerada
automaticamente](media/46711791667967275966c255016951ec.png)

Caixa de diálogo para seleção de parâmetros da ferramenta Janus

O algoritmo lê a imagem raster e a coverte num array (rows x cols). Inicia então
uma convolução com dois loopings encadeados percorrendo linhas e colunas do
array onde aplica a máscara de acordo com o tamanho da janela móvel selecionada
e chama uma função principal que calcula a métrica de complexidade selecionada
para os DN contidos na mácara e armazena o resultado num novo array na mesma
posição. Ao final da convolução o novo array é convertido numa imagem raster e
salvo no arquivo selecionado.

**3.3 CompPlex Chronos**

A ferramenta Chronos destina-se à uma análise multitemporal das métricas de
complexidade. Para isso é necessário imagens de diferentes épocas para um mesmo
local, com a mesma resolução espacial, o mesmo sistema de referência, para que
cada pixel nas diferentes imagens representem o mesmo local. Além disso, seria
interessante que o sensor também fosse o mesmo a fim de evitar-se problemas
decorrentes de calibração.

A utilização da ferramenta se dá por meio de uma caixa de diálogo que é acessada
pela ToolBar CompPlex Tools na área de ferramentas do QGIS.

Nesta caixa de diálogo tem-se acesso a dois botões que abrem uma caixa de
diálogo para escolher um diretório de dados de entrada e um diretório de dados
de saída. Nesta ferramenta basta selecionar o diretório onde estão contidos os
arquivos rasters do formato tiff, e todos os arquivos com este formatos serão
utilizados como dados de entrada para a análise multitemporal. Da mesma forma
para os arquivos de saída basta-se escolher o diretório e nele serão salvos
quatro arquivos tiff, um para cada métrica (He, He/HMax, SDL e LMC), como
resultados da análise.

![Interface gráfica do usuário, Aplicativo Descrição gerada
automaticamente](media/a1481d9d2a13056e244650ece2f99b04.png)

Caixa de diálogo para seleção de parâmetros da ferramenta Chronos

O algoritmo lê cada imagem raster de formato tiff contido no diretório
selecionado para os dados de entrada e as coverte em arrays (rows x cols) e os
salva numa lista de arrays. É como se tivéssemos um array de N dimensões sendo N
o número de imagens tiff contidas no diretório.

Com três loopings encadeados (rows, cols e arrays nas lista) percorre-se todos
os pixels das imagens, a cada rodada no looping os pixels da posição (row, col)
de todas as imagens são armazenados num vetor, que é repassado a uma função
principal que calcula as quatros métricas de complexidade com base nos valores
de DN contidos nesse vetor. O valor de cada uma delas é armazenado em um novo
array exatamente na posição (row x col). Ao final da convolução os quatro novos
arrays são convertidos em imagens e salvos no diretório de dados de saída
escolhido.

**4 Testando a ferramenta – Exercícios de aplicação**

No diretório que contém os arquivos do plugin foi colocado uma pasta (Examples)
com alguns arquivos para utilização neste exercício.

Convenções usadas no texto para execução dos exercícios práticos:

| DESCRIÇÃO                                                             | PROCEDIMENTO PARA EXECUÇÃO                                                                                                                                              |
|-----------------------------------------------------------------------|-------------------------------------------------------------------------------------------------------------------------------------------------------------------------|
| Aceitar, ligar, ativar ou selecionar algo.                            | Clicar o botão esquerdo (**BE**) ou o direito (**BD**) do mouse para aceitar, ligar, ativar ou selecionar algo. Ex.: **Ativar** (**BE**); **Selecionar** (**BE**); etc. |
| Ação ou recurso, os quais são indicados pelo nome e ícone.            | Clicar (BE) sobre o ícone ![CatalogWindowShow32](media/8e227e1a4c3387af9e7a9a2b438a427e.png) da barra standard                                                          |
| Entrada de nome ou valor pelo teclado, com borda.                     | Digitar nome ou valor, se ele estiver com bordas. Ex.: Nome: HIDRO                                                                                                      |
| Seleção de uma ação ou recurso, se sublinhado**.**                    | Quando a palavra estiver sublinhada, clicar com **BE** para selecioná-la. Exemplo: Gerenciador: Access (**BE**)                                                         |
| Nome da janela (Negrito, Itálico)                                     | Nome da caixa de diálogo aberta**.** Ex.: **Table of Contents**                                                                                                         |
| Nome do campo ou sub-caixa de diálogo dentro de uma janela (Itálico). | Indica sub-caixa de diálogo ou campo dentro de uma janela aberta. Ex: *Construction Tools*                                                                              |

Inicie o QGIS ![Desenho de um círculo Descrição gerada automaticamente com
confiança média](media/79da839be4b83a0ccab44b4a4a1e1590.jpeg)

Iremos acessar o diretório de plugins default do QGIS para colar o diretório do
plugin CompPlex Tools

Acesse o menu Settings

**Settings**

**User Profiles**

**Open Active Profile Folder**

![Interface gráfica do usuário, Aplicativo, Word Descrição gerada
automaticamente](media/70ecafafaa74e167f3ec76f62f36a6ea.png)

Isso irá abrir o Windows explorer na pasta de perfil padrão do QGIS para o
usuário do computador:

![Interface gráfica do usuário, Aplicativo Descrição gerada
automaticamente](media/ee03f8bd972a26986399a7f4fbd9747a.png)

Selecione a pasta **python**

![Interface gráfica do usuário, Aplicativo, Word Descrição gerada
automaticamente](media/adb7273413f91759a8e8a5805fabffa7.png)

E agora selecione a pasta **plugins**. Cole dentro deste diretório a pasta
**complexidade** que contem a ferramenta.

![Interface gráfica do usuário, Aplicativo, Word Descrição gerada
automaticamente](media/db6f28f9bce36532e8eba73e0d59d9ca.png)

Feche e reabra o QGIS para que o novo plugin fique disponível.

Acesse o menu Plugins

**Plugins**

**Manage and Install Plugins...**

![Interface gráfica do usuário, Aplicativo, Word Descrição gerada
automaticamente](media/9970b7f028f4c4ebd3ed8d1b94d4c6a0.png)

Selecione: Installed

Selecione na lista de plugins instalados: CompPlex Tools

![Interface gráfica do usuário, Texto, Aplicativo, chat ou mensagem de texto
Descrição gerada automaticamente](media/4ec9f35e12b9d163ea1c2f1b17105ecc.png)

Feche a caixa de diálogo de *Plugins:* **Close**

Se a barra de ferramentas CompPlex Tools ainda não ficar ativa na interface do
QGIS será necessário habilitá-la.

Acesse o menu View:

**View**

**Toolbars**

CompPlex Tools

![](media/85af18dbaa0b2068b64093d0bcc4ff56.png)

Agora a ferramenta está disponível para utilização.

Inicie um novo projeto no QGIS e salve-o na pasta Example.

Project

New

Project

Save

Selecione a pasta Example

Nome: EX_CompPlex

![Interface gráfica do usuário, Texto, Aplicativo Descrição gerada
automaticamente](media/0405f0230ee89c3e4824b2f6b46189eb.png)

Adicione os layers e raster para o exercício usando o painel Browser

No painel Browser: **(BD) Favorites**

Add a directory

![Interface gráfica do usuário, Aplicativo, Word Descrição gerada
automaticamente](media/8f6d650fde00afe5d59da8528ec0fb04.png)

Navegue até o diretório Examples e adicione-o aos favoritos

![Interface gráfica do usuário, Texto, Aplicativo Descrição gerada
automaticamente](media/1bfe4b053d80fe846f5ddefe51a571c7.png)

Assim você tem acesso aos arquivos dos layers que usaremos a seguir

Para adicioná-los ao projeto clique com BE sobre o layer, mantenha apertado o
BE, arraste o layer para dentro da área de mapas e solte o botão.

Salve novamente o projeto. Agora os layers aparecem no painel de layers.

![Interface gráfica do usuário, Texto, Aplicativo Descrição gerada
automaticamente](media/02a2d22b9f07a57d16b87a6b47c34a05.png)

Iniciaremos com uma análise das métricas de complexidade para regiões de
interesse usando a ferramenta CompPlex HeROI.

Umas das layers adicionadas ao projeto é um vetor de polígonos que delimitam
áreas específicas de quatro classes de uso diferentes que são água, pastagem,
área urbanizada e vegetação nativa. Uma imagem de sensoriamento desta região
será usada para o cálculo das medidas de complexidade para cada tipo de classe.

Na barra de ferramenta ComPlex Tools clique no botão (BE) ![Gráfico, Gráfico de
bolhas Descrição gerada
automaticamente](media/8e4e9e075dd06fa50c090f039d1661c8.png)

A caixa de diálogo da ferramenta **HeROI** é ativada e nela são selecionados os
parâmetros para a sua execução.

Raster Input: SJRP_5metros

ROI: ROIs

ROI Field: Classe

Output File: clique no botão (**BE**)
![](media/4ff2db44af5b35cb6d661ba2ce07d17b.png), navegue até o diretório do
exercício e use o nome analiseROI

![Interface gráfica do usuário, Aplicativo Descrição gerada
automaticamente](media/1163d656281004d41a933f68183e051f.png)

Clique no botão (**BE**) **OK** para executar

O arquivo analiseROI.csv é salvo no diretório escolhido e uma tabela é
adicionada ao painel de layers com o nome de Resultado.

![Interface gráfica do usuário, Texto, Aplicativo Descrição gerada
automaticamente](media/4ede03301a6096df0f6e6d1478d069bc.png)

![Interface gráfica do usuário, Aplicativo, Word Descrição gerada
automaticamente](media/544c6ce82a1904f1f1d0546206725bb0.png)

Para calcular uma imagem de uma das medidas de complexidade usamos a ferramenta
Janus.

Na barra de ferramenta ComPlex Tools clique no botão (BE) ![Desenho de um
cachorro Descrição gerada
automaticamente](media/6f28570a98e4dcb2e10f582823ece513.png)

A caixa de diálogo da ferramenta **Janus** é ativada e nela são selecionados os
parâmetros para a sua execução.

Raster Input: SJRP_5metros_Clip (usaremos aqui apenas um recorte da imagem por
questão de tempo de processamento)

Kernel: 3x3

Type of Metrics: He

Output File: clique no botão (**BE**)
![](media/4ff2db44af5b35cb6d661ba2ce07d17b.png), navegue até o diretório do
exercício e use o nome SJRP_He_3x3

![Interface gráfica do usuário, Aplicativo Descrição gerada
automaticamente](media/b3114507532c456b4564145f75ef7005.png)

Clique no botão (**BE**) **OK** para executar

A imagem calculada é salva no diretório selecionado, note que foi adicionado ao
nome do arquivo um sufixo \_B1 que corresponde a banda 1, caso a imagem de
entrada fosse multibanda teríamos uma imagem de saída para cada banda da imagem
de entrada e identificada pelo respectivo sufixo.

![Interface gráfica do usuário, Texto, Aplicativo Descrição gerada
automaticamente](media/154a768c94e76670ea0245a534dd155d.png)

Podemos adicionar a imagem resultante arrastando-a para dentro da área de mapas.

De forma análoga poderíamos calcular imagens para outras métricas ou com outros
tamanhos de janelas móveis alterando a seleção de parâmetros na caixa de diálogo
nas respectivas caixas de seleção.

Para realizar a análise multitemporal das medidas de complexidade usamos a
ferramenta Chronos.

Na barra de ferramenta ComPlex Tools clique no botão (**BE**) ![Desenho de um
cachorro Descrição gerada
automaticamente](media/5368c0e9e32db424b753b9d2ec4c075c.png)

A caixa de diálogo da ferramenta **Chronos** é ativada e nela são selecionados
os parâmetros para a sua execução.

A entrada neste caso é a indicação de um diretório onde se encontram as imagens
que serão utilizadas na análise. No painel Browser podemos ver a pasta Tiff
dentro do diretório Examples contendo uma série de imagens, e um diretório que
foi criado para receber as imagens resultantes da análise e ainda está vazio.

![Interface gráfica do usuário, Aplicativo Descrição gerada
automaticamente](media/85ea205dbfb5fabeb0b66f98b85f624c.png)

Selecione na caixa de diálogo os parâmetros

Raster Input Directory: Tiff

Raster Output Directory: Tiff/output

![Interface gráfica do usuário, Texto, Aplicativo, chat ou mensagem de texto
Descrição gerada automaticamente](media/9c6badab75edfb5898c0df64ce28a7cb.png)

Clique no botão (BE) **OK** para executar

As imagens calculadas são salvas no diretório selecionado, uma para cada métrica
de complexidade em que são adicionados os respectivos sufixos para
identificação.

![Interface gráfica do usuário, Texto, Aplicativo Descrição gerada
automaticamente](media/bc4710e0adc7517e60444df29318156b.png)
