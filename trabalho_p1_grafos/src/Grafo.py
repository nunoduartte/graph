from collections import deque
import statistics
import psutil
import os
import time


class Grafo:
    def __init__(self, num_vertices, representacao="lista"):
        self.num_vertices = num_vertices
        self.representacao = representacao
        self.num_arestas = 0

        if representacao == "lista":
            self.adjacencia = {i: [] for i in range(1, num_vertices + 1)}
        elif representacao == "matriz":
            self.adjacencia = [[0] * num_vertices for _ in range(num_vertices)]
        else:
            raise ValueError("Representação deve ser 'lista' ou 'matriz'.")

    @staticmethod
    def medir_memoria():
        process = psutil.Process(os.getpid())
        memoria = process.memory_info().rss / 1024 / 1024  # Converter para MB
        return memoria

    @staticmethod
    def comparar_memoria(arquivo):
        # Carregar grafo com lista de adjacência
        print("Carregando grafo usando lista de adjacência...")
        memoria_inicial = Grafo.medir_memoria()
        grafo_lista = Grafo.ler_grafo(arquivo, representacao="lista")
        memoria_lista = Grafo.medir_memoria() - memoria_inicial
        print(f"Memória usada pela lista de adjacência: {memoria_lista:.2f} MB")

        # Carregar grafo com matriz de adjacência
        #print("Carregando grafo usando matriz de adjacência...")
        #memoria_inicial = Grafo.medir_memoria()
        #grafo_matriz = Grafo.ler_grafo(arquivo, representacao="matriz")
        #memoria_matriz = Grafo.medir_memoria() - memoria_inicial
        #print(f"Memória usada pela matriz de adjacência: {memoria_matriz:.2f} MB")

    @staticmethod
    def medir_tempo_bfs(grafo, num_buscas=100):
        inicio_total = time.time()
        num_vertices = grafo.num_vertices

        for i in range(1, num_buscas + 1):
            vertice_inicial = (i % num_vertices) + 1  # Usar diferentes vértices como ponto de partida
            grafo.busca_largura(vertice_inicial)

        tempo_total = time.time() - inicio_total
        tempo_medio = tempo_total / num_buscas
        return tempo_medio

    @staticmethod
    def comparar_tempo_bfs(arquivo):
        # Carregar grafo com lista de adjacência e medir tempo
        print("Carregando grafo usando lista de adjacência...")
        grafo_lista = Grafo.ler_grafo(arquivo, representacao="lista")
        tempo_medio_lista = Grafo.medir_tempo_bfs(grafo_lista)
        print(f"Tempo médio de BFS com lista de adjacência: {tempo_medio_lista:.6f} segundos")

        ## Carregar grafo com matriz de adjacência e medir tempo
        print("Carregando grafo usando matriz de adjacência...")
        grafo_matriz = Grafo.ler_grafo(arquivo, representacao="matriz")
        tempo_medio_matriz = Grafo.medir_tempo_bfs(grafo_matriz)
        print(f"Tempo médio de BFS com matriz de adjacência: {tempo_medio_matriz:.6f} segundos")

    @staticmethod
    def medir_tempo_dfs(grafo, num_buscas=100):
        inicio_total = time.time()
        num_vertices = grafo.num_vertices

        for i in range(1, num_buscas + 1):
            vertice_inicial = (i % num_vertices) + 1  # Usar diferentes vértices como ponto de partida
            grafo.busca_profundidade(vertice_inicial)

        tempo_total = time.time() - inicio_total
        tempo_medio = tempo_total / num_buscas
        return tempo_medio

    @staticmethod
    def comparar_tempo_dfs(arquivo):
        # Carregar grafo com lista de adjacência e medir tempo
        print("Carregando grafo usando lista de adjacência...")
        grafo_lista = Grafo.ler_grafo(arquivo, representacao="lista")
        tempo_medio_lista = Grafo.medir_tempo_dfs(grafo_lista)
        print(f"Tempo médio de DFS com lista de adjacência: {tempo_medio_lista:.6f} segundos")

        # Carregar grafo com matriz de adjacência e medir tempo
        print("Carregando grafo usando matriz de adjacência...")
        grafo_matriz = Grafo.ler_grafo(arquivo, representacao="matriz")
        tempo_medio_matriz = Grafo.medir_tempo_dfs(grafo_matriz)
        print(f"Tempo médio de DFS com matriz de adjacência: {tempo_medio_matriz:.6f} segundos")

    @staticmethod
    def ler_grafo(arquivo, representacao="lista"):
        with open(arquivo, 'r') as f:
            num_vertices = int(f.readline().strip())
            grafo = Grafo(num_vertices, representacao)
            for linha in f:
                v1, v2 = map(int, linha.strip().split())
                grafo.adicionar_aresta(v1, v2)
                grafo.num_arestas += 1
        return grafo

    def adicionar_aresta(self, v1, v2):
        if self.representacao == "lista":
            self.adjacencia[v1].insert(0, v2)
            self.adjacencia[v2].insert(0, v1)
        elif self.representacao == "matriz":
            self.adjacencia[v1 - 1][v2 - 1] = 1
            self.adjacencia[v2 - 1][v1 - 1] = 1

    def exibir_grafo(self):
        if self.representacao == "lista":
            for vertice, vizinhos in self.adjacencia.items():
                print(f'{vertice}: {vizinhos}')
        elif self.representacao == "matriz":
            for i in range(self.num_vertices):
                print(f'{i + 1}: {self.adjacencia[i]}')

    def grau_vertices(self):
        if self.representacao == "lista":
            graus = {v: len(self.adjacencia[v]) for v in self.adjacencia}
        elif self.representacao == "matriz":
            graus = {i + 1: sum(self.adjacencia[i]) for i in range(self.num_vertices)}
        return graus

    def gerar_estatisticas(self):
        graus = list(self.grau_vertices().values())
        grau_minimo = min(graus)
        grau_maximo = max(graus)
        grau_medio = sum(graus) / len(graus)
        grau_mediana = statistics.median(graus)
        return {
            'num_vertices': self.num_vertices,
            'num_arestas': self.num_arestas,
            'grau_minimo': grau_minimo,
            'grau_maximo': grau_maximo,
            'grau_medio': grau_medio,
            'grau_mediana': grau_mediana
        }

    def buscar_vizinhos(self, vertice):
        """ Retorna os vizinhos de um vértice dependendo da representação. """
        if self.representacao == "lista":
            return self.adjacencia[vertice]
        elif self.representacao == "matriz":
            return [i + 1 for i, adj in enumerate(self.adjacencia[vertice - 1]) if adj == 1]

    def busca_largura(self, inicio):
        visitados = {v: False for v in range(1, self.num_vertices + 1)}
        nivel = {v: -1 for v in range(1, self.num_vertices + 1)}
        arvore = {v: None for v in range(1, self.num_vertices + 1)}
        fila = deque([inicio])

        visitados[inicio] = True
        nivel[inicio] = 0

        while fila:
            v = fila.popleft()
            for vizinho in self.buscar_vizinhos(v):  # Vizinhos ajustados
                if not visitados[vizinho]:
                    visitados[vizinho] = True
                    nivel[vizinho] = nivel[v] + 1
                    arvore[vizinho] = v
                    fila.append(vizinho)

        return arvore, nivel

    def busca_largura_para_distancia(self, inicio):
        visitados = {v: False for v in range(1, self.num_vertices + 1)}
        distancias = {v: -1 for v in range(1, self.num_vertices + 1)}
        fila = deque([inicio])

        visitados[inicio] = True
        distancias[inicio] = 0

        while fila:
            v = fila.popleft()
            for vizinho in self.adjacencia[v]:
                if not visitados[vizinho]:
                    visitados[vizinho] = True
                    distancias[vizinho] = distancias[v] + 1
                    fila.append(vizinho)

        return distancias

    def busca_profundidade(self, inicio):
        visitados = {v: False for v in range(1, self.num_vertices + 1)}
        nivel = {v: -1 for v in range(1, self.num_vertices + 1)}
        arvore = {v: None for v in range(1, self.num_vertices + 1)}

        def dfs(v, n):
            visitados[v] = True
            nivel[v] = n
            for vizinho in self.buscar_vizinhos(v):
                if vizinho in visitados and not visitados[vizinho]:  # Verifica se o vizinho é válido e não visitado
                    arvore[vizinho] = v
                    dfs(vizinho, n + 1)

        dfs(inicio, 0)
        return arvore, nivel

    def distancia_entre_vertices(self, inicio, destino):
        if inicio == destino:
            return 0

        visitados = {v: False for v in range(1, self.num_vertices + 1)}
        nivel = {v: -1 for v in range(1, self.num_vertices + 1)}
        fila = deque([inicio])

        visitados[inicio] = True
        nivel[inicio] = 0

        while fila:
            v = fila.popleft()
            for vizinho in self.buscar_vizinhos(v):  # Vizinhos ajustados
                if not visitados[vizinho]:
                    visitados[vizinho] = True
                    nivel[vizinho] = nivel[v] + 1
                    fila.append(vizinho)
                    if vizinho == destino:
                        return nivel[vizinho]

        return -1

    def calcular_diametro_exato(self):
        diametro = 0
        for vertice in range(1, self.num_vertices + 1):
            _, distancias = self.busca_largura(vertice)
            max_distancia = max(distancias.values())
            diametro = max(diametro, max_distancia)
        return diametro

    def calcular_diametro_aproximado(self):
        primeiro_vertice = 1
        _, distancias = self.busca_largura(primeiro_vertice)
        vertice_mais_distante = max(distancias, key=distancias.get)
        _, distancias = self.busca_largura(vertice_mais_distante)
        diametro_aproximado = max(distancias.values())
        return diametro_aproximado

    def calcular_diametro(self, modo_diametro="exato"):
        if modo_diametro == "exato":
            return self.calcular_diametro_exato()
        elif modo_diametro == "aproximado":
            return self.calcular_diametro_aproximado()
        else:
            raise ValueError("Modo de diâmetro deve ser 'exato' ou 'aproximado'.")


    def componentes_conexas(self):
        visitados = {v: False for v in range(1, self.num_vertices + 1)}
        componentes = []

        for vertice in range(1, self.num_vertices + 1):
            if not visitados[vertice]:
                componente_atual = []
                fila = deque([vertice])
                visitados[vertice] = True

                while fila:
                    v = fila.popleft()
                    componente_atual.append(v)
                    for vizinho in self.buscar_vizinhos(v):  # Vizinhos ajustados
                        if not visitados[vizinho]:
                            visitados[vizinho] = True
                            fila.append(vizinho)

                componentes.append(componente_atual)

        componentes.sort(key=len, reverse=True)
        return componentes

    def salvar_arvore_busca(self, inicio, tipo_busca, arquivo_saida):
        if tipo_busca == 'bfs':
            arvore, niveis = self.busca_largura(inicio)
        elif tipo_busca == 'dfs':
            arvore, niveis = self.busca_profundidade(inicio)
        else:
            raise ValueError("Tipo de busca deve ser 'bfs' ou 'dfs'.")

        with open(arquivo_saida, 'w') as f:
            f.write(f"Árvore de Busca ({tipo_busca.upper()}) a partir do vértice {inicio}:\n")
            f.write("Vértice\tPai\tNível\n")
            for vertice in range(1, self.num_vertices + 1):
                f.write(f"{vertice}\t{arvore[vertice]}\t{niveis[vertice]}\n")
        print(f"Resultados da busca ({tipo_busca.upper()}) foram salvos em: {arquivo_saida}")

    def salvar_estatisticas(self, arquivo_saida):
        estatisticas = self.gerar_estatisticas()
        componentes = self.componentes_conexas()
        with open(arquivo_saida, 'w') as f:
            f.write(f"Número de vértices: {estatisticas['num_vertices']}\n")
            f.write(f"Número de arestas: {estatisticas['num_arestas']}\n")
            f.write(f"Grau mínimo: {estatisticas['grau_minimo']}\n")
            f.write(f"Grau máximo: {estatisticas['grau_maximo']}\n")
            f.write(f"Grau médio: {estatisticas['grau_medio']:.2f}\n")
            f.write(f"Grau mediana: {estatisticas['grau_mediana']}\n")
            f.write(f"Componentes conexas: {len(componentes)}\n")
            for i, componente in enumerate(componentes, 1):
                f.write(f"Componente {i}: {componente} (Tamanho: {len(componente)})\n")

    ##### Estudo de Caso 4

    def buscar_pais_dos_vertices(self, vertices_iniciais, vertices_a_verificar):
        resultados = {}

        for vertice_inicial in vertices_iniciais:
            # Executar BFS
            arvore_bfs, _ = self.busca_largura(vertice_inicial)
            # Executar DFS
            arvore_dfs, _ = self.busca_profundidade(vertice_inicial)

            # Armazenar os pais dos vértices solicitados
            resultados[vertice_inicial] = {
                'BFS': {v: arvore_bfs[v] for v in vertices_a_verificar},
                'DFS': {v: arvore_dfs[v] for v in vertices_a_verificar}
            }

        return resultados

    def salvar_resultados_pais(self, resultados, vertices_iniciais, vertices_a_verificar, arquivo_saida):
        with open(arquivo_saida, 'w') as f:
            for vertice_inicial in vertices_iniciais:
                f.write(f"Resultados para busca a partir do vértice {vertice_inicial}:\n")
                f.write("Método\tVértice 10\tVértice 20\tVértice 30\n")
                f.write("BFS\t")
                f.write(f"{resultados[vertice_inicial]['BFS'].get(10, 'N/A')}\t")
                f.write(f"{resultados[vertice_inicial]['BFS'].get(20, 'N/A')}\t")
                f.write(f"{resultados[vertice_inicial]['BFS'].get(30, 'N/A')}\n")
                f.write("DFS\t")
                f.write(f"{resultados[vertice_inicial]['DFS'].get(10, 'N/A')}\t")
                f.write(f"{resultados[vertice_inicial]['DFS'].get(20, 'N/A')}\t")
                f.write(f"{resultados[vertice_inicial]['DFS'].get(30, 'N/A')}\n")
                f.write("\n")
        print(f"Resultados salvos no arquivo: {arquivo_saida}")

    #### Estudo de Caso 5
    def calcular_distancias_pares(self, pares_vertices):
        resultados = {}
        for v1, v2 in pares_vertices:
            dist = self.distancia_entre_vertices(v1, v2)
            resultados[(v1, v2)] = dist
        return resultados

    def salvar_resultados_distancias(self, resultados, arquivo_saida):
        with open(arquivo_saida, 'w') as f:
            f.write("Distâncias entre pares de vértices:\n")
            for (v1, v2), dist in resultados.items():
                f.write(f"Distância entre {v1} e {v2}: {dist}\n")
        print(f"Resultados de distância salvos no arquivo: {arquivo_saida}")

    #### Estudo de Caso 6
    def analisar_componentes_conexas(self):
        componentes = self.componentes_conexas()
        num_componentes = len(componentes)
        maior_componente = max(componentes, key=len)
        menor_componente = min(componentes, key=len)

        return {
            "num_componentes": num_componentes,
            "tamanho_maior_componente": len(maior_componente),
            "tamanho_menor_componente": len(menor_componente),
            "componentes": componentes
        }

    def salvar_componentes_conexas(self, resultado, arquivo_saida):
        with open(arquivo_saida, 'w') as f:
            f.write(f"Número de componentes conexas: {resultado['num_componentes']}\n")
            f.write(f"Tamanho da maior componente: {resultado['tamanho_maior_componente']}\n")
            f.write(f"Tamanho da menor componente: {resultado['tamanho_menor_componente']}\n")
            f.write("\nComponentes:\n")
            for i, componente in enumerate(resultado['componentes'], 1):
                f.write(f"Componente {i}: {componente} (Tamanho: {len(componente)})\n")
        print(f"Resultados de componentes conexas salvos no arquivo: {arquivo_saida}")

    ####Estudo de Caso 7

    def salvar_diametro(self, diametro, arquivo_saida):
        with open(arquivo_saida, 'w') as f:
            f.write(f"Diâmetro do grafo: {diametro}\n")
        print(f"Resultado do diâmetro salvo no arquivo: {arquivo_saida}")