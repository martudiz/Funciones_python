"""
Funciones_guia1.py

Archivo autocontenido (sin imports locales) que resuelve TODOS los ejercicios de la
guía 1 a partir de una única entrada: un mensaje emitido por la fuente.

De ese mensaje se derivan el alfabeto y las probabilidades (para los ejercicios de
fuente de memoria nula) y la matriz de transición (para los de Markov), y con eso se
resuelven los ejercicios 1, 2, 8, 10, 11, 14 y 15.

Parámetros secundarios (no forman parte del mensaje, pero hacen falta para algunos
ejercicios) quedan como constantes editables al principio de main(): longitud de las
cadenas a simular, orden de la extensión y tolerancia para el vector estacionario.
"""
from math import log
import random


# =============================================================================
# FUENTES DE MEMORIA NULA
# =============================================================================

def informacion_bits(probabilidades, base=2):
    """
    Pre: 'probabilidades' es una lista de floats en (0, 1] con la distribución
         de probabilidad de los símbolos de una fuente de memoria nula.
         'base' es la base del logaritmo: 2 para bits (valor por defecto), o el
         tamaño del alfabeto código (r) si se quiere expresar la información en
         unidades r-arias (como en los ejercicios de códigos).
    Post: retorna una lista paralela con la cantidad de información de cada
          símbolo (I = log_base(1/p)). No modifica la lista recibida.
    """
    return [log(1 / p, base) for p in probabilidades]


def entropia(probabilidades, informacion=None):
    """
    Pre: 'probabilidades' es una lista de floats que suman 1 (o prácticamente 1).
         'informacion' es opcional: si no se pasa, se calcula internamente con
         informacion_bits().
    Post: retorna la entropía de la fuente en bits/símbolo (H = sum(p * I)).
    """
    if informacion is None:
        informacion = informacion_bits(probabilidades)
    return sum(p * i for p, i in zip(probabilidades, informacion))


def alfabeto_y_frecuencias(lista_simbolos):
    """
    Pre: 'lista_simbolos' es una lista de strings, los símbolos emitidos por una fuente.
    Post: retorna dos listas paralelas: el alfabeto (símbolos únicos, en orden de
          primera aparición) y la frecuencia absoluta de cada símbolo.
    """
    alfabeto = []
    frecuencias = []
    for simbolo in lista_simbolos:
        if simbolo in alfabeto:
            frecuencias[alfabeto.index(simbolo)] += 1
        else:
            alfabeto.append(simbolo)
            frecuencias.append(1)
    return alfabeto, frecuencias


def frecuencias_a_probabilidades(frecuencias, total):
    """
    Pre: 'frecuencias' es una lista de enteros >= 0; 'total' es la cantidad total
         de símbolos observados (> 0).
    Post: retorna la lista de probabilidades (frecuencia / total), paralela a
          'frecuencias'.
    """
    return [f / total for f in frecuencias]


def frecuencia_acumulada(probabilidades):
    """
    Pre: 'probabilidades' es una lista de floats que suman 1.
    Post: retorna la lista de frecuencias acumuladas (F(i) = p1 + ... + pi),
          usada para simular por el método de la transformada inversa.
    """
    acumulada = []
    acumulado = 0
    for p in probabilidades:
        acumulado += p
        acumulada.append(acumulado)
    return acumulada


def elegir_simbolo(alfabeto, fiacum):
    """
    Pre: 'alfabeto' es una lista de símbolos; 'fiacum' es su frecuencia acumulada
         (misma longitud que 'alfabeto', último valor ≈ 1).
    Post: retorna un símbolo elegido al azar, respetando la distribución de
          probabilidad representada por 'fiacum'.
    """
    r = random.random()
    x = 0
    while x < len(fiacum) - 1 and r > fiacum[x]:
        x += 1
    return alfabeto[x]


def generar_cadena(alfabeto, probabilidades, n):
    """
    Pre: 'alfabeto' y 'probabilidades' son listas paralelas de una fuente de
         memoria nula; 'n' es la longitud deseada de la cadena (entero > 0).
    Post: retorna una lista de n símbolos generados simulando la fuente.
    """
    fiacum = frecuencia_acumulada(probabilidades)
    return [elegir_simbolo(alfabeto, fiacum) for _ in range(n)]


def extension_orden_n(alfabeto, probabilidades, n):
    """
    Pre: 'alfabeto' y 'probabilidades' son listas paralelas de una fuente de
         memoria nula de orden 1; 'n' es un entero >= 1.
    Post: retorna (alfabeto_extendido, probabilidades_extendidas), la extensión
          de orden n de la fuente (S^n), calculada recursivamente.
    """
    if n == 1:
        return list(alfabeto), list(probabilidades)

    alfabeto_anterior, prob_anterior = extension_orden_n(alfabeto, probabilidades, n - 1)

    nuevo_alfabeto = []
    nueva_prob = []
    for i in range(len(alfabeto_anterior)):
        for j in range(len(alfabeto)):
            nuevo_alfabeto.append(alfabeto_anterior[i] + alfabeto[j])
            nueva_prob.append(prob_anterior[i] * probabilidades[j])
    return nuevo_alfabeto, nueva_prob


# =============================================================================
# FUENTES DE MARKOV
# =============================================================================

def vector_estacionario(matriz, tolerancia=0.001, max_iteraciones=1000):
    """
    Pre: 'matriz' es una matriz de transición cuadrada (lista de listas) donde
         cada fila suma 1. 'tolerancia' es la diferencia máxima admitida entre
         iteraciones sucesivas para considerar que convergió.
    Post: retorna el vector estacionario (lista de floats) de la cadena de
          Markov representada por 'matriz', calculado por iteración de potencias.
          Si no converge en 'max_iteraciones', retorna la última aproximación.
    """
    n = len(matriz)
    vector_actual = [1 / n] * n
    for _ in range(max_iteraciones):
        vector_nuevo = [
            sum(vector_actual[i] * matriz[i][j] for i in range(n)) for j in range(n)
        ]
        if all(abs(vector_nuevo[k] - vector_actual[k]) <= tolerancia for k in range(n)):
            return vector_nuevo
        vector_actual = vector_nuevo
    return vector_actual


def entropia_markov(matriz, tolerancia=0.001):
    """
    Pre: 'matriz' es una matriz de transición válida (cada fila suma 1).
    Post: retorna la entropía de la fuente de Markov, ponderando la entropía de
          cada fila (estado) por su probabilidad estacionaria.
    """
    v_est = vector_estacionario(matriz, tolerancia)
    n = len(matriz)
    total = 0
    for i in range(n):
        fila_entropia = sum(
            matriz[i][j] * log(1 / matriz[i][j], 2) for j in range(n) if matriz[i][j] > 0
        )
        total += v_est[i] * fila_entropia
    return total


def matriz_transicion_desde_mensaje(lista_simbolos, alfabeto):
    """
    Pre: 'lista_simbolos' es la secuencia de símbolos emitidos por la fuente;
         'alfabeto' es la lista de símbolos únicos (por ejemplo, obtenida con
         alfabeto_y_frecuencias).
    Post: retorna la matriz de transición (lista de listas) estimada contando
          pares consecutivos en 'lista_simbolos' y normalizando cada fila para
          que sume 1. Si un símbolo nunca aparece como origen, su fila queda en 0.
    """
    n = len(alfabeto)
    matriz = [[0] * n for _ in range(n)]
    for i in range(len(lista_simbolos) - 1):
        fila = alfabeto.index(lista_simbolos[i])
        columna = alfabeto.index(lista_simbolos[i + 1])
        matriz[fila][columna] += 1
    for fila in matriz:
        total_fila = sum(fila)
        if total_fila > 0:
            for j in range(n):
                fila[j] /= total_fila
    return matriz


def generar_cadena_markov(alfabeto, matriz, n):
    """
    Pre: 'alfabeto' y 'matriz' representan una fuente de Markov (matriz cuadrada,
         cada fila suma 1); 'n' es la longitud deseada de la cadena (>= 1).
    Post: retorna una lista de n símbolos generados simulando la fuente. El
          primer símbolo se elige según el vector estacionario; los siguientes,
          según la fila de la matriz correspondiente al símbolo anterior.
    """
    v_est = vector_estacionario(matriz)
    fiacum_est = frecuencia_acumulada(v_est)
    simbolo_actual = elegir_simbolo(alfabeto, fiacum_est)
    cadena = [simbolo_actual]
    for _ in range(n - 1):
        fila = matriz[alfabeto.index(simbolo_actual)]
        fiacum_fila = frecuencia_acumulada(fila)
        simbolo_actual = elegir_simbolo(alfabeto, fiacum_fila)
        cadena.append(simbolo_actual)
    return cadena


def es_memoria_nula(matriz, tolerancia):
    """
    Pre: 'matriz' es una matriz de transición cuadrada (lista de listas);
         'tolerancia' es la diferencia máxima admitida entre filas para
         considerarlas "iguales".
    Post: retorna True si todas las filas de 'matriz' son iguales entre sí
          (dentro de la tolerancia) -> fuente de memoria nula.
          Retorna False si al menos una fila difiere -> fuente con memoria.
    """
    n = len(matriz)
    if n <= 1:
        return True
    fila_ref = matriz[0]
    for i in range(1, n):
        for j in range(n):
            if abs(matriz[i][j] - fila_ref[j]) > tolerancia:
                return False
    return True


# =============================================================================
# MAIN: resuelve todos los ejercicios a partir de un único mensaje
# =============================================================================

def main():
    # ------------------------------------------------------------------
    # --- ÚNICA ENTRADA: el mensaje emitido por la fuente ---
    mensaje = "A B C C A B B A C A".split()

    # --- Parámetros secundarios (editables) ---
    longitud_simulacion = 10   # longitud de las cadenas a generar (ej. 2b y 15b)
    orden_extension = 3        # orden N de la extensión (ej. 10)
    tolerancia = 0.01          # tolerancia para vector estacionario / memoria nula (ej. 14 y 15c)
    # ------------------------------------------------------------------

    print("=" * 70)
    print("MENSAJE DE ENTRADA:", " ".join(mensaje))
    print("=" * 70)

    # --- Ejercicio 2a: alfabeto y probabilidades a partir del mensaje ---
    alfabeto, frecuencias = alfabeto_y_frecuencias(mensaje)
    probabilidades = frecuencias_a_probabilidades(frecuencias, len(mensaje))
    print("\n--- Ejercicio 2a: alfabeto y probabilidades ---")
    print("Alfabeto:", alfabeto)
    print("Probabilidades:", probabilidades)

    # --- Ejercicio 1: información en bits y entropía ---
    print("\n--- Ejercicio 1: información en bits y entropía ---")
    info = informacion_bits(probabilidades)
    print("Información en bits:", info)
    h = entropia(probabilidades, info)
    print("Entropía H(S):", h)

    # --- Ejercicio 2b: simulación de una cadena de memoria nula ---
    print(f"\n--- Ejercicio 2b: simulación de una cadena de longitud {longitud_simulacion} ---")
    cadena = generar_cadena(alfabeto, probabilidades, longitud_simulacion)
    print("Cadena generada:", cadena)

    # --- Ejercicio 8: fórmula binaria (solo tiene sentido si el alfabeto es binario) ---
    print("\n--- Ejercicio 8: entropía por fórmula binaria ---")
    if len(alfabeto) == 2:
        w = probabilidades[0]
        print(f"w = {w} -> entropía = {entropia([w, 1 - w])} (debería coincidir con H(S) de arriba)")
    else:
        print(f"No aplica: la fórmula del ej. 8 es para fuentes binarias y el alfabeto tiene {len(alfabeto)} símbolos.")

    # --- Ejercicio 10: extensión de orden N ---
    print(f"\n--- Ejercicio 10: extensión de orden {orden_extension} ---")
    alfabeto_n, prob_n = extension_orden_n(alfabeto, probabilidades, orden_extension)
    print(f"Alfabeto de orden {orden_extension}:", alfabeto_n)
    print(f"Probabilidades de orden {orden_extension}:", prob_n)

    # --- Ejercicio 11: extensiones de orden 2 y 3, verificación H(Sn) = n * H(S) ---
    print("\n--- Ejercicio 11: extensiones de orden 2 y 3 ---")
    for n in (2, 3):
        alfabeto_n2, prob_n2 = extension_orden_n(alfabeto, probabilidades, n)
        h_n = entropia(prob_n2)
        print(f"Orden {n}: H(S_{n}) = {h_n:.4f} bits | n * H(S) = {n * h:.4f} bits (verificación)")

    # --- Ejercicio 15a: matriz de transición a partir del mismo mensaje ---
    print("\n--- Ejercicio 15a: matriz de transición a partir del mensaje ---")
    matriz = matriz_transicion_desde_mensaje(mensaje, alfabeto)
    for fila in matriz:
        print([round(v, 3) for v in fila])

    # --- Ejercicio 14: vector estacionario y entropía de Markov ---
    print("\n--- Ejercicio 14: vector estacionario y entropía de Markov ---")
    v_est = vector_estacionario(matriz, tolerancia=tolerancia)
    print("Vector estacionario:", [round(v, 4) for v in v_est])
    h_markov = entropia_markov(matriz, tolerancia=tolerancia)
    print("Entropía de Markov:", h_markov)

    # --- Ejercicio 15b: simulación de una cadena de Markov ---
    print(f"\n--- Ejercicio 15b: simulación de una cadena de Markov de longitud {longitud_simulacion} ---")
    cadena_markov = generar_cadena_markov(alfabeto, matriz, longitud_simulacion)
    print("Cadena simulada:", " ".join(cadena_markov))

    # --- Ejercicio 15c: ¿memoria nula o con memoria? ---
    print("\n--- Ejercicio 15c: ¿memoria nula o con memoria? ---")
    if es_memoria_nula(matriz, tolerancia):
        print("La fuente (según la matriz derivada del mensaje) es de memoria nula.")
    else:
        print("La fuente (según la matriz derivada del mensaje) es con memoria.")


if __name__ == "__main__":
    main()