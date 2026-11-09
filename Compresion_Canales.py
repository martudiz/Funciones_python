"""HACE TODO: fuente, codigos optimos (Huffman / Shannon-Fano), compresion, RLC,
codigos de bloque (Hamming, paridad) y canales.
Solo usa la libreria math.

Las funciones son las de la guia (TP4, TP5, TP6). Lo unico que se toco esta
marcado con  # CORREGIDO  o  # NUEVO  y explicado al lado."""

import math


# =============================================================================
# UTILIDADES
# =============================================================================

# NUEVO (unifica): en la guia habia 3 versiones distintas de matriz_ceros con
# parametros distintos. En un solo archivo la ultima pisa a las otras y rompe
# las anteriores, asi que queda una sola: recibe filas y columnas.
def matriz_ceros(filas, columnas):
    matriz = []
    for i in range(filas):
        fila = []
        for j in range(columnas):
            fila.append(0)
        matriz.append(fila)
    return matriz


# =============================================================================
# FUENTE
# =============================================================================

def alfabeto_y_frecuencias(lista_simbolos):
    """Pre: lista (o string) de simbolos emitidos.
    Post: dos listas paralelas: alfabeto (orden de aparicion) y frecuencias."""
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
    return [f / total for f in frecuencias]


def calcular_r(palabras_codigo):
    # r = cantidad de caracteres distintos que aparecen en las palabras codigo
    alfabeto = []
    for palabra in palabras_codigo:
        for caracter in palabra:
            if caracter not in alfabeto:
                alfabeto.append(caracter)
    return len(alfabeto)


def calcular_entropia(probabilidades, r=2):
    # H_r(S) = suma de p * log_r(1/p)
    entropia = 0.0
    for p in probabilidades:
        if p > 0:   # CORREGIDO: con p = 0 tiraba ZeroDivisionError (por convencion 0*log(1/0) = 0)
            info = math.log(1.0 / p) / math.log(r)
            entropia = entropia + p * info
    return entropia


def calcular_probs_extension(probabilidades, N):
    # P(bloque de N simbolos) = producto de las probabilidades (memoria nula)
    probs_extension = [1.0]
    for paso in range(N):
        nuevas = []
        for prob_bloque in probs_extension:
            for p in probabilidades:
                nuevas.append(prob_bloque * p)
        probs_extension = nuevas
    return probs_extension


# NUEVO: arma los nombres de los simbolos de la extension, en el MISMO orden
# que calcular_probs_extension (sin itertools, solo math y random permitidas).
def calcular_alfabeto_extension(alfabeto, N):
    alfabeto_extension = [""]
    for paso in range(N):
        nuevos = []
        for bloque in alfabeto_extension:
            for simbolo in alfabeto:
                nuevos.append(bloque + simbolo)
        alfabeto_extension = nuevos
    return alfabeto_extension


def calcular_longitud_media(probabilidades, palabras_codigo):
    # L = suma de P(simbolo) * longitud de su palabra codigo
    L = 0.0
    for i in range(len(probabilidades)):
        L = L + probabilidades[i] * len(palabras_codigo[i])
    return L


def cumple_primer_teorema_shannon(probabilidades, palabras_codigo, N):
    # H_r(S) <= Ln/N < H_r(S) + 1/N
    r = calcular_r(palabras_codigo)
    entropia = calcular_entropia(probabilidades, r)
    probs_extension = calcular_probs_extension(probabilidades, N)
    Ln = calcular_longitud_media(probs_extension, palabras_codigo)
    longitud_por_simbolo = Ln / N

    tolerancia = 1e-9
    cota_inferior_ok = longitud_por_simbolo >= entropia - tolerancia
    cota_superior_ok = longitud_por_simbolo < entropia + 1.0 / N
    return cota_inferior_ok and cota_superior_ok


def calcular_rendimiento(entropia_r, longitud_media):
    return entropia_r / longitud_media


def calcular_redundancia(rendimiento):
    return 1.0 - rendimiento


def analizar_rendimiento_y_redundancia(probabilidades, codificacion):
    r = calcular_r(codificacion)
    if r < 2:
        raise ValueError("El codigo debe tener al menos 2 simbolos distintos.")
    entropia_r = calcular_entropia(probabilidades, r)
    L = calcular_longitud_media(probabilidades, codificacion)
    rendimiento = calcular_rendimiento(entropia_r, L)
    redundancia = calcular_redundancia(rendimiento)
    return rendimiento, redundancia
# =============================================================================
# CLASIFICACION DE CODIGOS
# =============================================================================
def es_instantaneo(codigo_ordenado):
    """
    Entra: 'codigo_ordenado', la lista de palabras código YA ORDENADA
           alfabéticamente (ordenar antes alcanza para detectar cualquier
           relación de prefijo: toda palabra que extiende a otra queda
           inmediatamente después de ella al ordenar).
    Sale: True si ninguna palabra es prefijo de otra (código instantáneo /
          libre de prefijos); False si encuentra alguna.
    """
    for i in range(len(codigo_ordenado) - 1):
        actual = codigo_ordenado[i]
        siguiente = codigo_ordenado[i + 1]
        if siguiente.startswith(actual):
            return False
    return True

def no_singular(codigo):
    """
    Entra: 'codigo', lista de strings con las palabras código.
    Sale: (es_no_singular, conjunto), donde es_no_singular es True si todas las
          palabras son distintas entre sí, y 'conjunto' es el set de palabras
          código (se reutiliza en es_univocamente_decodificable()).
    """
    conjunto = set(codigo)
    return len(conjunto) == len(codigo), conjunto

def es_univocamente_decodificable(codigo, conjunto=None):
    """
    Entra: 'codigo' (lista de palabras código, en cualquier orden) y 'conjunto'
           (opcional; el set de esas palabras, ver no_singular()).
    Sale: True si el código es unívocamente decodificable, False si no.
    
    Implementa el algoritmo de Sardinas-Patterson con impresiones paso a paso.
    """
    if conjunto is None:
        conjunto = set(codigo)

    print(f"\n--- Iniciando Teorema de Sardinas-Patterson ---")
    print(f"Generación 0 (Código original): {conjunto}")

    # Generación 1: sufijos pendientes de comparar cada palabra código con las
    # demás palabras código de las que es prefijo.
    generacion = set()
    for a in codigo:
        for b in codigo:
            if a != b and b.startswith(a):
                generacion.add(b[len(a):])
                
    print(f"Generación 1 (Sufijos iniciales): {generacion}")
    
    # Verificamos conflicto en la Generación 1
    interseccion = generacion.intersection(conjunto)
    if interseccion:
        print(f"En la Generación 1 el sufijo {interseccion} pertenece al código original.")
        return False

    generaciones_vistas = [frozenset(generacion)]
    numero_generacion = 2

    # Ciclo de generaciones siguientes
    while generacion:
        siguiente = set()
        for sufijo in generacion:
            for c in codigo:
                if sufijo.startswith(c):
                    siguiente.add(sufijo[len(c):])
                elif c.startswith(sufijo):
                    siguiente.add(c[len(sufijo):])
                    
        print(f"Generación {numero_generacion}: {siguiente}")
        
        # 1. Condición de parada por conflicto (no unívocamente decodificable)
        interseccion = siguiente.intersection(conjunto)
        if interseccion:
            print(f"En la Generación {numero_generacion} el sufijo {interseccion} pertenece al código original.")
            return False
            
        # 2. Condición de parada por ciclo infinito sin conflictos (unívocamente decodificable)
        if frozenset(siguiente) in generaciones_vistas:
            print(f"La Generación {numero_generacion} ya apareció antes. ¡Se detectó un ciclo cerrado sin conflictos!")
            return True
            
        generaciones_vistas.append(frozenset(siguiente))
        generacion = siguiente
        numero_generacion += 1

    # 3. Condición de parada por conjunto vacío (unívocamente decodificable)
    print("✅ No se generaron más sufijos nuevos (Generación vacía).")
    return True

# =============================================================================
# CODIGOS OPTIMOS Y COMPRESION
# =============================================================================

def crea_huffman(probs):
    cant = len(probs)
    codigos = []
    probs_grupos = []
    simbolos_grupos = []
    for i in range(cant):
        codigos.append("")
        probs_grupos.append(probs[i])
        simbolos_grupos.append([i])

    while len(probs_grupos) > 1:
        # indice del menor
        a = 0
        for k in range(len(probs_grupos)):
            if probs_grupos[k] < probs_grupos[a]:
                a = k
        # indice del segundo menor (salteando a)
        b = -1
        for k in range(len(probs_grupos)):
            if k != a:
                if b == -1 or probs_grupos[k] < probs_grupos[b]:
                    b = k

        for s in simbolos_grupos[a]:
            codigos[s] = "0" + codigos[s]
        for s in simbolos_grupos[b]:
            codigos[s] = "1" + codigos[s]

        nueva_prob = probs_grupos[a] + probs_grupos[b]
        nuevos_simbolos = simbolos_grupos[a] + simbolos_grupos[b]

        # borrar primero el indice mayor
        if a > b:
            probs_grupos.pop(a)
            simbolos_grupos.pop(a)
            probs_grupos.pop(b)
            simbolos_grupos.pop(b)
        else:
            probs_grupos.pop(b)
            simbolos_grupos.pop(b)
            probs_grupos.pop(a)
            simbolos_grupos.pop(a)

        probs_grupos.append(nueva_prob)
        simbolos_grupos.append(nuevos_simbolos)

    return codigos


def ordena_indices(probs):
    indices = []
    for i in range(len(probs)):
        indices.append(i)
    # burbuja de mayor a menor, estable
    for i in range(len(indices) - 1):
        for k in range(len(indices) - 1 - i):
            if probs[indices[k]] < probs[indices[k + 1]]:
                auxiliar = indices[k]
                indices[k] = indices[k + 1]
                indices[k + 1] = auxiliar
    return indices


def divide(indices, probs, codigos):
    if len(indices) == 1:
        return

    total = 0
    for i in indices:
        total = total + probs[i]

    corte = 1
    mejor_dif = total + 1
    izq = 0
    for k in range(1, len(indices)):
        izq = izq + probs[indices[k - 1]]
        diferencia = abs(total - 2 * izq)
        if diferencia < mejor_dif:
            mejor_dif = diferencia
            corte = k

    izquierda = indices[:corte]
    derecha = indices[corte:]

    for i in izquierda:
        codigos[i] = codigos[i] + "0"
    for i in derecha:
        codigos[i] = codigos[i] + "1"

    divide(izquierda, probs, codigos)
    divide(derecha, probs, codigos)


def crea_shannon_fano(probs):
    codigos = []
    for i in range(len(probs)):
        codigos.append("")
    indices = ordena_indices(probs)
    divide(indices, probs, codigos)
    return codigos


def codificar(alfabeto, codigos, mensaje):
    # primer byte = cantidad de bits de relleno agregados al final
    bits = ""
    for caracter in mensaje:
        posicion = alfabeto.index(caracter)
        bits = bits + codigos[posicion]

    relleno = (8 - len(bits) % 8) % 8
    for i in range(relleno):
        bits = bits + "0"

    resultado = bytearray()
    resultado.append(relleno)
    for k in range(0, len(bits), 8):
        grupo = bits[k:k + 8]
        resultado.append(int(grupo, 2))
    return resultado


def decodificar(alfabeto, codigos, secuencia):
    relleno = secuencia[0]
    bits = ""
    for k in range(1, len(secuencia)):
        bits = bits + format(secuencia[k], "08b")
    bits = bits[:len(bits) - relleno]

    mensaje = ""
    acumulado = ""
    for bit in bits:
        acumulado = acumulado + bit
        if acumulado in codigos:
            posicion = codigos.index(acumulado)
            mensaje = mensaje + alfabeto[posicion]
            acumulado = ""
    return mensaje


# NUEVO: es la funcion que el main de tu amiga llamaba y no existia.
# Igual que decodificar() pero desde una cadena de bits ("0101...") en vez
# de un bytearray. Sirve si en el parcial te dan el mensaje ya codificado.
# Asume codigo instantaneo (como Huffman / Shannon-Fano).
def decodificar_bits(bits, alfabeto, codigos):
    simbolos = []
    acumulado = ""
    for bit in bits:
        acumulado = acumulado + bit
        if acumulado in codigos:
            simbolos.append(alfabeto[codigos.index(acumulado)])
            acumulado = ""
    if acumulado != "":
        print(f"  AVISO: sobraron bits sin decodificar: '{acumulado}'")
    return simbolos


def tasa_compresion(mensaje, secuencia):
    # cada simbolo original ocupa 1 byte (ASCII)
    tamano_original = len(mensaje)
    tamano_comprimido = len(secuencia)
    return tamano_original / tamano_comprimido


def comprimir_rlc(mensaje):
    # 1 byte para el caracter (ASCII) + 1 byte para la cantidad de repeticiones
    resultado = bytearray()
    if len(mensaje) == 0:          # CORREGIDO: con mensaje vacio tiraba IndexError
        return resultado
    caracter_actual = mensaje[0]
    contador = 1

    for k in range(1, len(mensaje)):
        # CORREGIDO: un byte guarda hasta 255; si la racha es mas larga se corta en varias
        if mensaje[k] == caracter_actual and contador < 255:
            contador = contador + 1
        else:
            resultado.append(ord(caracter_actual))
            resultado.append(contador)
            caracter_actual = mensaje[k]
            contador = 1

    resultado.append(ord(caracter_actual))
    resultado.append(contador)
    return resultado


# =============================================================================
# CODIGOS DE BLOQUE: HAMMING Y PARIDAD
# =============================================================================

def distancia_hamming(palabra1, palabra2):
    distancia = 0
    for i in range(len(palabra1)):
        if palabra1[i] != palabra2[i]:
            distancia += 1
    return distancia


def analizar_capacidades_codigo(codificacion):
    # d_min = menor distancia entre todos los pares de palabras
    d_min = float('inf')
    cantidad_palabras = len(codificacion)
    for i in range(cantidad_palabras):
        for j in range(i + 1, cantidad_palabras):
            dist = distancia_hamming(codificacion[i], codificacion[j])
            if dist < d_min:
                d_min = dist
    if d_min == float('inf'):
        d_min = 0

    errores_detectables = max(0, d_min - 1)          # d - 1
    errores_corregibles = max(0, (d_min - 1) // 2)   # piso((d - 1) / 2)
    return d_min, errores_detectables, errores_corregibles


def codificar_con_paridad(caracter):
    # ASCII 7 bits corrido 1 a la izquierda, paridad PAR en el bit menos significativo
    ascii_val = ord(caracter) & 0x7F
    cantidad_unos = bin(ascii_val).count('1')
    bit_paridad = cantidad_unos % 2
    byte_codificado = (ascii_val << 1) | bit_paridad
    return byte_codificado


def verificar_byte(byte):
    cantidad_unos = bin(byte).count('1')
    if cantidad_unos % 2 == 0:
        return True, chr(byte >> 1)
    else:
        return False, None


def generar_matriz_bits(cadena):
    # cada fila: 7 bits ASCII + 1 bit de paridad par (VRC)
    matriz = []
    for i in range(len(cadena)):
        ascii_val = ord(cadena[i])
        fila_bits = []
        for j in range(6, -1, -1):
            fila_bits.append((ascii_val >> j) & 1)

        suma_bits = 0
        for j in range(len(fila_bits)):
            suma_bits += fila_bits[j]
        fila_bits.append(suma_bits % 2)
        matriz.append(fila_bits)
    return matriz


def calcular_lrc(matriz):
    # paridad de cada columna (LRC); el bit 8 es la paridad cruzada
    lrc = [0] * 8
    for col in range(8):
        suma_columna = 0
        for fila in range(len(matriz)):
            suma_columna += matriz[fila][col]
        lrc[col] = suma_columna % 2
    return lrc


def codificar_mensaje(cadena):
    matriz = generar_matriz_bits(cadena)
    fila_lrc = calcular_lrc(matriz)
    matriz.insert(0, fila_lrc)         # LRC como PRIMER byte

    secuencia_bytes = bytearray()
    for i in range(len(matriz)):
        byte_val = 0
        for j in range(8):
            byte_val = (byte_val << 1) | matriz[i][j]
        secuencia_bytes.append(byte_val)
    return secuencia_bytes


def bytes_a_matriz(secuencia_bytes):
    matriz = []
    for i in range(len(secuencia_bytes)):
        byte_val = secuencia_bytes[i]
        fila_bits = []
        for j in range(7, -1, -1):
            fila_bits.append((byte_val >> j) & 1)
        matriz.append(fila_bits)
    return matriz


def decodificar_y_corregir(secuencia_bytes):
    matriz = bytes_a_matriz(secuencia_bytes)
    # CORREGIDO: se saco el print(matriz) de depuracion

    # filas de datos con paridad mal (la ultima fila es el LRC)
    error_fila = -1
    conteo_errores_filas = 0
    for i in range(1, len(matriz)):
        suma_fila = 0
        for j in range(8):
            suma_fila += matriz[i][j]
        if suma_fila % 2 != 0:
            error_fila = i
            conteo_errores_filas += 1

    # columnas con paridad mal (incluye la fila LRC)
    error_columna = -1
    conteo_errores_cols = 0
    for j in range(8):
        suma_columna = 0
        for i in range(len(matriz)):
            suma_columna += matriz[i][j]
        if suma_columna % 2 != 0:
            error_columna = j
            conteo_errores_cols += 1

    if conteo_errores_filas == 1 and conteo_errores_cols == 1:
        # 1 error: esta en el cruce de la fila y la columna marcadas
        matriz[error_fila][error_columna] ^= 1
    elif conteo_errores_filas == 0 and conteo_errores_cols <= 1:
        # sin errores, o 1 error en la propia fila LRC (los datos estan bien)
        pass
    else:
        # CORREGIDO: antes el caso "1 fila mal y 0 columnas mal" (que solo pasa
        # con 2+ errores) devolvia el mensaje SIN corregir como si estuviera bien
        return ""

    mensaje = ""
    for i in range(1, len(matriz)):
        ascii_val = 0
        for j in range(7):          # los 7 bits ASCII, se descarta el de paridad
            ascii_val = (ascii_val << 1) | matriz[i][j]
        mensaje += chr(ascii_val)
    return mensaje


# =============================================================================
# CANALES   (matriz del canal: fila = entrada a_i, columna = salida b_j,
#            elemento = P(b_j / a_i))
# =============================================================================

def calcular_alfabeto(cadena):
    alfabeto = []
    for caracter in cadena:
        if caracter not in alfabeto:
            alfabeto.append(caracter)
    return len(alfabeto), alfabeto


def normalizar_filas(matriz):
    for i in range(len(matriz)):
        total_fila = 0
        for j in range(len(matriz[i])):
            total_fila = total_fila + matriz[i][j]
        for j in range(len(matriz[i])):
            matriz[i][j] = matriz[i][j] / total_fila
    return matriz


def genero_matriz(cadena_entrada, cadena_salida):
    cant_filas, alfabeto_entrada = calcular_alfabeto(cadena_entrada)
    cant_columnas, alfabeto_salida = calcular_alfabeto(cadena_salida)
    matriz = matriz_ceros(cant_filas, cant_columnas)
    for k in range(len(cadena_entrada)):
        fila = alfabeto_entrada.index(cadena_entrada[k])
        columna = alfabeto_salida.index(cadena_salida[k])
        matriz[fila][columna] = matriz[fila][columna] + 1
    matriz = normalizar_filas(matriz)
    return matriz


def crea_matriz_eventos_simultaneos(matriz_canal, prob_priori):
    # P(a_i, b_j) = P(a_i) * P(b_j / a_i)
    cant_filas = len(matriz_canal)
    cant_columnas = len(matriz_canal[0])
    matriz_eventos_simultaneos = matriz_ceros(cant_filas, cant_columnas)
    for i in range(cant_filas):
        for j in range(cant_columnas):
            matriz_eventos_simultaneos[i][j] = matriz_canal[i][j] * prob_priori[i]
    return matriz_eventos_simultaneos


def crea_prob_simbolos_salida(matriz_eventos_simultaneos):
    # P(b_j) = suma de la columna j de simultaneos
    cant_columnas = len(matriz_eventos_simultaneos[0])
    cant_filas = len(matriz_eventos_simultaneos)
    prob_simbolos_salida = []
    for j in range(cant_columnas):
        prob_simbolos_salida.append(0)
    for j in range(cant_columnas):
        for i in range(cant_filas):
            prob_simbolos_salida[j] += matriz_eventos_simultaneos[i][j]
    return prob_simbolos_salida


def crea_matriz_posteriori(matriz_eventos_simultaneos, prob_simbolos_salida):
    # P(a_i / b_j) = P(a_i, b_j) / P(b_j)
    cant_filas = len(matriz_eventos_simultaneos)
    cant_columnas = len(matriz_eventos_simultaneos[0])
    matriz_posteriori = matriz_ceros(cant_filas, cant_columnas)
    for i in range(cant_filas):
        for j in range(cant_columnas):
            if prob_simbolos_salida[j] > 0:   # CORREGIDO: si una salida nunca ocurre, dividia por 0
                matriz_posteriori[i][j] = matriz_eventos_simultaneos[i][j] / prob_simbolos_salida[j]
    return matriz_posteriori


def calcula_valores(prob_priori, matriz_canal):
    simultaneos = crea_matriz_eventos_simultaneos(matriz_canal, prob_priori)
    prob_salida = crea_prob_simbolos_salida(simultaneos)
    posteriori = crea_matriz_posteriori(simultaneos, prob_salida)
    return posteriori, prob_salida, simultaneos


# OJO: esta recibe (matriz_canal, prob_priori), al reves que las demas (asi esta en la guia)
def calcular_entropia_posteriori(matriz_canal, prob_priori):
    # H(A / b_j) = suma_i P(a_i / b_j) * log2(1 / P(a_i / b_j)), una por columna
    simultaneos = crea_matriz_eventos_simultaneos(matriz_canal, prob_priori)
    prob_salida = crea_prob_simbolos_salida(simultaneos)
    posteriori = crea_matriz_posteriori(simultaneos, prob_salida)
    cant_columnas = len(posteriori[0])
    entropia_posteriori = []
    for j in range(cant_columnas):
        entropia_posteriori.append(0)
    for j in range(cant_columnas):
        for i in range(len(posteriori)):
            if posteriori[i][j] > 0:
                entropia_posteriori[j] = entropia_posteriori[j] + posteriori[i][j] * math.log2(1 / posteriori[i][j])
    return entropia_posteriori


def calcular_equivocacion(prob_priori, matriz_canal):
    # H(A/B) = suma P(a,b) * log2(1 / P(a/b))   (equivocacion o ruido)
    posteriori, prob_salida, simultaneos = calcula_valores(prob_priori, matriz_canal)
    acumulado = 0
    for i in range(len(posteriori)):
        for j in range(len(posteriori[0])):
            if simultaneos[i][j] > 0:
                acumulado = acumulado + simultaneos[i][j] * math.log2(1 / posteriori[i][j])
    return acumulado


def calcular_perdida(prob_priori, matriz_canal):
    # H(B/A) = suma P(a,b) * log2(1 / P(b/a))
    posteriori, prob_salida, simultaneos = calcula_valores(prob_priori, matriz_canal)
    acumulado = 0
    for i in range(len(posteriori)):
        for j in range(len(posteriori[0])):
            if simultaneos[i][j] > 0:
                acumulado = acumulado + simultaneos[i][j] * math.log2(1 / matriz_canal[i][j])
    return acumulado


def calcular_entropia_afin(prob_priori, matriz_canal):
    # H(A,B) = suma P(a,b) * log2(1 / P(a,b))
    posteriori, prob_salida, simultaneos = calcula_valores(prob_priori, matriz_canal)
    acumulado = 0
    for i in range(len(posteriori)):
        for j in range(len(posteriori[0])):
            if simultaneos[i][j] > 0:
                acumulado = acumulado + simultaneos[i][j] * math.log2(1 / simultaneos[i][j])
    return acumulado


def calcular_informacion_mutua(prob_priori, matriz_canal):
    # I(A;B) = suma P(a,b) * log2( P(a,b) / (P(a) P(b)) )
    posteriori, prob_salida, simultaneos = calcula_valores(prob_priori, matriz_canal)
    acumulado = 0
    for i in range(len(posteriori)):
        for j in range(len(posteriori[0])):
            if simultaneos[i][j] > 0:
                acumulado = acumulado + simultaneos[i][j] * math.log2(simultaneos[i][j] / (prob_priori[i] * prob_salida[j]))
    return acumulado


def canal_sin_ruido(matriz_canal):
    # sin ruido: en cada columna hay a lo sumo UN elemento distinto de 0
    cant_filas = len(matriz_canal)
    cant_columnas = len(matriz_canal[0])
    for j in range(cant_columnas):
        cont = 0
        for i in range(cant_filas):
            if matriz_canal[i][j] != 0:
                cont += 1
        if cont > 1:   # CORREGIDO: era  cont != 1, y una columna toda en 0 (salida que nunca sale) lo hacia dar False
            return False
    return True


def canal_determinante(matriz_canal):
    # determinante: en cada fila hay exactamente UN elemento distinto de 0
    cant_filas = len(matriz_canal)
    cant_columnas = len(matriz_canal[0])
    for i in range(cant_filas):
        cont = 0
        for j in range(cant_columnas):
            if matriz_canal[i][j] != 0:
                cont += 1
        if cont != 1:
            return False
    return True


def verifica_canal_uniforme(matriz_canal):
    # uniforme: todas las filas son permutaciones de la primera
    primera_fila = sorted(matriz_canal[0])
    for i in range(len(matriz_canal)):
        if sorted(matriz_canal[i]) != primera_fila:
            return False
    return True


def calcular_capacidad(matriz_canal):
    # devuelve -1 si el canal no es de ninguno de los 3 tipos
    if canal_determinante(matriz_canal):
        # CORREGIDO: era log2(cant. columnas); si hay una columna toda en 0 esa salida
        # no se puede alcanzar, asi que se cuentan solo las columnas usadas
        cols_usadas = 0
        for j in range(len(matriz_canal[0])):
            usada = False
            for i in range(len(matriz_canal)):
                if matriz_canal[i][j] != 0:
                    usada = True
            if usada:
                cols_usadas += 1
        return math.log2(cols_usadas)
    elif canal_sin_ruido(matriz_canal):
        return math.log2(len(matriz_canal))
    elif verifica_canal_uniforme(matriz_canal):
        cant_columnas = len(matriz_canal[0])
        entropia = 0
        for j in range(cant_columnas):
            if matriz_canal[0][j] > 0:
                entropia = entropia + matriz_canal[0][j] * math.log2(1 / matriz_canal[0][j])
        return math.log2(cant_columnas) - entropia
    else:
        return -1


def estimativo_capacidad_canal(matriz_canal_binario, valor_paso):
    # barre P(a1) = paso, 2*paso, ... y se queda con la mayor I(A;B)
    cant_pasos = round(1 / valor_paso)
    capacidad = -1
    prob_capacidad = 0
    for k in range(1, cant_pasos):
        p = k * valor_paso
        prob_priori = [p, 1 - p]
        inf_mutua = calcular_informacion_mutua(prob_priori, matriz_canal_binario)
        if inf_mutua > capacidad:
            capacidad = inf_mutua
            prob_capacidad = p
    return capacidad, prob_capacidad


def calcula_prob_error(prob_priori, matriz_canal):
    # regla de MAXIMA POSIBILIDAD: para cada salida b_j se decide la entrada con
    # mayor P(b_j / a_i) (maximo de la columna de la matriz del canal)
    cant_columnas = len(matriz_canal[0])
    cant_filas = len(matriz_canal)
    matriz_eventos_simultaneos = crea_matriz_eventos_simultaneos(matriz_canal, prob_priori)
    prob_aciertos = 0
    for j in range(cant_columnas):
        maximo = -1
        fila_max = 0
        for i in range(cant_filas):
            if matriz_canal[i][j] > maximo:
                maximo = matriz_canal[i][j]
                fila_max = i
        prob_aciertos = prob_aciertos + matriz_eventos_simultaneos[fila_max][j]
    return 1 - prob_aciertos


def crea_matriz_compuesta(matriz_canal1, matriz_canal2):
    # canales en serie: producto de matrices
    cant_filas = len(matriz_canal1)
    cant_columnas = len(matriz_canal2[0])
    dimension = len(matriz_canal2)
    matriz_compuesta = matriz_ceros(cant_filas, cant_columnas)
    for i in range(cant_filas):
        for k in range(cant_columnas):
            for j in range(dimension):
                matriz_compuesta[i][k] = matriz_compuesta[i][k] + matriz_canal1[i][j] * matriz_canal2[j][k]
    return matriz_compuesta


def se_puede_combinar(matriz_canal, col1, col2):
    # reduccion suficiente: las dos columnas tienen que ser proporcionales
    cant_filas = len(matriz_canal)
    r = -1
    for i in range(cant_filas):
        if matriz_canal[i][col2] != 0 and r == -1:
            r = i
    if r == -1:
        return True
    for i in range(cant_filas):
        izquierda = matriz_canal[i][col1] * matriz_canal[r][col2]
        derecha = matriz_canal[r][col1] * matriz_canal[i][col2]
        if abs(izquierda - derecha) > 0.000001:
            return False
    return True


def genera_matriz_determinante(matriz_canal, col1, col2):
    if not se_puede_combinar(matriz_canal, col1, col2):
        return False
    cant_columnas_original = len(matriz_canal[0])
    cant_columnas = cant_columnas_original - 1
    matriz_determinante = matriz_ceros(cant_columnas_original, cant_columnas)
    menor = min(col1, col2)
    mayor = max(col1, col2)
    for j in range(cant_columnas_original):
        if j < mayor:
            destino = j
        elif j == mayor:
            destino = menor
        else:
            destino = j - 1
        matriz_determinante[j][destino] = 1
    return matriz_determinante


def busca_par_combinable(matriz_canal):
    cant_columnas = len(matriz_canal[0])
    for col1 in range(cant_columnas):
        for col2 in range(col1 + 1, cant_columnas):
            if se_puede_combinar(matriz_canal, col1, col2):
                return col1, col2
    return -1, -1


def reduce_canal(matriz_canal):
    while True:
        col1, col2 = busca_par_combinable(matriz_canal)
        if col1 == -1:
            break
        determinante = genera_matriz_determinante(matriz_canal, col1, col2)
        matriz_canal = crea_matriz_compuesta(matriz_canal, determinante)
    return matriz_canal


# =============================================================================
# INFORMES (solo imprimen; todo el calculo esta en las funciones de arriba)
# =============================================================================

DECIMALES = 4    # nunca menos de 3


def fmt(x):
    return f"{x:.{DECIMALES}f}"


def imprimo_titulo(texto):
    print("\n" + "=" * 60)
    print(" " + texto)
    print("=" * 60)


def imprimo_vector(etiquetas, valores):
    ancho = max(9, max(len(str(e)) for e in etiquetas) + 2)
    print("".join(f"{str(e):>{ancho}}" for e in etiquetas))
    print("".join(f"{fmt(v):>{ancho}}" for v in valores))


def imprimo_matriz(matriz, etiq_filas, etiq_cols):
    ancho = max(9, max(len(str(e)) for e in etiq_filas + etiq_cols) + 2)
    print(" " * ancho + "".join(f"{str(c):>{ancho}}" for c in etiq_cols))
    for etiqueta, fila in zip(etiq_filas, matriz):
        print(f"{str(etiqueta):>{ancho}}" + "".join(f"{fmt(x):>{ancho}}" for x in fila))


def imprimo_tabla_codigo(alfabeto, probs, codigos):
    ancho = max(8, max(len(str(s)) for s in alfabeto) + 2)
    print(f"  {'Simbolo':>{ancho}} {'P(s)':>10} {'Codigo':>14} {'l':>4}")
    for s, p, c in zip(alfabeto, probs, codigos):
        print(f"  {str(s):>{ancho}} {fmt(p):>10} {c:>14} {len(c):>4}")


def informe_fuente(alfabeto, probs, N):
    imprimo_titulo("FUENTE")
    print(f"  {'Simbolo':>8} {'P(s)':>10} {'I(s) bits':>10}")
    for s, p in zip(alfabeto, probs):
        print(f"  {str(s):>8} {fmt(p):>10} {fmt(math.log2(1 / p)):>10}")
    H = calcular_entropia(probs, 2)
    print(f"\n  H(S)    = {fmt(H)} bits/simbolo")
    print(f"  log2(q) = {fmt(math.log2(len(alfabeto)))} bits (maximo, si fuera equiprobable)")
    if N > 1:
        H_ext = calcular_entropia(calcular_probs_extension(probs, N), 2)
        print(f"  H(S^{N})  = {fmt(H_ext)}   {N}*H(S) = {fmt(N * H)}   (tienen que coincidir)")


def informe_codigo_optimo(alfabeto, probs, N, metodo):
    alfabeto_ext = calcular_alfabeto_extension(alfabeto, N)
    probs_ext = calcular_probs_extension(probs, N)
    if metodo == "huffman":
        codigos = crea_huffman(probs_ext)
    else:
        codigos = crea_shannon_fano(probs_ext)

    nombre = "LA FUENTE" if N == 1 else f"LA EXTENSION DE ORDEN {N}"
    imprimo_titulo(f"CODIGO {metodo.upper()} PARA {nombre}")
    imprimo_tabla_codigo(alfabeto_ext, probs_ext, codigos)

    r = calcular_r(codigos)
    H = calcular_entropia(probs, r)
    H_ext = calcular_entropia(probs_ext, r)
    L = calcular_longitud_media(probs_ext, codigos)
    rendimiento = calcular_rendimiento(H_ext, L)
    print(f"\n  r = {r}")
    if N == 1:
        print(f"  L = {fmt(L)} simbolos codigo / simbolo fuente")
    else:
        print(f"  L_{N} = {fmt(L)}   L_{N}/{N} = {fmt(L / N)} por simbolo ORIGINAL")
    print(f"  Rendimiento  = H/L = {fmt(rendimiento)}   Redundancia = {fmt(calcular_redundancia(rendimiento))}")

    cumple = cumple_primer_teorema_shannon(probs, codigos, N)
    print(f"\n  Primer Teorema de Shannon: H_r(S) <= L_N/N < H_r(S) + 1/N")
    print(f"  {fmt(H)} <= {fmt(L / N)} < {fmt(H + 1 / N)}  -> {'CUMPLE' if cumple else 'NO CUMPLE'}")
    return codigos


def informe_compresion(alfabeto, codigos, mensaje):
    imprimo_titulo("COMPRESION")
    secuencia = codificar(alfabeto, codigos, mensaje)
    bits = ""
    for s in mensaje:
        bits = bits + codigos[alfabeto.index(s)]
    print(f"  Mensaje       : {''.join(mensaje)}  ({len(mensaje)} simbolos = {len(mensaje)} bytes)")
    print(f"  Bits          : {bits}  ({len(bits)} bits)")
    print(f"  Bytes         : {list(secuencia)}  ({len(secuencia)} bytes, el 1ro = bits de relleno = {secuencia[0]})")
    print(f"  En binario    : {' '.join(format(b, '08b') for b in secuencia)}")
    print(f"  Tasa de compresion = {len(mensaje)} / {len(secuencia)} = {fmt(tasa_compresion(mensaje, secuencia))} : 1")
    recuperado = decodificar(alfabeto, codigos, secuencia)
    print(f"  Decodificado  : {recuperado}  (coincide: {recuperado == ''.join(mensaje)})")


def informe_codigo_dado(probs, codigo, N):
    imprimo_titulo("CODIGO DADO: RENDIMIENTO, REDUNDANCIA Y PRIMER TEOREMA")
    es_no_singular, conjunto = no_singular(codigo)
    print(f"  No singular: {es_no_singular}")
    if es_no_singular:
        instantaneo = es_instantaneo(sorted(codigo))
        print(f"  Instantaneo: {instantaneo}")
        if instantaneo:
            print("  Univocamente decodificable: True (todo instantaneo lo es)")
        else:
            print(f"  Univocamente decodificable: {es_univocamente_decodificable(codigo, conjunto)}")
    else:
        print("  Instantaneo: False   Univocamente decodificable: False")
    r = calcular_r(codigo)
    if N == 1:
        rendimiento, redundancia = analizar_rendimiento_y_redundancia(probs, codigo)
        print(f"  Codigo: {codigo}   r = {r}")
        print(f"  H_r(S) = {fmt(calcular_entropia(probs, r))}   L = {fmt(calcular_longitud_media(probs, codigo))}")
        print(f"  Rendimiento = {fmt(rendimiento)} ({rendimiento * 100:.2f}%)   Redundancia = {fmt(redundancia)} ({redundancia * 100:.2f}%)")
    else:
        print(f"  Codigo para la extension de orden {N}: {codigo}   r = {r}")
    H = calcular_entropia(probs, r)
    Ln = calcular_longitud_media(calcular_probs_extension(probs, N), codigo)
    cumple = cumple_primer_teorema_shannon(probs, codigo, N)
    print(f"  Primer Teorema (N = {N}): {fmt(H)} <= {fmt(Ln / N)} < {fmt(H + 1 / N)}  -> {'CUMPLE' if cumple else 'NO CUMPLE'}")


def informe_hamming(codigo):
    imprimo_titulo("DISTANCIA DE HAMMING")
    d, detectables, corregibles = analizar_capacidades_codigo(codigo)
    print(f"  Codigo: {codigo}")
    print(f"  Distancia minima d = {d}")
    print(f"  Detecta  d - 1        = {detectables} errores")
    print(f"  Corrige  (d - 1) // 2 = {corregibles} errores")


def informe_paridad(texto, errores, bytes_recibidos):
    imprimo_titulo("PARIDAD")
    print("  Paridad simple (LSB, paridad par):")
    for c in texto:
        byte = codificar_con_paridad(c)
        ok, _ = verificar_byte(byte)
        print(f"    '{c}'  ASCII {ord(c):>3} = {format(ord(c), '07b')}  ->  {format(byte, '08b')} ({byte})  verifica: {ok}")

    secuencia = codificar_mensaje(texto)
    print("\n  Paridad cruzada (cada fila = 7 bits + VRC, primer byte = LRC):")
    for b in secuencia:
        print(f"    {format(b, '08b')}  ({b})")
    print(f"  Bytes: {list(secuencia)}")
    print(f"  Decodificado sin errores: '{decodificar_y_corregir(secuencia)}'")

    if len(errores) > 0:
        con_error = bytearray(secuencia)
        for byte_i, bit_i in errores:      # bit_i: 0 = el de la izquierda, 7 = el de paridad
            con_error[byte_i] ^= (1 << (7 - bit_i))
        print(f"\n  Con errores en (byte, bit) {errores}: {list(con_error)}")
        print(f"  Decodificado: '{decodificar_y_corregir(con_error)}'  ('' = no se puede corregir)")

    if len(bytes_recibidos) > 0:
        print(f"\n  Secuencia recibida: {bytes_recibidos}")
        print(f"  Decodificado: '{decodificar_y_corregir(bytearray(bytes_recibidos))}'  ('' = no se puede corregir)")


def informe_rlc(mensaje):
    imprimo_titulo("RLC")
    secuencia = comprimir_rlc(mensaje)
    print(f"  Mensaje: {mensaje}  ({len(mensaje)} bytes)")
    print(f"  Bytes  : {list(secuencia)}  ({len(secuencia)} bytes)")
    print(f"  Tasa de compresion = {fmt(tasa_compresion(mensaje, secuencia))} : 1")


def informe_canal(prob_priori, matriz, alf_entrada, alf_salida, paso):
    imprimo_titulo("CANAL")
    print("  Matriz del canal P(b/a)  (fila = entrada, columna = salida):")
    imprimo_matriz(matriz, alf_entrada, alf_salida)
    print("\n  Probabilidades a priori P(a):")
    imprimo_vector(alf_entrada, prob_priori)

    posteriori, prob_salida, simultaneos = calcula_valores(prob_priori, matriz)
    print("\n  Probabilidades de salida P(b):")
    imprimo_vector(alf_salida, prob_salida)
    print("\n  Eventos simultaneos P(a,b):")
    imprimo_matriz(simultaneos, alf_entrada, alf_salida)
    print("\n  A posteriori P(a/b):")
    imprimo_matriz(posteriori, alf_entrada, alf_salida)
    print("\n  Entropias a posteriori H(A/b):")
    imprimo_vector(alf_salida, calcular_entropia_posteriori(matriz, prob_priori))

    HA = calcular_entropia(prob_priori, 2)
    HB = calcular_entropia(prob_salida, 2)
    equivocacion = calcular_equivocacion(prob_priori, matriz)
    perdida = calcular_perdida(prob_priori, matriz)
    afin = calcular_entropia_afin(prob_priori, matriz)
    mutua = calcular_informacion_mutua(prob_priori, matriz)
    print(f"\n  H(A)   (a priori)            = {fmt(HA)}")
    print(f"  H(B)                         = {fmt(HB)}")
    print(f"  H(A/B) equivocacion o ruido  = {fmt(equivocacion)}")
    print(f"  H(B/A) perdida               = {fmt(perdida)}")
    print(f"  H(A,B) entropia afin         = {fmt(afin)}")
    print(f"  I(A;B) informacion mutua     = {fmt(mutua)}")
    print(f"  Control: H(A)-H(A/B) = {fmt(HA - equivocacion)}   H(B)-H(B/A) = {fmt(HB - perdida)}   H(A)+H(B)-H(A,B) = {fmt(HA + HB - afin)}")

    print(f"\n  Sin ruido: {canal_sin_ruido(matriz)}   Determinante: {canal_determinante(matriz)}   Uniforme: {verifica_canal_uniforme(matriz)}")
    capacidad = calcular_capacidad(matriz)
    if capacidad != -1:
        print(f"  Capacidad (formula) = {fmt(capacidad)}")
    else:
        print("  Capacidad (formula): no es determinante, sin ruido ni uniforme")
    if len(matriz) == 2:
        cap, p = estimativo_capacidad_canal(matriz, paso)
        print(f"  Capacidad estimada (paso {paso}) = {fmt(cap)}  con P(a1) = {fmt(p)}, P(a2) = {fmt(1 - p)}")

    print(f"\n  P(error), maxima posibilidad = {fmt(calcula_prob_error(prob_priori, matriz))}")

    reducida = reduce_canal(matriz)
    if len(reducida[0]) < len(matriz[0]):
        print("\n  Canal reducido (todas las reducciones suficientes posibles):")
        etiquetas = [f"b{j + 1}'" for j in range(len(reducida[0]))]
        imprimo_matriz(reducida, alf_entrada, etiquetas)
    else:
        print("\n  No admite reducciones suficientes")


def informe_canales_serie(matriz1, matriz2):
    imprimo_titulo("CANALES EN SERIE")
    compuesta = crea_matriz_compuesta(matriz1, matriz2)
    etiq_filas = [f"a{i + 1}" for i in range(len(compuesta))]
    etiq_cols = [f"c{k + 1}" for k in range(len(compuesta[0]))]
    imprimo_matriz(compuesta, etiq_filas, etiq_cols)

# =============================================================================
# MAIN
# =============================================================================

"""CONFIGURACION: lo unico que hay que tocar en el parcial.
Cada bloque se saltea si su dato principal esta vacio ("" o [])."""

# --- 1) FUENTE -> CODIGO OPTIMO -> COMPRESION ---
mensaje = ""   # msj de la consigna
# si en vez del msj te dan los BITS y el codigo, completar esto (tiene prioridad):
mensajeBits = ""
codigoDado = []
alfabetoDado = []  # ([] => S1, S2, ...)
ordenExt = 1                         # orden de la extension
metodo = "huffman"                   # "huffman" o "shannon-fano"

# --- 2) CODIGO DADO CON SUS PROBABILIDADES ---
codigoAnalizar = ["0", "10", "11"]                  # ej: ['(]', ']', '[)', ')', '([']
probsCodigo = [0.5 , 0.25 , 0.25]                     # ej: [0.15, 0.25, 0.05, 0.45, 0.10]
ordenCodigo = 1                      # N si el codigo es para la extension de orden N

# --- 3) HAMMING ---
codigoBloque = []                    # ej: ["00000", "01101", "10110", "11011"]

# --- 4) PARIDAD ---
textoParidad = ""                    # ej: "HOLA"
errores = []                         # ej: [(1, 5)] = invierte el bit 5 del byte 1 (el byte 0 es el LRC)
bytesRecibidos = []                  # si te dan la secuencia recibida, ej: [20, 144, 159, 153, 130]

# --- 5) RLC ---
mensajeRLC = ""                      # ej: "AAAABBBCCDAA"

# --- 6) CANAL ---
entradaCanal = ""
salidaCanal = ""
probPriori = []                      # si te dan la matriz, ej: [0.5, 0.5]
matrizCanal = []                     #                     ej: [[0.9, 0.1], [0.1, 0.9]]
pasoCapacidad = 0.001
matrizCanal2 = []                    # segundo canal, solo si piden canales en serie

if __name__ == "__main__":

    # 1) FUENTE -> CODIGO OPTIMO -> COMPRESION
    if mensajeBits != "":
        alfabeto = alfabetoDado if len(alfabetoDado) > 0 else [f"S{i + 1}" for i in range(len(codigoDado))]
        simbolos = decodificar_bits(mensajeBits, alfabeto, codigoDado)
        print("Mensaje decodificado:", " ".join(simbolos))
        frecs = [simbolos.count(s) for s in alfabeto]
        probs = frecuencias_a_probabilidades(frecs, len(simbolos))
        mensaje = simbolos
    elif mensaje != "":
        alfabeto, frecs = alfabeto_y_frecuencias(mensaje)
        probs = frecuencias_a_probabilidades(frecs, len(mensaje))

    if mensajeBits != "" or mensaje != "":
        informe_fuente(alfabeto, probs, ordenExt)
        codigos = informe_codigo_optimo(alfabeto, probs, 1, metodo)
        if ordenExt > 1:
            informe_codigo_optimo(alfabeto, probs, ordenExt, metodo)
        informe_compresion(alfabeto, codigos, mensaje)

    # 2) CODIGO DADO
    if len(codigoAnalizar) > 0:
        informe_codigo_dado(probsCodigo, codigoAnalizar, ordenCodigo)

    # 3) HAMMING
    if len(codigoBloque) > 0:
        informe_hamming(codigoBloque)

    # 4) PARIDAD
    if textoParidad != "" or len(bytesRecibidos) > 0:
        informe_paridad(textoParidad, errores, bytesRecibidos)

    # 5) RLC
    if mensajeRLC != "":
        informe_rlc(mensajeRLC)

    # 6) CANAL
    if entradaCanal != "":
        matrizCanal = genero_matriz(entradaCanal, salidaCanal)
        _, alfEntrada = calcular_alfabeto(entradaCanal)
        _, alfSalida = calcular_alfabeto(salidaCanal)
        frecs = [entradaCanal.count(a) for a in alfEntrada]
        probPriori = frecuencias_a_probabilidades(frecs, len(entradaCanal))
        informe_canal(probPriori, matrizCanal, alfEntrada, alfSalida, pasoCapacidad)
    elif len(matrizCanal) > 0:
        alfEntrada = [f"a{i + 1}" for i in range(len(matrizCanal))]
        alfSalida = [f"b{j + 1}" for j in range(len(matrizCanal[0]))]
        informe_canal(probPriori, matrizCanal, alfEntrada, alfSalida, pasoCapacidad)

    if len(matrizCanal2) > 0:
        informe_canales_serie(matrizCanal, matrizCanal2)