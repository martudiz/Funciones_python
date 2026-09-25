"""
Funciones_guia2.py

Archivo autocontenido (sin imports locales) que resuelve TODOS los ejercicios de la
guía 2 (códigos) a partir de una única entrada: el código (palabras código) y las
probabilidades de cada símbolo.

Con esa entrada se resuelven los ejercicios 6, 9, 11, 14 y 16 de la guía.

Parámetro secundario (no forma parte del código/probabilidades, pero hace falta para
el ejercicio 16): la cantidad N de símbolos a generar, editable en main().
"""
from math import log, ceil
import random


# =============================================================================
# Ejercicio 6: no singular / instantáneo / unívocamente decodificable
# =============================================================================

def no_singular(codigo):
    """
    Entra: 'codigo', lista de strings con las palabras código.
    Sale: (es_no_singular, conjunto), donde es_no_singular es True si todas las
          palabras son distintas entre sí, y 'conjunto' es el set de palabras
          código (se reutiliza en es_univocamente_decodificable()).
    """
    conjunto = set(codigo)
    return len(conjunto) == len(codigo), conjunto


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


def es_univocamente_decodificable(codigo, conjunto=None):
    """
    Entra: 'codigo' (lista de palabras código, en cualquier orden) y 'conjunto'
           (opcional; el set de esas palabras, ver no_singular()).
    Sale: True si el código es unívocamente decodificable, False si no.

    Implementa el algoritmo de Sardinas-Patterson completo (no solo la primera
    "generación"): arma el conjunto de sufijos pendientes comparando las
    palabras código entre sí (generación 1); si algún sufijo coincide con una
    palabra código, hay conflicto. Si no, arma la siguiente generación de
    sufijos comparando los sufijos ya encontrados contra el código, y repite.
    Si una generación se repite (ciclo) sin haber encontrado nunca un
    conflicto, el código es unívocamente decodificable.

    (Esto reemplaza la versión anterior, que solo evaluaba la primera
    generación: eso hacía que códigos como {0, 01, 10} -que sí son ambiguos,
    porque "010" se puede leer como '01'+'0' o como '0'+'10'- dieran
    incorrectamente 'True'. Con las generaciones siguientes, el conflicto en
    ese ejemplo aparece recién en la segunda: comparando el sufijo '1' contra
    '10' se obtiene el sufijo '0', que sí es una palabra del código.)
    """
    if conjunto is None:
        conjunto = set(codigo)

    # Generación 1: sufijos pendientes de comparar cada palabra código con las
    # demás palabras código de las que es prefijo.
    generacion = set()
    for a in codigo:
        for b in codigo:
            if a != b and b.startswith(a):
                generacion.add(b[len(a):])
    if any(sufijo in conjunto for sufijo in generacion):
        return False

    generaciones_vistas = [frozenset(generacion)]
    while generacion:
        siguiente = set()
        for sufijo in generacion:
            for c in codigo:
                if sufijo.startswith(c):
                    siguiente.add(sufijo[len(c):])
                elif c.startswith(sufijo):
                    siguiente.add(c[len(sufijo):])
        if any(s in conjunto for s in siguiente):
            return False
        if frozenset(siguiente) in generaciones_vistas:
            return True
        generaciones_vistas.append(frozenset(siguiente))
        generacion = siguiente
    return True


def analizar_codigo(codigo):
    """
    Entra: 'codigo', lista de strings con las palabras código (en cualquier orden).
    Sale: diccionario con 'ordenado' (código ordenado), 'no_singular',
          'instantaneo' y 'univocamente_decodificable' (estos dos últimos en
          None si el código es singular, porque no tiene sentido evaluarlos).
          Combina los incisos a, b y c del ejercicio 6 en una sola llamada.
    """
    codigo_ordenado = sorted(codigo)
    singular_ok, conjunto = no_singular(codigo)
    resultado = {
        "ordenado": codigo_ordenado,
        "no_singular": singular_ok,
        "instantaneo": None,
        "univocamente_decodificable": None,
    }
    if singular_ok:
        resultado["instantaneo"] = es_instantaneo(codigo_ordenado)
        resultado["univocamente_decodificable"] = es_univocamente_decodificable(
            codigo_ordenado, conjunto
        )
    return resultado


# =============================================================================
# Ejercicio 9: alfabeto código, longitudes y desigualdad de Kraft
# =============================================================================

def alfabeto_desde_codigo(codigo):
    """
    Entra: 'codigo', lista de strings con las palabras código.
    Sale: lista con los símbolos únicos usados en el código (el "alfabeto
          código"), en orden de primera aparición.
    """
    alfabeto = []
    for palabra in codigo:
        for caracter in palabra:
            if caracter not in alfabeto:
                alfabeto.append(caracter)
    return alfabeto


def longitudes_palabras(codigo):
    """
    Entra: 'codigo', lista de strings con las palabras código.
    Sale: lista paralela con la longitud de cada palabra.
    """
    return [len(palabra) for palabra in codigo]


def kraft(longitudes, r):
    """
    Entra: 'longitudes' (lista de longitudes de las palabras código) y 'r' (el
           tamaño del alfabeto código, radix).
    Sale: la sumatoria de la inecuación de Kraft: sum(r^(-l)) para cada
          longitud l. Si el resultado es <= 1, el código puede ser
          instantáneo; si es > 1, no puede serlo.
    """
    return sum(r ** (-l) for l in longitudes)


# =============================================================================
# Ejercicio 11: entropía de la fuente y longitud media del código
# =============================================================================

def informacion_bits(probabilidades, base=2):
    """
    Entra: 'probabilidades' (lista de floats en (0,1], suman 1) y 'base' (base
           del logaritmo; para códigos se usa r = tamaño del alfabeto código).
    Sale: lista paralela con la información de cada símbolo: I = log_base(1/p).
    """
    return [log(1 / p, base) for p in probabilidades]


def entropia(probabilidades, informacion=None):
    """
    Entra: 'probabilidades' (lista de floats que suman 1) e 'informacion'
           (opcional; si no se pasa, se calcula con informacion_bits()).
    Sale: la entropía de la fuente: H = sum(p * I).
    """
    if informacion is None:
        informacion = informacion_bits(probabilidades)
    return sum(p * i for p, i in zip(probabilidades, informacion))


def longitud_media(longitudes, probabilidades):
    """
    Entra: 'longitudes' y 'probabilidades', listas paralelas (longitud de cada
           palabra código y probabilidad de cada símbolo, suman 1).
    Sale: la longitud media del código: sum(longitud_i * probabilidad_i).
    """
    return sum(l * p for l, p in zip(longitudes, probabilidades))


# =============================================================================
# Ejercicio 14: código compacto
# =============================================================================

def longitudes_ideales(probabilidades, r):
    """
    Entra: 'probabilidades' (suman 1) y 'r' (tamaño del alfabeto código).
    Sale: lista paralela con la longitud ideal (de Shannon) de cada símbolo,
          redondeada hacia arriba: techo(log_r(1/p)).
    Precondición: r >= 2 (un alfabeto código de un solo símbolo no permite
    calcular longitudes de Shannon, porque el logaritmo en base 1 no existe).
    """
    if r < 2:
        raise ValueError(
            f"El alfabeto código debe tener al menos 2 símbolos (r >= 2); se recibió r={r}."
        )
    return [ceil(i) for i in informacion_bits(probabilidades, base=r)]


def es_compacto(codigo, probabilidades):
    """
    Entra: 'codigo' y 'probabilidades', listas paralelas (palabras código y su
           probabilidad, suman 1). El tamaño del alfabeto código (r) se infiere
           de los símbolos usados en 'codigo'.
    Sale: True si el código es compacto: es unívocamente decodificable y
          ninguna palabra es más larga que su longitud ideal de Shannon.
          False si el código es singular, no es unívocamente decodificable, o
          alguna palabra excede su longitud ideal.
    """
    singular_ok, conjunto = no_singular(codigo)
    if not singular_ok:
        return False
    codigo_ordenado = sorted(codigo)
    if not es_univocamente_decodificable(codigo_ordenado, conjunto):
        return False
    r = len(alfabeto_desde_codigo(codigo))
    longitudes = longitudes_palabras(codigo)
    ideales = longitudes_ideales(probabilidades, r)
    return all(l <= li for l, li in zip(longitudes, ideales))


# =============================================================================
# Ejercicio 16: generar un mensaje aleatorio codificado
# =============================================================================

def generar_mensaje_codificado(n, codigo, probabilidades):
    """
    Entra:
      - n: entero > 0, cantidad de símbolos a generar.
      - codigo: lista con las palabras código de la fuente (strings).
      - probabilidades: lista paralela a 'codigo', con la probabilidad de cada
        palabra código (deben sumar 1).
    Sale:
      - un string con los 'n' símbolos generados al azar (respetando las
        probabilidades dadas) ya codificados y concatenados uno atrás del otro,
        simulando un mensaje codificado emitido por la fuente.
    """
    frecuencia_acumulada = []
    acumulado = 0
    for p in probabilidades:
        acumulado += p
        frecuencia_acumulada.append(acumulado)

    mensaje_codificado = ""
    for _ in range(n):
        r = random.random()
        indice = 0
        while indice < len(frecuencia_acumulada) - 1 and r > frecuencia_acumulada[indice]:
            indice += 1
        mensaje_codificado += codigo[indice]

    return mensaje_codificado


# =============================================================================
# MAIN: resuelve todos los ejercicios a partir de un único código + probabilidades
# =============================================================================

def main():
    # ------------------------------------------------------------------
    # --- ÚNICA ENTRADA: el código y la probabilidad de cada palabra ---
    codigo = ['0', '10', '110', '111']
    probabilidades = [0.4, 0.3, 0.2, 0.1]

    # --- Parámetro secundario (editable) ---
    n_a_generar = 15   # cantidad de símbolos a generar (ej. 16)
    # ------------------------------------------------------------------

    print("=" * 70)
    print("CÓDIGO DE ENTRADA:", codigo)
    print("PROBABILIDADES:", probabilidades)
    print("=" * 70)

    # --- Ejercicio 6: no singular, instantáneo, unívocamente decodificable ---
    print("\n--- Ejercicio 6: no singular / instantáneo / unívocamente decodificable ---")
    resultado_6 = analizar_codigo(codigo)
    print("Código ordenado:", resultado_6["ordenado"])
    print("Inciso A - Es no singular:", resultado_6["no_singular"])
    if resultado_6["no_singular"]:
        print("Inciso B - Es instantáneo:", resultado_6["instantaneo"])
        print("Inciso C - Es unívocamente decodificable:", resultado_6["univocamente_decodificable"])
    else:
        print("Al ser singular, no tiene sentido evaluar si es instantáneo o unívocamente decodificable.")

    # --- Ejercicio 9: alfabeto código, longitudes y desigualdad de Kraft ---
    print("\n--- Ejercicio 9: alfabeto código, longitudes y desigualdad de Kraft ---")
    alfabeto = alfabeto_desde_codigo(codigo)
    print("Inciso A - Alfabeto código:", sorted(alfabeto))
    longitudes = longitudes_palabras(codigo)
    print("Inciso B - Longitudes:", longitudes)
    r = len(alfabeto)
    print("Inciso C - Sumatoria de Kraft:", kraft(longitudes, r))

    # --- Ejercicio 11: entropía de la fuente y longitud media del código ---
    print("\n--- Ejercicio 11: entropía de la fuente y longitud media ---")
    info = informacion_bits(probabilidades, base=r)
    print("Información por símbolo:", info)
    print("Inciso A - Entropía de la fuente:", entropia(probabilidades, info))
    print("Inciso B - Longitud media del código:", longitud_media(longitudes, probabilidades))

    # --- Ejercicio 14: código compacto ---
    print("\n--- Ejercicio 14: ¿es compacto? ---")
    print("¿Es compacto?:", es_compacto(codigo, probabilidades))

    # --- Ejercicio 16: generar un mensaje codificado al azar ---
    print(f"\n--- Ejercicio 16: mensaje codificado de {n_a_generar} símbolos ---")
    mensaje = generar_mensaje_codificado(n_a_generar, codigo, probabilidades)
    print("Mensaje codificado:", mensaje)


if __name__ == "__main__":
    main()