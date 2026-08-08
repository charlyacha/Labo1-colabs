"""
lab1_utils.py — funciones de análisis de Laboratorio 1
Departamento de Física, FCEN-UBA

Reúne todas las funciones construidas a lo largo de los Colabs 01 a 12,
más las de presentación y lectura de archivos que usan los cuadernos 07, 11, 12,
S1 y la plantilla de Práctica Especial.
Uso en Colab:

    !wget -q https://raw.githubusercontent.com/charlyacha/Labo1-colabs/main/lab1_utils.py
    from lab1_utils import *

o simplemente pegando este archivo en una celda.
"""

from math import floor, log10

import numpy as np
import matplotlib.pyplot as plt
from scipy import stats
from scipy.optimize import curve_fit

__all__ = [
    # presentación
    "estilo_lab1", "guardar_figura",
    # reporte
    "redondear_con_error", "formatear", "reportar",
    # una variable
    "estadisticos", "bins_scott", "bins_freedman_diaconis",
    "sigma_resolucion", "incerteza_de_s", "chauvenet",
    # comparación
    "compatibilidad", "promedio_ponderado",
    # bondad de ajuste
    "chi2_reducido", "matriz_correlacion", "R2",
    # ajuste
    "ajustar", "grafico_con_residuos",
    # señales y propagación
    "tiempos_de_flanco", "periodo_por_ajuste", "propagar",
    # lectura de archivos
    "leer_datos", "leer_tracker",
]


# ---------------------------------------------------------------- reporte

def redondear_con_error(x, dx, dos_cifras_si_empieza_en_1_o_2=True):
    """Redondea (valor, error) según la convención estándar.

    El error se lleva a una cifra significativa (dos si empieza en 1) y el
    valor se redondea a la misma posición decimal.
    """
    if dx <= 0:
        raise ValueError("La incerteza debe ser positiva.")
    orden = floor(log10(abs(dx)))
    primera = int(dx / 10**orden)
    cifras = 2 if (dos_cifras_si_empieza_en_1_o_2 and primera in (1, 2)) else 1
    dec = -(orden - (cifras - 1))
    return round(x, dec), round(dx, dec), max(dec, 0)


def formatear(x, dx, unidad=""):
    """'(12,73 ± 0,30) cm' a partir de valor e incerteza."""
    x_r, dx_r, dec = redondear_con_error(x, dx)
    return f"({x_r:.{dec}f} ± {dx_r:.{dec}f}) {unidad}".strip()


# ------------------------------------------------------------ una variable

def estadisticos(x, verbose=True):
    """Promedio, desviación estándar muestral (ddof=1) y error de la media.

    Devuelve (media, s, sem). El ddof=1 es deliberado: np.std usa por defecto
    ddof=0, que calcula la desviación estándar POBLACIONAL y subestima la
    dispersión de una muestra experimental.
    """
    x = np.asarray(x, float)
    n = len(x)
    media = x.mean()
    s = np.std(x, ddof=1)
    sem = s / np.sqrt(n)
    if verbose:
        print(f"N = {n}")
        print(f"x̄   = {media:.6g}")
        print(f"s   = {s:.6g}   (dispersión de una medición individual)")
        print(f"SEM = {sem:.6g}   (incerteza del promedio)")
    return media, s, sem


def bins_scott(x):
    """Cantidad de bins por la regla de Scott (Biometrika 66(3), 605, 1979)."""
    x = np.asarray(x, float)
    h = 3.49 * np.std(x, ddof=1) / len(x)**(1/3)
    return max(1, int(np.ceil((x.max() - x.min()) / h)))


def bins_freedman_diaconis(x):
    """Cantidad de bins por Freedman-Diaconis (robusta frente a outliers)."""
    x = np.asarray(x, float)
    iqr = np.percentile(x, 75) - np.percentile(x, 25)
    h = 2 * iqr / len(x)**(1/3)
    return max(1, int(np.ceil((x.max() - x.min()) / h)))


# ------------------------------------------------------------ comparación

def compatibilidad(x1, dx1, x2, dx2, etiquetas=("medición 1", "medición 2"), verbose=True):
    """Discrepancia normalizada z entre dos mediciones independientes."""
    z = abs(x1 - x2) / np.sqrt(dx1**2 + dx2**2)
    if verbose:
        if z < 1:     v = "compatibles"
        elif z < 2:   v = "compatibles dentro de lo esperable"
        elif z < 3:   v = "en tensión — revisar"
        else:         v = "INCOMPATIBLES — sistemático o incerteza subestimada"
        print(f"{etiquetas[0]}: {x1:.6g} ± {dx1:.3g}")
        print(f"{etiquetas[1]}: {x2:.6g} ± {dx2:.3g}")
        print(f"z = {z:.2f}  ->  {v}")
    return z


def promedio_ponderado(x, sigma, verbose=True):
    """Combina mediciones independientes pesando por 1/σ².

    Devuelve (valor, incerteza, chi2_reducido_de_consistencia). El último
    diagnostica si las mediciones combinadas son compatibles entre sí:
    si es mucho mayor que 1, no tenía sentido combinarlas.
    """
    x, sigma = np.asarray(x, float), np.asarray(sigma, float)
    w = 1 / sigma**2
    xp = np.sum(w * x) / np.sum(w)
    sp = 1 / np.sqrt(np.sum(w))
    c2r = np.sum(w * (x - xp)**2) / (len(x) - 1) if len(x) > 1 else np.nan
    if verbose:
        print(f"{formatear(xp, sp)}   (χ²_ν de consistencia = {c2r:.2f})")
    return xp, sp, c2r


# --------------------------------------------------------- bondad de ajuste

def chi2_reducido(y, y_modelo, sigma, n_parametros, verbose=True):
    """χ² reducido y su p-valor.

    sigma debe ser la incerteza EXPERIMENTAL de cada punto, estimada de forma
    independiente del ajuste. Usar una sigma derivada del propio ajuste vuelve
    el indicador circular y sin valor diagnóstico.
    """
    y, y_modelo, sigma = (np.asarray(v, float) for v in (y, y_modelo, sigma))
    chi2 = np.sum(((y - y_modelo) / sigma)**2)
    nu = len(y) - n_parametros
    c2r = chi2 / nu
    p = stats.chi2.sf(chi2, nu)
    if verbose:
        if c2r > 2:      lect = "modelo inadecuado o incertezas subestimadas"
        elif c2r < 0.5:  lect = "incertezas probablemente sobreestimadas"
        else:            lect = "consistente con un buen ajuste"
        print(f"χ² = {chi2:.2f}   ν = {nu}   χ²_ν = {c2r:.3f}   p = {p:.3f}")
        print(f"  -> {lect}")
    return c2r, p


def matriz_correlacion(pcov):
    """Matriz de correlación de los parámetros a partir de la covarianza."""
    d = np.sqrt(np.diag(pcov))
    return pcov / np.outer(d, d)


def R2(y, y_modelo):
    """Coeficiente de determinación convencional.

    ADVERTENCIA: no es un indicador válido de bondad de ajuste fuera de una
    regresión lineal ordinaria. Depende del rango muestreado y no distingue
    un modelo correcto de uno incorrecto (ver Spiess & Neumeyer, BMC
    Pharmacology 10:6, 2010). Se incluye solo para poder mostrar su falla.
    """
    y = np.asarray(y, float)
    return 1 - np.sum((y - y_modelo)**2) / np.sum((y - y.mean())**2)


# -------------------------------------------------------------- ajuste

def ajustar(modelo, x, y, yerr=None, p0=None, nombres=None, verbose=True, **kwargs):
    """curve_fit con las opciones correctas y un informe completo.

    Si se pasan barras de error se usa absolute_sigma=True, que es lo que
    corresponde cuando las incertezas son físicas y no meros pesos relativos.

    Devuelve (popt, perr, pcov).
    """
    x, y = np.asarray(x, float), np.asarray(y, float)
    if yerr is None:
        popt, pcov = curve_fit(modelo, x, y, p0=p0, **kwargs)
    else:
        yerr = np.asarray(yerr, float)
        popt, pcov = curve_fit(modelo, x, y, p0=p0, sigma=yerr,
                               absolute_sigma=True, **kwargs)
    perr = np.sqrt(np.diag(pcov))

    if verbose:
        if nombres is None:
            nombres = [f"p{i}" for i in range(len(popt))]
        for n, v, e in zip(nombres, popt, perr):
            rel = f"({100*abs(e/v):.2f} %)" if v != 0 else ""
            print(f"{n:<12}{v:>13.6g} ± {e:<12.4g} {rel}")
        if yerr is not None:
            chi2_reducido(y, modelo(x, *popt), yerr, len(popt))
        if len(popt) > 1:
            M = matriz_correlacion(pcov)
            fuera = [(nombres[i], nombres[j], M[i, j])
                     for i in range(len(popt)) for j in range(i+1, len(popt))
                     if abs(M[i, j]) > 0.5]
            if fuera:
                print("  parámetros correlacionados (|ρ| > 0,5):")
                for a, b, r in fuera:
                    print(f"    ρ({a}, {b}) = {r:+.2f}")
    return popt, perr, pcov


def grafico_con_residuos(x, y, modelo, popt, yerr=None, normalizar_residuos=False,
                         xlabel="x", ylabel="y", titulo="", etiqueta_modelo="ajuste"):
    """Figura estándar del curso: datos y modelo arriba, residuos abajo.

    Todo ajuste de Laboratorio 1 se informa con esta figura. El panel de
    residuos detecta modelos mal especificados que ningún estadístico único
    resume.
    """
    x, y = np.asarray(x, float), np.asarray(y, float)
    fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(7, 5.6), sharex=True,
                                   gridspec_kw={"height_ratios": [3, 1]})
    xx = np.linspace(x.min(), x.max(), 500)

    if yerr is None:
        ax1.plot(x, y, "o", ms=5, label="datos")
    else:
        ax1.errorbar(x, y, yerr=yerr, fmt="o", ms=5, capsize=3, label="datos")
    ax1.plot(xx, modelo(xx, *popt), "crimson", lw=1.8, label=etiqueta_modelo)
    ax1.set_ylabel(ylabel)
    ax1.grid(alpha=0.3)
    ax1.legend()
    if titulo:
        ax1.set_title(titulo)

    res = y - modelo(x, *popt)
    if normalizar_residuos and yerr is not None:
        ax2.plot(x, res / np.asarray(yerr, float), "o", ms=5)
        ax2.set_ylabel("residuo / σ")
        for s in (-2, 2):
            ax2.axhline(s, color="gray", lw=0.7, ls=":")
    elif yerr is None:
        ax2.plot(x, res, "o", ms=5)
        ax2.set_ylabel("residuo")
    else:
        ax2.errorbar(x, res, yerr=yerr, fmt="o", ms=5, capsize=3)
        ax2.set_ylabel("residuo")
    ax2.axhline(0, color="crimson", lw=1.2)
    ax2.set_xlabel(xlabel)
    ax2.grid(alpha=0.3)
    fig.subplots_adjust(hspace=0.08, left=0.13, right=0.97, top=0.93, bottom=0.11)
    return fig, (ax1, ax2)


# ------------------------------------------------- agregados de la revisión 260805

def reportar(x, dx, unidad="", nombre=None):
    """Formatea un resultado como '(valor ± error) unidad'.

    Tiene dos modos, según cómo lo usan los cuadernos:

    - Sin ``nombre`` (Colabs 02–05): devuelve el string y no imprime nada, de
      modo que se puede intercalar dentro de un f-string o un print.
    - Con ``nombre`` (Colabs 10, 12): imprime la línea etiquetada
      ``nombre = (valor ± error) unidad`` y además devuelve el string.
    """
    s = formatear(x, dx, unidad)
    if nombre is not None:
        print(f"{nombre} = {s}")
    return s


def sigma_resolucion(delta):
    """Desviación estándar del error de cuantización de un instrumento de resolución delta.

    El error de resolución se distribuye uniformemente en un intervalo de ancho delta,
    cuya varianza es delta**2/12.  (Colab 01, sección 3; Colab 03, sección 7.)
    """
    return delta / np.sqrt(12)


def incerteza_de_s(m):
    """Incerteza RELATIVA de una desviación estándar estimada con m datos.

    sigma_s / s ~ 1/sqrt(2(m-1)).  Con m=3 es 50 %; con m=8, 27 %.  (Colab 03, sección 6.)
    """
    return 1.0 / np.sqrt(2 * (m - 1))


def chauvenet(x, verbose=True):
    """Máscara de los datos que SOBREVIVEN al criterio de Chauvenet.

    Se descarta el dato cuyo numero esperado de ocurrencias tan extremas, en una muestra
    de tamaño N, sea menor que 0,5.  Supone distribución normal.

    Reglas de conducta (Colab 03, sección 9):
      1. el criterio se fija antes de mirar los datos;
      2. se aplica una sola vez;
      3. todo descarte se informa.
    """
    from scipy.stats import norm
    x = np.asarray(x, float)
    n = len(x)
    z = np.abs(x - x.mean()) / np.std(x, ddof=1)
    esperados = n * 2 * norm.sf(z)
    mantener = esperados >= 0.5
    if verbose:
        for i in np.where(~mantener)[0]:
            print(f"  descartar x[{i}] = {x[i]:.6g}   z = {z[i]:.2f}   "
                  f"esperados = {esperados[i]:.3f}")
        print(f"  sobreviven {mantener.sum()} de {n} datos")
    return mantener


def tiempos_de_flanco(t, V, umbral, tipo="bajada"):
    """Instantes en que V(t) cruza un umbral, con interpolación lineal entre muestras.

    tipo: 'bajada' (la señal entra al nivel bajo) o 'subida'.
    Usar SIEMPRE el mismo tipo dentro de una misma medición.  (Colab 04, sección 2.)

    Nota: si el flanco es más rápido que el intervalo de muestreo, la interpolación no
    aporta información y el piso de incerteza sigue siendo sigma_resolucion(1/f_s).
    """
    t, V = np.asarray(t, float), np.asarray(V, float)
    b = (V < umbral).astype(int)
    d = np.diff(b)
    idx = np.where(d == (1 if tipo == "bajada" else -1))[0]
    V0, V1 = V[idx], V[idx + 1]
    t0, t1 = t[idx], t[idx + 1]
    return t0 + (umbral - V0) / (V1 - V0) * (t1 - t0)


def periodo_por_ajuste(t_eventos, eventos_por_periodo=1, sigma_t=None, verbose=True):
    """Período por ajuste de t_n vs. n, que usa TODOS los eventos.

    Promediar diferencias consecutivas telescopea y usa solo el primero y el último.
    (Colab 04, sección 2.3; Colab 09, sección 4.)

    eventos_por_periodo: 2 para un péndulo que corta el haz de ida y de vuelta.
    """
    t_eventos = np.asarray(t_eventos, float)[::eventos_por_periodo]
    n = np.arange(len(t_eventos))
    if sigma_t is None:
        popt, pcov = curve_fit(lambda k, T, t0: T * k + t0, n, t_eventos)
    else:
        popt, pcov = curve_fit(lambda k, T, t0: T * k + t0, n, t_eventos,
                               sigma=np.full(len(n), float(sigma_t)),
                               absolute_sigma=True)
    T, dT = popt[0], np.sqrt(pcov[0, 0])
    if verbose:
        print(f"T = {formatear(T, dT, 's')}   (a partir de {len(t_eventos)} eventos)")
    return T, dT


def propagar(expr, variables, valores, errores, nombre="f", verbose=True):
    """Propagación de primer orden sobre una expresión de SymPy, con tabla de contribuciones.

    expr      : expresión de sympy
    variables : lista de símbolos
    valores   : dict {símbolo: valor}
    errores   : dict {símbolo: incerteza}

    La tabla de contribuciones es el resultado útil: dice qué medición conviene mejorar.
    (Colab 02, sección 4.)
    """
    import sympy as sp
    valor = float(expr.subs(valores))
    contrib = {}
    for v in variables:
        d = float(sp.diff(expr, v).subs(valores))
        contrib[v] = (d * errores[v]) ** 2
    var_total = sum(contrib.values())
    sigma = np.sqrt(var_total)
    if verbose:
        print(f"{nombre} = {valor:.6g}  ±  {sigma:.3g}   "
              f"({100*sigma/abs(valor):.2f} %)")
        print(f"{'variable':>10s} {'σ rel.':>9s} {'contrib. a σ²':>15s}")
        print("-" * 38)
        for v in variables:
            print(f"{str(v):>10s} {100*errores[v]/float(valores[v]):8.2f}% "
                  f"{100*contrib[v]/var_total:14.1f}%")
    return valor, sigma


# ============================================================ presentación
#
# Estas dos funciones no existían en la primera versión del módulo y las
# invocaban seis cuadernos (07, 10, 11, 12, S1 y la plantilla PE), que por lo
# tanto fallaban en su primera celda. Se agregaron para cerrar esa brecha.

def estilo_lab1():
    """Estilo de figura uniforme para todo el curso.

    Ajusta los parámetros de Matplotlib para que todas las figuras salgan con
    el mismo tamaño de fuente, grilla tenue y resolución razonable. Se llama
    una vez, al principio del cuaderno, justo después de importar el módulo.
    No devuelve nada: modifica el estado global de Matplotlib.
    """
    plt.rcParams.update({
        "figure.figsize": (7, 4.5),
        "figure.dpi": 110,
        "savefig.dpi": 200,
        "font.size": 12,
        "axes.titlesize": 13,
        "axes.labelsize": 12,
        "axes.grid": True,
        "grid.alpha": 0.3,
        "grid.linewidth": 0.6,
        "lines.markersize": 5,
        "errorbar.capsize": 3,
        "legend.frameon": False,
        "axes.spines.top": False,
        "axes.spines.right": False,
    })


def guardar_figura(nombre, fig=None, formato="pdf"):
    """Guarda la figura actual como archivo vectorial, listo para el informe.

    Por defecto exporta en PDF, que es vectorial y no pixela al ampliarlo — la
    forma correcta de llevar una figura a un informe. Se le pasa el nombre sin
    extensión. Si no se le pasa una figura explícita, guarda la activa.
    """
    if fig is None:
        fig = plt.gcf()
    archivo = f"{nombre}.{formato}"
    fig.savefig(archivo, bbox_inches="tight")
    print(f"Figura guardada: {archivo}")
    return archivo


# ======================================================= lectura de archivos
#
# leer_datos es un lector robusto de tablas "sucias" reales: descarta líneas de
# comentario y de metadatos, ubica la fila de encabezado, autodetecta el
# separador (tab, ';', ',' o espacios) y la coma decimal, y convierte los
# faltantes (NaN, celdas vacías) en np.nan sin fallar. leer_tracker es la
# especialización para los .txt que exporta Tracker. (Colab 07.)

def _detectar_separador(linea):
    """Elige el separador más probable de una línea de datos."""
    for sep in ["\t", ";"]:
        if sep in linea:
            return sep
    # coma como separador solo si hay más de una y no es la coma decimal:
    # heurística — si hay comas rodeadas de dígitos a ambos lados en un patrón
    # que se repite, es decimal; si separa campos, hay comas seguidas de espacio
    # o de texto. Se resuelve con el conteo relativo en _leer_core.
    if "," in linea:
        return ","
    return None  # espacios en blanco


def _leer_core(archivo):
    """Motor común de lectura. Devuelve (nombres, datos_2D)."""
    with open(archivo, "r", encoding="utf-8", errors="replace") as f:
        crudas = [ln.rstrip("\n") for ln in f]

    # 1) descartar comentarios y líneas vacías
    lineas = [ln for ln in crudas
              if ln.strip() and not ln.lstrip().startswith("#")]
    if not lineas:
        raise ValueError(f"{archivo}: no hay datos después de quitar comentarios.")

    # 2) ubicar el encabezado: primera línea que contiene letras (nombres de
    #    columna) en lugar de solo números/separadores
    def es_numerica(ln):
        # para CLASIFICAR una línea (dato vs. encabezado) tratamos coma,
        # punto y coma y tabulador como separadores. Así '0,0000' (coma
        # decimal) da los tokens '0' y '0000' —ambos numéricos— y '0.0,0.02'
        # (coma separadora) da '0.0' y '0.02'; en los dos casos la línea se
        # reconoce como de datos. Una línea con letras (encabezado) no.
        for ch in (",", ";", "\t"):
            ln = ln.replace(ch, " ")
        toks = ln.split()
        if not toks:
            return False
        for tok in toks:
            try:
                float(tok)
            except ValueError:
                if tok.lower() not in ("nan", "inf", "-inf"):
                    return False
        return True

    # localizar la primera fila de datos: todo lo anterior es metadatos +
    # encabezado. Los archivos de Tracker traen dos líneas de metadatos
    # ('Tracker 6.1.4', 'masa A') antes del encabezado real, así que hay que
    # saltar TODAS las líneas no numéricas iniciales, no solo una.
    idx_datos = 0
    while idx_datos < len(lineas) and not es_numerica(lineas[idx_datos]):
        idx_datos += 1
    if idx_datos == len(lineas):
        raise ValueError(f"{archivo}: no se encontró ninguna fila de datos numéricos.")

    nombres = None
    if idx_datos > 0:
        # el encabezado es la última línea no numérica antes de los datos;
        # las anteriores (si las hay) son metadatos y se descartan
        cab = lineas[idx_datos - 1]
        sep_h = _detectar_separador(cab)
        nombres = [c.strip() for c in (cab.split(sep_h) if sep_h else cab.split())]

    cuerpo = lineas[idx_datos:]
    if not cuerpo:
        raise ValueError(f"{archivo}: encabezado sin filas de datos.")

    # 3) separador de los datos, mirando la primera fila de datos
    sep = _detectar_separador(cuerpo[0])

    # 4) coma decimal: si el separador NO es ',', pero aparecen comas dentro de
    #    los campos, entonces la coma es decimal y hay que reemplazarla.
    coma_decimal = (sep != ",") and any("," in ln for ln in cuerpo)

    filas = []
    for ln in cuerpo:
        s = ln
        if coma_decimal:
            s = s.replace(",", ".")
        campos = s.split(sep) if sep else s.split()
        fila = []
        for c in campos:
            c = c.strip()
            if c == "" or c.lower() == "nan":
                fila.append(np.nan)
            else:
                try:
                    fila.append(float(c))
                except ValueError:
                    fila.append(np.nan)
        filas.append(fila)

    ancho = max(len(f) for f in filas)
    for f in filas:
        f += [np.nan] * (ancho - len(f))
    datos = np.array(filas, dtype=float)

    if nombres is None:
        nombres = [f"col{i}" for i in range(ancho)]
    return nombres, datos


def leer_datos(archivo):
    """Lector robusto de una tabla de datos con formato imperfecto.

    Descarta líneas de comentario (que empiezan con '#') y de metadatos,
    ubica la fila de encabezado, autodetecta el separador (tabulador, ';',
    ',' o espacios) y la coma decimal, y convierte los faltantes en np.nan.

    Devuelve (nombres, datos), con nombres la lista de columnas y datos un
    arreglo 2D de floats. (Colab 07.)
    """
    return _leer_core(archivo)


def leer_tracker(archivo):
    """Lee un .txt exportado por Tracker y devuelve (t, x, y).

    Tolera las tres patologías típicas de esos archivos: las dos primeras
    líneas de metadatos (versión y nombre de la masa), la coma decimal, y los
    'NaN' de los cuadros no trackeados. Los NaN se conservan a propósito: el
    recorte del tramo útil con máscaras booleanas es parte del ejercicio del
    Colab 07, no algo que la función deba resolver sola.

    Devuelve las tres primeras columnas como arreglos separados.
    """
    _, datos = _leer_core(archivo)
    if datos.shape[1] < 3:
        raise ValueError(
            f"{archivo}: se esperaban al menos 3 columnas (t, x, y), "
            f"se encontraron {datos.shape[1]}.")
    return datos[:, 0], datos[:, 1], datos[:, 2]
