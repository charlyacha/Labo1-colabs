"""
lab1_utils.py — funciones de análisis de Laboratorio 1
Departamento de Física, FCEN-UBA

Reúne las funciones que se construyen a lo largo de los Colabs 01 a 12.
No es una caja negra: cada función que está acá se escribió antes, a mano,
en el cuaderno de la clase correspondiente. El módulo existe para que en la
Práctica Especial no haya que reescribirlas.

Uso en Colab:

    !wget -q -O lab1_utils.py https://raw.githubusercontent.com/charlyacha/Labo1-colabs/main/lab1_utils.py
    import lab1_utils as lab

Uso local (Anaconda, Spyder, VS Code):

    poné este archivo en la misma carpeta que tu cuaderno y hacé
    import lab1_utils as lab
"""

from math import floor, log10

import numpy as np
import matplotlib.pyplot as plt
from scipy import stats
from scipy.optimize import curve_fit

__all__ = [
    "estilo_lab1", "guardar_figura",
    "redondear_con_error", "formatear", "reportar",
    "estadisticos", "bins_scott", "bins_freedman_diaconis", "chauvenet",
    "compatibilidad", "promedio_ponderado",
    "chi2_reducido", "matriz_correlacion", "R2",
    "ajustar", "grafico_con_residuos",
    "leer_datos", "leer_tracker",
]


# ============================================================ presentación

def estilo_lab1():
    """Fija tamaños de fuente y resolución razonables para informes.

    Los valores por defecto de matplotlib producen figuras cuya tipografía
    resulta ilegible al insertarlas en un informe o proyectarlas.
    """
    plt.rcParams.update({
        "figure.figsize": (7, 4.5),
        "figure.dpi": 110,
        "savefig.dpi": 200,
        "savefig.bbox": "tight",
        "font.size": 12,
        "axes.labelsize": 13,
        "axes.titlesize": 13,
        "legend.fontsize": 11,
        "xtick.labelsize": 11,
        "ytick.labelsize": 11,
        "axes.grid": True,
        "grid.alpha": 0.3,
        "errorbar.capsize": 3,
    })


def guardar_figura(nombre, fig=None):
    """Guarda la figura en PNG (200 dpi) y en PDF vectorial.

    El PDF es el que conviene insertar en el informe: no se pixela al
    ampliarlo ni al imprimirlo.
    """
    fig = fig or plt.gcf()
    fig.savefig(f"{nombre}.png", dpi=200, bbox_inches="tight")
    fig.savefig(f"{nombre}.pdf", bbox_inches="tight")
    print(f"Guardado: {nombre}.png y {nombre}.pdf")


# ================================================================ reporte

def redondear_con_error(x, dx, dos_cifras_si_empieza_bajo=True):
    """Redondea el par (valor, incerteza) según la convención estándar.

    La incerteza se lleva a una cifra significativa, o a dos si la primera
    es 1 o 2 (donde redondear a una sola cambiaría la incerteza en más de
    un 25 %). El valor se redondea a la misma posición decimal.
    Referencia: Taylor secc. 2.5; GUM (JCGM 100:2008) secc. 7.2.6.

    Devuelve (valor_redondeado, incerteza_redondeada, decimales).
    """
    if dx <= 0:
        raise ValueError("La incerteza tiene que ser positiva.")
    orden = floor(log10(abs(dx)))
    primera = int(dx / 10**orden)
    cifras = 2 if (dos_cifras_si_empieza_bajo and primera in (1, 2)) else 1
    dec = -(orden - (cifras - 1))
    return round(x, dec), round(dx, dec), max(dec, 0)


def formatear(x, dx, unidad=""):
    """Devuelve el string '(12,73 ± 0,30) cm' a partir de valor e incerteza."""
    x_r, dx_r, dec = redondear_con_error(x, dx)
    texto = f"({x_r:.{dec}f} ± {dx_r:.{dec}f}) {unidad}".strip()
    return texto.replace(".", ",")


def reportar(x, dx, unidad="", nombre=""):
    """Imprime el resultado ya redondeado y lo devuelve como string."""
    texto = formatear(x, dx, unidad)
    print(f"{nombre + ' = ' if nombre else ''}{texto}")
    return texto


# =========================================================== una variable

def estadisticos(x, verbose=True):
    """Promedio, desviación estándar muestral (ddof=1) y error de la media.

    Devuelve (media, s, sem). El ddof=1 es deliberado: np.std usa por
    defecto ddof=0, que calcula la desviación estándar POBLACIONAL y
    subestima la dispersión cuando lo que se tiene es una muestra.
    """
    x = np.asarray(x, float)
    n = len(x)
    media = x.mean()
    s = np.std(x, ddof=1)
    sem = s / np.sqrt(n)
    if verbose:
        print(f"N   = {n}")
        print(f"x̄   = {media:.6g}")
        print(f"s   = {s:.6g}    (dispersión de UNA medición)")
        print(f"SEM = {sem:.6g}    (incerteza DEL PROMEDIO)")
    return media, s, sem


def bins_scott(x):
    """Cantidad de bins por la regla de Scott (Biometrika 66(3), 605, 1979)."""
    x = np.asarray(x, float)
    h = 3.49 * np.std(x, ddof=1) / len(x)**(1/3)
    return max(1, int(np.ceil((x.max() - x.min()) / h)))


def bins_freedman_diaconis(x):
    """Cantidad de bins por Freedman-Diaconis (robusta frente a datos anómalos)."""
    x = np.asarray(x, float)
    iqr = np.percentile(x, 75) - np.percentile(x, 25)
    h = 2 * iqr / len(x)**(1/3)
    return max(1, int(np.ceil((x.max() - x.min()) / h)))


def chauvenet(x, verbose=True):
    """Criterio de Chauvenet: marca los datos sospechosos de una muestra.

    Devuelve un array de booleanos, True para los datos que el criterio
    CONSERVA. Ojo: el criterio se fija ANTES de mirar los datos, se aplica
    UNA sola vez, y todo descarte se informa en el informe.
    """
    x = np.asarray(x, float)
    n = len(x)
    z = np.abs(x - x.mean()) / np.std(x, ddof=1)
    n_esperados = n * 2 * stats.norm.sf(z)
    conservar = n_esperados >= 0.5
    if verbose:
        malos = np.where(~conservar)[0]
        if len(malos) == 0:
            print("Chauvenet no marca ningún dato.")
        else:
            for i in malos:
                print(f"  índice {i}: x = {x[i]:.6g}  ({z[i]:.2f} s del promedio, "
                      f"esperados {n_esperados[i]:.2f} datos así en N = {n})")
    return conservar


# ============================================================= comparación

def compatibilidad(x1, dx1, x2, dx2, etiquetas=("medición 1", "medición 2"),
                   verbose=True):
    """Discrepancia normalizada z entre dos mediciones independientes."""
    z = abs(x1 - x2) / np.sqrt(dx1**2 + dx2**2)
    if verbose:
        if z < 1:
            v = "compatibles"
        elif z < 2:
            v = "compatibles dentro de lo esperable"
        elif z < 3:
            v = "en tensión, conviene revisar"
        else:
            v = "INCOMPATIBLES: hay un sistemático o una incerteza subestimada"
        print(f"{etiquetas[0]}: {x1:.6g} ± {dx1:.3g}")
        print(f"{etiquetas[1]}: {x2:.6g} ± {dx2:.3g}")
        print(f"z = {z:.2f}   ->   {v}")
    return z


def promedio_ponderado(x, sigma, verbose=True):
    """Combina mediciones independientes pesando por 1/sigma².

    Devuelve (valor, incerteza, chi2_de_consistencia). El último diagnostica
    si las mediciones que se están combinando son compatibles entre sí: si
    da mucho mayor que 1, combinarlas no tenía sentido.
    """
    x, sigma = np.asarray(x, float), np.asarray(sigma, float)
    w = 1 / sigma**2
    xp = np.sum(w * x) / np.sum(w)
    sp = 1 / np.sqrt(np.sum(w))
    c2r = np.sum(w * (x - xp)**2) / (len(x) - 1) if len(x) > 1 else np.nan
    if verbose:
        print(f"{formatear(xp, sp)}    (χ²_ν de consistencia = {c2r:.2f})")
    return xp, sp, c2r


# ======================================================= bondad de ajuste

def chi2_reducido(y, y_modelo, sigma, n_parametros, verbose=True):
    """Chi cuadrado reducido y su p-valor.

    sigma tiene que ser la incerteza EXPERIMENTAL de cada punto, estimada de
    forma independiente del ajuste. Si se usa una sigma sacada del propio
    ajuste, el indicador se vuelve circular y no diagnostica nada.
    """
    y, y_modelo, sigma = (np.asarray(v, float) for v in (y, y_modelo, sigma))
    chi2 = np.sum(((y - y_modelo) / sigma)**2)
    nu = len(y) - n_parametros
    c2r = chi2 / nu
    p = stats.chi2.sf(chi2, nu)
    if verbose:
        if p < 0.01:
            lect = "modelo inadecuado o incertezas subestimadas"
        elif p < 0.05:
            lect = "en tensión: mirá los residuos antes de darlo por bueno"
        elif p > 0.99:
            lect = "demasiado bueno: incertezas sobreestimadas (o datos retocados)"
        else:
            lect = "consistente con un buen ajuste"
        print(f"χ² = {chi2:.2f}    ν = {nu}    χ²_ν = {c2r:.3f}    p = {p:.4f}")
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
    un modelo correcto de uno incorrecto (Spiess y Neumeyer, BMC
    Pharmacology 10:6, 2010). Está acá para poder mostrar cómo falla.
    """
    y = np.asarray(y, float)
    return 1 - np.sum((y - y_modelo)**2) / np.sum((y - y.mean())**2)


# ================================================================= ajuste

def ajustar(modelo, x, y, yerr=None, p0=None, nombres=None, verbose=True,
            **kwargs):
    """curve_fit con las opciones correctas y un informe completo.

    Si se pasan barras de error, usa absolute_sigma=True. Ese es el valor
    que corresponde cuando las incertezas son magnitudes físicas y no meros
    pesos relativos. Con el default (False) scipy reescala la covarianza por
    el χ²_ν del propio ajuste y destruye el diagnóstico de bondad de ajuste.

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
                     for i in range(len(popt)) for j in range(i + 1, len(popt))
                     if abs(M[i, j]) > 0.5]
            if fuera:
                print("  parámetros correlacionados (|ρ| > 0,5):")
                for a, b, r in fuera:
                    print(f"    ρ({a}, {b}) = {r:+.2f}")
    return popt, perr, pcov


def grafico_con_residuos(x, y, modelo, popt, yerr=None,
                         normalizar_residuos=False,
                         xlabel="x", ylabel="y", titulo="",
                         etiqueta_modelo="ajuste"):
    """Figura estándar del curso: datos y modelo arriba, residuos abajo.

    Todo ajuste de Laboratorio 1 se informa con esta figura. El panel de
    residuos detecta modelos mal especificados que ningún estadístico único
    alcanza a resumir.
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
    fig.subplots_adjust(hspace=0.08, left=0.13, right=0.97, top=0.93,
                        bottom=0.11)
    return fig, (ax1, ax2)


# ================================================== lectura de archivos

def _detectar_formato(texto):
    """Adivina separador de columnas y separador decimal."""
    if "\t" in texto:
        sep = "\t"
    elif ";" in texto:
        sep = ";"
    else:
        sep = None
    if sep is None:
        coma_decimal = False
        sep = ","
    else:
        coma_decimal = "," in texto
    return sep, coma_decimal


def leer_datos(archivo, sep=None, coma_decimal=None, verbose=True):
    """Lee una tabla numérica de texto sin pelearse con el formato.

    Resuelve los tres problemas que rompen np.loadtxt en la práctica:
    encabezados de largo variable, separador desconocido y coma decimal
    (que es lo que exporta cualquier programa configurado en español).

    Devuelve (nombres_de_columna, matriz) con la matriz de forma (N, ncol).
    """
    with open(archivo, "r", encoding="utf-8", errors="replace") as f:
        lineas = f.read().splitlines()
    texto = "\n".join(lineas)
    sep_auto, coma_auto = _detectar_formato(texto)
    sep = sep or sep_auto
    coma_decimal = coma_auto if coma_decimal is None else coma_decimal

    def a_numeros(linea):
        campos = [c.strip() for c in linea.split(sep)]
        if not any(campos):
            return None
        salida = []
        for c in campos:
            if c == "":
                salida.append(np.nan)
                continue
            if coma_decimal:
                c = c.replace(",", ".")
            try:
                salida.append(float(c))
            except ValueError:
                return None
        return salida

    filas, nombres, ancho = [], None, None
    for linea in lineas:
        fila = a_numeros(linea)
        if fila is None:
            if not filas:
                nombres = [c.strip() for c in linea.split(sep) if c.strip()]
            continue
        if ancho is None:
            ancho = len(fila)
        if len(fila) == ancho:
            filas.append(fila)

    if not filas:
        raise ValueError(f"No encontré datos numéricos en {archivo}.")
    datos = np.array(filas, float)
    if nombres is None or len(nombres) != datos.shape[1]:
        nombres = [f"col{i}" for i in range(datos.shape[1])]
    if verbose:
        print(f"{archivo}: {datos.shape[0]} filas x {datos.shape[1]} columnas")
        print(f"  separador: {sep!r}   coma decimal: {coma_decimal}")
        print(f"  columnas: {nombres}")
    return nombres, datos


def leer_tracker(archivo, columnas="txy", verbose=True):
    """Lee un archivo exportado por Tracker.

    Tracker hace el seguimiento del punto adentro del programa y exporta una
    tabla numérica: acá no hay ningún análisis de imagen. Lo que sí trae
    problemas es el formato del archivo (coma decimal si el sistema está en
    español, encabezados de metadatos, NaN en los cuadros no trackeados).

    columnas="txy" devuelve (t, x, y) recortando los NaN de los extremos.
    columnas="todas" devuelve (nombres, matriz).

    ADVERTENCIA: si el archivo trae columnas de velocidad o aceleración
    calculadas por Tracker, son diferencias finitas y arrastran todo el
    ruido amplificado. Para determinar una aceleración conviene ajustar
    x(t), no promediar la columna a.
    """
    nombres, datos = leer_datos(archivo, verbose=verbose)
    if columnas == "todas":
        return nombres, datos
    if datos.shape[1] < 3:
        raise ValueError("Esperaba al menos tres columnas (t, x, y).")
    t, x, y = datos[:, 0], datos[:, 1], datos[:, 2]
    valido = ~(np.isnan(t) | np.isnan(x) | np.isnan(y))
    if valido.sum() == 0:
        raise ValueError("Todas las filas tienen NaN.")
    i0, i1 = np.where(valido)[0][[0, -1]]
    recorte = slice(i0, i1 + 1)
    if verbose and (i0 > 0 or i1 < len(t) - 1):
        print(f"  recorté {i0} filas al principio y "
              f"{len(t) - 1 - i1} al final por NaN")
    return t[recorte], x[recorte], y[recorte]
