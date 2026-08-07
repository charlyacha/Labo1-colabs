"""
lab1_utils.py — funciones de análisis de Laboratorio 1
Departamento de Física, FCEN-UBA

Reúne las funciones construidas a lo largo de los Colabs 01 a 10.
Uso en Colab:

    !wget -q https://raw.githubusercontent.com/TU-USUARIO/labo1-colabs/main/lab1_utils.py
    from lab1_utils import *

o simplemente pegando este archivo en una celda.
"""

from math import floor, log10

import numpy as np
import matplotlib.pyplot as plt
from scipy import stats
from scipy.optimize import curve_fit

__all__ = [
    "redondear_con_error", "formatear",
    "estadisticos", "bins_scott", "bins_freedman_diaconis",
    "compatibilidad", "promedio_ponderado",
    "chi2_reducido", "matriz_correlacion", "R2",
    "ajustar", "grafico_con_residuos",
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

def reportar(x, dx, unidad=""):
    """Alias de formatear(). Nombre usado en los Colabs 02 en adelante."""
    return formatear(x, dx, unidad)


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
