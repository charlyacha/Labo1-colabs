# Laboratorio 1 — Cuadernos de Colab

**Departamento de Física, FCEN-UBA · Cátedra Acha · 2.º cuatrimestre 2026**

Repositorio de los cuadernos de análisis de datos del curso. Cada cuaderno corresponde a una clase y se abre directamente en Google Colab desde el badge que trae arriba.

---

## Contenido

| Archivo | Clase | Pregunta / contenido |
|---|---|---|
| `01_Primeros_pasos_medir_y_reportar.ipynb` | 1 | ¿Por qué no medimos todos lo mismo? |
| `02_Estadistica_de_una_variable.ipynb` | 2 | ¿Cuánto tardo en reaccionar, y qué número reporto como incerteza? |
| `02b_Estadistica_de_una_variable_INCOMPLETO.ipynb` | 2 | Versión con celdas a completar en clase |
| `03_Gaussiana_TCL_y_compatibilidad.ipynb` | 3 | ¿Cuántas veces vale la pena medir? |
| `04_Propagacion_de_incertezas.ipynb` | 4 | ¿De qué material está hecho este cuerpo? |
| `05_Cuadrados_minimos_y_el_coeficiente_R.ipynb` | 5 | ¿Este resorte cumple la ley de Hooke, y hasta dónde? |
| `06_Ajuste_ponderado_covarianza_y_chi2.ipynb` | 6 | ¿Cambia *k* si tomo en serio las barras de error? |
| `07_Datos_reales_adquisicion_y_derivadas.ipynb` | 7 | ¿Tres instrumentos distintos miden la misma velocidad? |
| `08_Determinacion_de_g_y_comparacion_de_modelos.ipynb` | 8 | ¿Puedo determinar el valor de *g* en este laboratorio? |
| `09_Senales_periodicas_y_determinacion_del_periodo.ipynb` | 9 | ¿Puedo medir una masa desconocida sin usar una balanza? |
| `10_Ajuste_no_lineal_y_por_que_R2_miente.ipynb` | 10 | ¿Cuánto tarda en morirse esta oscilación? |
| `11_Viscosidad_por_caida_de_esferas.ipynb` | 11 | ¿Puedo medir la viscosidad de un fluido dejando caer una esfera? |
| `12_Choques_y_conservacion.ipynb` | 12 | ¿Se conserva realmente el momento en un choque? |
| `S1_El_experimento_repetido_y_la_distribucion_de_chi2.ipynb` | — | Optativo: de dónde sale la distribución de χ² |
| `PE_Plantilla_analisis_Practica_Especial.ipynb` | 14–16 | Plantilla de análisis para la Práctica Especial |
| `lab1_utils.py` | — | Módulo con las funciones construidas a lo largo del curso |

Tres cuadernos tienen celdas incompletas **a propósito** y por lo tanto no corren de punta a punta: `02b`, `11` y la plantilla `PE`. Los otros doce sí.

---

## Cómo publicarlo

El repositorio previsto es `https://github.com/charlyacha/Labo1-colabs`, rama `main`. Los badges de los cuadernos ya apuntan ahí.

```bash
git clone https://github.com/charlyacha/Labo1-colabs.git
cd Labo1-colabs
cp /ruta/a/los/archivos/*.ipynb .
cp /ruta/a/los/archivos/lab1_utils.py .
cp /ruta/a/los/archivos/README.md .
git add .
git commit -m "Cuadernos del 2.º cuatrimestre 2026"
git push
```

Si en algún momento cambiás el nombre del usuario o del repositorio, hay que actualizar dos cosas en cada cuaderno: la URL del badge y la URL del `wget` de la celda de preparación. Se hace con una línea:

```bash
sed -i 's|charlyacha/Labo1-colabs|NUEVO-USUARIO/NUEVO-REPO|g' *.ipynb
```

### Por qué GitHub y no Drive

Con GitHub, el badge "Abrir en Colab" genera una copia limpia para cada estudiante sin que puedan tocar el original, y vos tenés historial de versiones de cada cambio. Con un Drive compartido, o das permiso de edición —y alguien rompe el cuaderno para todos— o das solo lectura y cada estudiante tiene que hacer el "Guardar una copia" manualmente igual. GitHub resuelve las dos cosas y además te deja corregir un error a mitad de cuatrimestre sin volver a repartir links.

---

## Instrucciones para los estudiantes

1. Abrí el cuaderno con el badge **Abrir en Colab**.
2. **Archivo → Guardar una copia en Drive.** Esto es lo primero, antes de escribir nada: lo que hagas sin copiar primero no se guarda en ningún lado.
3. Ejecutá la celda de preparación (la primera con código). Baja `lab1_utils.py` automáticamente.
4. Recorré el cuaderno **de arriba hacia abajo**. Las celdas no son independientes: si algo tira `NameError`, casi siempre es porque salteaste una.
5. Para subir tus datos, usá el ícono de carpeta del panel izquierdo.
6. Antes de entregar: **Entorno de ejecución → Reiniciar y ejecutar todo**. Ése, y no otro, es el estado en el que se entrega un cuaderno.

---

## El módulo `lab1_utils.py`

No es una caja negra. Cada función que está ahí se escribe antes, a mano, en el cuaderno de la clase correspondiente; el módulo existe para no tener que reescribirlas en la Práctica Especial.

**Presentación:** `estilo_lab1()`, `guardar_figura()`
**Reporte:** `redondear_con_error()`, `formatear()`, `reportar()`
**Una variable:** `estadisticos()`, `bins_scott()`, `bins_freedman_diaconis()`, `chauvenet()`
**Comparación:** `compatibilidad()`, `promedio_ponderado()`
**Bondad de ajuste:** `chi2_reducido()`, `matriz_correlacion()`, `R2()`
**Ajuste:** `ajustar()`, `grafico_con_residuos()`
**Archivos:** `leer_datos()`, `leer_tracker()`

Tres decisiones de diseño que conviene conocer, porque contradicen los valores por defecto de las librerías:

- `estadisticos()` usa **`ddof=1`**. El default de `np.std` es `ddof=0`, que calcula la desviación estándar poblacional y subestima la dispersión de una muestra. No da error: da un número plausible y equivocado.
- `ajustar()` usa **`absolute_sigma=True`** siempre que se pasen barras de error. El default de `curve_fit` es `False`, que reescala la covarianza por el χ²ᵥ del propio ajuste y destruye el diagnóstico de bondad de ajuste por construcción. La trampa gemela está en `np.polyfit(..., cov=True)`, que necesita `cov='unscaled'`.
- `R2()` está incluida **solo para poder mostrar cómo falla**. No es un indicador válido de bondad de ajuste fuera de la regresión lineal ordinaria (Spiess y Neumeyer, *BMC Pharmacology* 10:6, 2010).

Para usarlo fuera de Colab, alcanza con poner el archivo en la misma carpeta que el cuaderno y hacer `import lab1_utils as lab`.

---

## Bibliografía de referencia

- J. R. Taylor, *An Introduction to Error Analysis*, 2.ª ed., University Science Books, 1997. Secc. 2.5 (cifras significativas), cap. 6 (Chauvenet), cap. 8 (cuadrados mínimos).
- P. R. Bevington y D. K. Robinson, *Data Reduction and Error Analysis for the Physical Sciences*, 3.ª ed., McGraw-Hill, 2003. Cap. 11 (bondad de ajuste).
- JCGM 100:2008, *Guide to the Expression of Uncertainty in Measurement* (GUM). Secc. 7.2.6 (reporte de la incerteza).
- A.-N. Spiess y N. Neumeyer, "An evaluation of R² as an inadequate measure for nonlinear models…", *BMC Pharmacology* **10**:6, 2010.
- F. J. Anscombe, "Graphs in Statistical Analysis", *The American Statistician* **27**(1), 17, 1973.
- D. W. Scott, "On optimal and data-based histograms", *Biometrika* **66**(3), 605, 1979.
- L. Lyons, *Statistics for Nuclear and Particle Physicists*, Cambridge University Press, 1986.
