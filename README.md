# Laboratorio 1 — Cuadernos de Colab

**Departamento de Física, FCEN-UBA · Cátedra Acha · 2.º cuatrimestre 2026**

Repositorio de los cuadernos de análisis de datos del curso. Cada clase experimental tiene
un cuaderno con el número correspondiente; se abre directamente en Google Colab desde el
enlace de la tabla. Todos son autocontenidos: generan sus propios datos de ejemplo con
semilla fija, en celdas marcadas `# DATOS DE EJEMPLO — REEMPLAZAR POR LOS PROPIOS`, y no
usan ninguna construcción de Python que no haya sido presentada en un cuaderno anterior.

Los archivos están verificados: todos los cuadernos ejecutables se corren de punta a punta
con `nbclient` y terminan sin errores. El único que no corre completo es `03b`, que es
incompleto **a propósito** (ver la tabla).

---

## Contenido

| Archivo | Clase | Pregunta / contenido | Abrir |
|---|---|---|---|
| `01_Primeros_pasos_medir_y_reportar.ipynb` | 1 | ¿Por qué no medimos todos lo mismo? Entorno, nonio y micrómetro, error de cero, resolución y Δ/√12, arrays, primer gráfico | [Colab](https://colab.research.google.com/github/charlyacha/Labo1-colabs/blob/main/01_Primeros_pasos_medir_y_reportar.ipynb) |
| `02_Mediciones_indirectas_y_propagacion.ipynb` | 2 | ¿De qué material está hecho este cuerpo? SymPy, tabla de contribuciones, redondeo y `reportar()`, sistemático instrumental | [Colab](https://colab.research.google.com/github/charlyacha/Labo1-colabs/blob/main/02_Mediciones_indirectas_y_propagacion.ipynb) |
| `03_Estadistica_gaussiana_y_compatibilidad.ipynb` | 3 | ¿Cuántas veces vale la pena medir, y qué número reporto? Estimador, `ddof=1`, *s* vs. SEM, 1/√N vs. 1/N, gaussiana predicha, TCL, compatibilidad, Chauvenet | [Colab](https://colab.research.google.com/github/charlyacha/Labo1-colabs/blob/main/03_Estadistica_gaussiana_y_compatibilidad.ipynb) |
| `03b_..._INCOMPLETO.ipynb` | 3 | **Variante con celdas incompletas para completar en clase.** No corre hasta llenar los huecos: es a propósito | [Colab](https://colab.research.google.com/github/charlyacha/Labo1-colabs/blob/main/03b_Estadistica_gaussiana_y_compatibilidad_INCOMPLETO.ipynb) |
| `04_Buscando_la_ley_pendulo_photogate_y_ajuste.ipynb` | 4 | ¿Puedo descubrir de qué depende el período sin suponerlo? Extracción de flancos, σ_t desde f_s, cuadrados mínimos, exponente por log-log, residuos, primera *g* | [Colab](https://colab.research.google.com/github/charlyacha/Labo1-colabs/blob/main/04_Buscando_la_ley_pendulo_photogate_y_ajuste.ipynb) |
| `05_Resorte_R_residuos_y_rango_de_validez.ipynb` | 5 | ¿Hasta dónde vale la ley de Hooke? Ordenada al origen con criterio, tensión inicial de bobinado, límites de R², Anscombe, rango de validez, histéresis | [Colab](https://colab.research.google.com/github/charlyacha/Labo1-colabs/blob/main/05_Resorte_R_residuos_y_rango_de_validez.ipynb) |
| `06_Ajuste_ponderado_covarianza_y_chi2.ipynb` | 6 | ¿Cambia *k* si tomo en serio las barras de error? `absolute_sigma=True`, `pcov`, χ²ᵥ y p-valor, promedio ponderado como ajuste a constante | [Colab](https://colab.research.google.com/github/charlyacha/Labo1-colabs/blob/main/06_Ajuste_ponderado_covarianza_y_chi2.ipynb) |
| `07_Datos_reales_adquisicion_y_derivadas.ipynb` | 7 | ¿Tres instrumentos distintos miden la misma velocidad? Tracker/Pasco/photogate, archivos sucios, máscaras, ruido en la derivada | [Colab](https://colab.research.google.com/github/charlyacha/Labo1-colabs/blob/main/07_Datos_reales_adquisicion_y_derivadas.ipynb) |
| `08_Determinacion_de_g_y_comparacion_de_modelos.ipynb` | 8 | ¿Puedo determinar *g* en este laboratorio? Ajuste cuadrático ponderado, modelos anidados, sistemático en el resultado | [Colab](https://colab.research.google.com/github/charlyacha/Labo1-colabs/blob/main/08_Determinacion_de_g_y_comparacion_de_modelos.ipynb) |
| `09_Senales_periodicas_y_determinacion_del_periodo.ipynb` | 9 | ¿Puedo medir una masa sin usar una balanza? `find_peaks`, por qué promediar intervalos usa solo dos picos, período por ajuste | [Colab](https://colab.research.google.com/github/charlyacha/Labo1-colabs/blob/main/09_Senales_periodicas_y_determinacion_del_periodo.ipynb) |
| `10_Ajuste_no_lineal_y_por_que_R2_miente.ipynb` | 10 | ¿Cuánto tarda en morirse esta oscilación, y cómo sé que mi ajuste es correcto? `curve_fit` no lineal, semillas físicas, correlación de parámetros, fallas de R² | [Colab](https://colab.research.google.com/github/charlyacha/Labo1-colabs/blob/main/10_Ajuste_no_lineal_y_por_que_R2_miente.ipynb) |
| `11_estudiante_Viscosidad_por_caida_de_esferas.ipynb` | 11 | ¿Puedo medir la viscosidad dejando caer una esfera? **Cuaderno casi vacío a propósito** (ensayo general de la Práctica Especial): trae la física, las consignas y la autoevaluación; el código lo escribe el estudiante | [Colab](https://colab.research.google.com/github/charlyacha/Labo1-colabs/blob/main/11_estudiante_Viscosidad_por_caida_de_esferas.ipynb) |
| `11_docente_Viscosidad_por_caida_de_esferas.ipynb` | 11 | **Versión resuelta (no se reparte).** El análisis completo del que el anterior es el andamiaje: Stokes, corrección de Ladenburg, ajuste ponderado, χ²ᵥ, Reynolds. Para consultar la solución en el aula | [Colab](https://colab.research.google.com/github/charlyacha/Labo1-colabs/blob/main/11_docente_Viscosidad_por_caida_de_esferas.ipynb) |
| `12_Choques_y_conservacion.ipynb` | 12 | ¿Se conserva realmente el momento en un choque? Δp con incerteza, energía remanente, restitución por ajuste exponencial | [Colab](https://colab.research.google.com/github/charlyacha/Labo1-colabs/blob/main/12_Choques_y_conservacion.ipynb) |
| `S1_El_experimento_repetido_y_la_distribucion_de_chi2.ipynb` | — | **Optativo.** De dónde sale la distribución de χ² que después se usa como test. Posterior a la Clase 6 | [Colab](https://colab.research.google.com/github/charlyacha/Labo1-colabs/blob/main/S1_El_experimento_repetido_y_la_distribucion_de_chi2.ipynb) |
| `PE_Plantilla_analisis_Practica_Especial.ipynb` | 14–16 | **Plantilla** de análisis para la Práctica Especial: la estructura, para llenar con los datos propios | [Colab](https://colab.research.google.com/github/charlyacha/Labo1-colabs/blob/main/PE_Plantilla_analisis_Practica_Especial.ipynb) |
| `lab1_utils.py` | — | Módulo con todas las funciones del curso, listas para reutilizar | — |

---

## Cómo lo usa el estudiante

1. Abrí el cuaderno con el enlace **Colab** de la tabla.
2. **Archivo → Guardar una copia en Drive.** Esto es lo primero, antes de escribir nada:
   lo que hagas sin copiar primero no se guarda en ningún lado, y el original queda intacto
   para el resto del curso.
3. Ejecutá la primera celda de código. Baja `lab1_utils.py` automáticamente desde este
   repositorio.
4. Recorré el cuaderno **de arriba hacia abajo.** Las celdas no son independientes: si algo
   tira `NameError`, casi siempre es porque salteaste una.
5. Para subir tus datos, usá el ícono de carpeta del panel izquierdo.
6. Antes de entregar: **Entorno de ejecución → Reiniciar y ejecutar todo.** Ése, y no otro,
   es el estado en el que se entrega un cuaderno.

---

## El módulo `lab1_utils.py`

No es una caja negra. Cada función se escribe antes, a mano, en el cuaderno de la clase
correspondiente; el módulo existe para no tener que reescribirlas en la Práctica Especial.
En la Clase 15 se arma con los estudiantes juntando lo que cada uno programó a lo largo del
cuatrimestre.

**Presentación:** `estilo_lab1()`, `guardar_figura()`
**Reporte:** `redondear_con_error()`, `formatear()`, `reportar()`
**Una variable:** `estadisticos()`, `bins_scott()`, `bins_freedman_diaconis()`, `sigma_resolucion()`, `incerteza_de_s()`, `chauvenet()`
**Comparación:** `compatibilidad()`, `promedio_ponderado()`
**Bondad de ajuste:** `chi2_reducido()`, `matriz_correlacion()`, `R2()`
**Ajuste:** `ajustar()`, `grafico_con_residuos()`
**Señales y propagación:** `tiempos_de_flanco()`, `periodo_por_ajuste()`, `propagar()`
**Lectura de archivos:** `leer_datos()`, `leer_tracker()`

Tres decisiones de diseño que conviene conocer, porque contradicen los valores por defecto
de las librerías y los tres fallan **en silencio** —no dan error, devuelven un número
plausible y equivocado:

- **`estadisticos()` usa `ddof=1`.** El default de `np.std` es `ddof=0`, que calcula la
  desviación estándar poblacional y subestima la dispersión de una muestra. Con tres datos,
  la subestimación es del 18 %.
- **`ajustar()` usa `absolute_sigma=True`** siempre que se le pasen barras de error. El
  default de `curve_fit` es `False`, que reescala la covarianza por el χ²ᵥ del propio ajuste
  y **destruye el diagnóstico de bondad de ajuste por construcción.** La trampa gemela está
  en `np.polyfit(..., cov=True)`, que necesita `cov='unscaled'`.
- **`R2()` está incluida solo para poder mostrar cómo falla.** No es un indicador válido de
  bondad de ajuste fuera de la regresión lineal ordinaria (Spiess y Neumeyer, *BMC
  Pharmacology* 10:6, 2010). El indicador correcto es χ²ᵥ con p-valor.

Los lectores `leer_datos()` y `leer_tracker()` toleran las patologías de los archivos
reales: encabezados de metadatos de largo variable, separador variable (tabulador, `;`, `,`
o espacios), coma decimal, y `NaN` en los cuadros no trackeados. `leer_tracker()` conserva
los `NaN` a propósito: recortar el tramo útil con máscaras booleanas es parte del ejercicio
del Colab 07.

Para usarlo fuera de Colab, alcanza con poner el archivo en la misma carpeta que el cuaderno
y hacer `import lab1_utils as lab`.

---

## Notas para el docente

- **Ejecutá el cuaderno de la semana el día anterior.** Están todos probados de punta a
  punta, pero Colab actualiza versiones de librerías sin aviso.
- **El Colab 03 conviene darlo con celdas incompletas** (`03b`) si el grupo viene peleando
  con la sintaxis. Es la clase de mayor densidad conceptual y cae cuando recién conocen el
  entorno; si gastan la atención en `matplotlib`, no les queda para el concepto de estimador.
- **El Colab 10 es el cuello de botella técnico.** El ajuste no lineal puede no converger si
  las semillas son malas. Tené a mano un conjunto de datos de respaldo con solución conocida.
- **El Colab 11 es el ensayo general de la Práctica Especial.** Está pensado para que el
  grupo arme el análisis con la mínima intervención: se responden preguntas conceptuales,
  no se escribe el código por ellos. Es el mejor predictor de qué grupos van a poder solos
  con la Práctica Especial.
- **Auditá cualquier Colab que adaptes de otra fuente** buscando `np.std(` sin `ddof=1` y
  `curve_fit` sin `absolute_sigma=True`. Son errores silenciosos.

---

## Parámetros del laboratorio a confirmar

Los cuadernos generan datos de ejemplo con semilla fija, pero tres dependen de parámetros
reales que conviene verificar antes de empezar. Los tres se comprueban en minutos y los tres
invalidan una sección entera si no se cumplen.

| Cuaderno | Supuesto | Cómo se verifica |
|---|---|---|
| `03` | Faro con T₀ ≈ 1,44 s y σ ≈ 0,17 s en el cronometrado manual | Cronometrar veinte destellos. Si el período es seteable, conviene fijarlo y anotarlo |
| `04` | Photogate a f_s = 1000 Hz, entregando V(t) muestreado | Mirar la configuración de adquisición. Si entrega tiempos de flanco directos, la sección 2 se adapta |
| `05` | Resorte con tensión inicial de bobinado F₀ ≈ 0,15 N | Colgar 5, 10 y 20 g y ver si elonga desde el primer gramo o si hay un umbral |

---

## Publicación y mantenimiento

El repositorio es `https://github.com/charlyacha/Labo1-colabs`, rama `main`. Los enlaces de
Colab de la tabla abren cada cuaderno **en modo playground**: el estudiante lo ejecuta y lo
modifica, pero los cambios no van a ningún lado salvo que haga `Guardar una copia en Drive`.

Es GitHub y no un Drive compartido por tres razones: el enlace genera una copia limpia para
cada estudiante sin que puedan tocar el original; hay historial de versiones de cada cambio;
y se puede corregir un error a mitad de cuatrimestre sin volver a repartir enlaces.

Los archivos de datos van al mismo repositorio y los cuadernos los leen con

```python
!wget -q https://raw.githubusercontent.com/charlyacha/Labo1-colabs/main/datos/histograma_ejemplo.txt
```

Si cambia el usuario o el nombre del repositorio, hay que actualizar en cada cuaderno la URL
del enlace y la del `wget`, con una línea:

```bash
sed -i 's|charlyacha/Labo1-colabs|NUEVO-USUARIO/NUEVO-REPO|g' *.ipynb
```

---

## Bibliografía de referencia

- J. R. Taylor, *An Introduction to Error Analysis*, 2.ª ed., University Science Books, 1997.
  §2.5 (cifras significativas), §4.2 y §4.4 (*s* vs. SEM), cap. 6 (Chauvenet), cap. 8
  (cuadrados mínimos), cap. 12 (χ²).
- P. R. Bevington y D. K. Robinson, *Data Reduction and Error Analysis for the Physical
  Sciences*, 3.ª ed., McGraw-Hill, 2003. Cap. 11 (bondad de ajuste).
- JCGM 100:2008, *Guide to the Expression of Uncertainty in Measurement* (GUM). §7.2.6
  (reporte de la incerteza), §F.2.2.1 (resolución uniforme).
- A.-N. Spiess y N. Neumeyer, "An evaluation of R² as an inadequate measure for nonlinear
  models…", *BMC Pharmacology* **10**:6, 2010.
- F. J. Anscombe, "Graphs in Statistical Analysis", *The American Statistician* **27**(1),
  17–21, 1973.
- D. W. Scott, "On optimal and data-based histograms", *Biometrika* **66**(3), 605–610, 1979.
- L. Lyons, *A Practical Guide to Data Analysis for Physical Science Students*, Cambridge
  University Press, 1991.
- Q. Kong, T. Siauw y A. Bayen, *Python Programming and Numerical Methods*, Academic Press,
  2020 (acceso abierto).
