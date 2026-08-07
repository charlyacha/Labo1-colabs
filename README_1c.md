# Colabs de Laboratorio 1 — Departamento de Física, FCEN-UBA

Diez notebooks que acompañan el cronograma de 16 clases del 2.º cuatrimestre 2026, más un módulo
de funciones reutilizables.

## Contenido

| Archivo | Clase | Qué introduce |
|---|---|---|
| `01_Primeros_pasos_medir_y_reportar.ipynb` | 1 | Entorno, nonio y micrómetro, error de cero, cuantización y $\Delta/\sqrt{12}$, arrays, primer gráfico |
| `02_Mediciones_indirectas_y_propagacion.ipynb` | 2 | SymPy, tabla de contribuciones, redondeo y `reportar()`, sistemático instrumental |
| `03_Estadistica_gaussiana_y_compatibilidad.ipynb` | 3 | Promedio como estimador, `ddof=1`, $s$ vs. SEM, Scott, $1/\sqrt{N}$ vs. $1/N$, gaussiana predicha, TCL, compatibilidad, Chauvenet |
| `04_Buscando_la_ley_pendulo_photogate_y_ajuste.ipynb` | 4 | Extracción de flancos, $\sigma_t$ desde $f_s$, cuadrados mínimos, exponente por log-log, residuos, primera determinación de $g$ |
| `05_Resorte_R_residuos_y_rango_de_validez.ipynb` | 5 | Ordenada al origen con criterio, tensión inicial de bobinado, límites de $R^2$, Anscombe, rango de validez, histéresis |
| `06_Ajuste_ponderado_covarianza_y_chi2.ipynb` | 6 | `absolute_sigma=True`, `pcov`, $\chi^2_\nu$ y p-valor, promedio ponderado como ajuste a constante |
| `07_Datos_reales_adquisicion_y_derivadas.ipynb` | 7 | Tracker/Pasco/photogate, máscaras, ruido en la derivada numérica |
| `08_Determinacion_de_g_y_comparacion_de_modelos.ipynb` | 8 | Ajuste cuadrático ponderado, modelos anidados, sistemático en el resultado |
| `09_Senales_periodicas_y_determinacion_del_periodo.ipynb` | 9 | `find_peaks`, por qué promediar intervalos usa solo dos picos |
| `10_Ajuste_no_lineal_y_por_que_R2_miente.ipynb` | 10 | `curve_fit` no lineal, semillas físicas, correlación de parámetros, fallas de $R^2$ |
| `lab1_utils.py` | 14 | Todas las funciones del curso, listas para reutilizar |

Cada notebook declara sus requisitos previos y no usa ninguna construcción de Python que no haya
sido presentada antes. Todos son autocontenidos: generan sus propios datos de ejemplo con semilla
fija, en celdas marcadas `# --- DATOS DE EJEMPLO (reemplazar por los propios) ---`.


---

## Revisión del 5 de agosto de 2026: reordenamiento de los Bloques I y II

El cronograma se reordenó y con él la numeración de los cuadernos. **Las fechas y la cantidad de
clases no cambiaron**: sólo se movió el contenido, y sólo entre las Clases 1 y 6.

### Correspondencia con la versión anterior

| Nuevo | Reemplaza a | Qué cambió |
|---|---|---|
| `01_Primeros_pasos_medir_y_reportar` | mismo nombre | **Salió** el redondeo y `reportar()` (se fueron al 02). **Entró** el nonio del calibre, el tambor y el trinquete del micrómetro, la medición y corrección del error de cero, y la simulación de cuantización que muestra que el error de resolución es uniforme con $\sigma = \Delta/\sqrt{12}$ |
| `02_Mediciones_indirectas_y_propagacion` | `04_Propagacion_de_incertezas` | La densidad se escribe en función del **diámetro** y no del radio. Datos de ejemplo corregidos: los anteriores daban $\rho \approx 16\,600$ kg/m³, que no corresponde a ningún material. **Entró** el redondeo y `reportar()`. **Salió** la simulación de "promediar no corrige el sistemático" (se fue al 03) |
| `03_Estadistica_gaussiana_y_compatibilidad` | `02_Estadistica_de_una_variable` + `03_Gaussiana_TCL_y_compatibilidad` | Fusión de los dos. **Entró** el faro luminoso, el contraste $1/\sqrt{N}$ vs. $1/N$, la incerteza de $s$ misma, el criterio de Chauvenet, la simulación del sistemático, y la actividad opcional del contraste de formas. **Salió** el promedio ponderado (se fue al 06) |
| `04_Buscando_la_ley_pendulo_photogate_y_ajuste` | *(nuevo)* | Toma la maquinaria de ajuste del ex-05 y la idea de extracción de período del ex-09 |
| `05_Resorte_R_residuos_y_rango_de_validez` | `05_Cuadrados_minimos_y_el_coeficiente_R` | Reformulado: la medición ya **no sale del rango elástico**. La falla del modelo viene de la tensión inicial de bobinado, del lado de las masas chicas, y es no destructiva |
| `06_Ajuste_ponderado_covarianza_y_chi2` | mismo nombre | Datos de ejemplo cambiados por los del resorte de la Clase 5 con repeticiones. **Entró** la sección del promedio ponderado como ajuste a modelo constante, la verificación de que el cociente de errores es $\sqrt{\chi^2_\nu}$, y la advertencia sobre comparar dos análisis del mismo conjunto |
| `07`, `08`, `10` | mismo nombre | Sin cambios de contenido |
| `09_Senales_periodicas_y_determinacion_del_periodo` | mismo nombre | **Entró** la sección 1bis que conecta con el Colab 04: mismo problema (marca temporal reproducible sobre una señal periódica) con dos soluciones según la forma de la señal |

La variante con celdas incompletas para completar en clase pasa de llamarse `02b` a **`03b`**.

### Corrección de formato en los diez cuadernos

Los diez `.ipynb` de la versión anterior guardaban el campo `source` de cada celda como una lista de
líneas **sin el salto de línea final**. El formato `nbformat` exige que cada elemento de esa lista
termine en `\n` (salvo el último), porque el intérprete concatena la lista con `''`. Sin los saltos,
una celda como

```python
import numpy as np
import matplotlib.pyplot as plt
```

se convierte al ejecutarse en `import numpy as npimport matplotlib.pyplot as plt`, y **la primera
celda de código de cada cuaderno falla con `SyntaxError`**.

Los diez archivos de esta entrega están normalizados y verificados: se ejecutaron de punta a punta
con `nbclient` y los diez terminan sin errores. Si venís usando los archivos anteriores, conviene
reemplazarlos aunque no te interese el reordenamiento.

### Nuevas funciones en `lab1_utils.py`

| Función | Origen |
|---|---|
| `reportar(x, dx, unidad)` | alias de `formatear()`, que es el nombre usado en los Colabs 02 en adelante |
| `sigma_resolucion(delta)` | $\Delta/\sqrt{12}$ (Colabs 01 y 03) |
| `incerteza_de_s(m)` | $1/\sqrt{2(m-1)}$: la incerteza relativa de una desviación estándar estimada con $m$ datos (Colab 03) |
| `chauvenet(x)` | máscara de los datos que sobreviven al criterio (Colab 03) |
| `tiempos_de_flanco(t, V, umbral, tipo)` | extracción de eventos de una señal de dos niveles (Colab 04) |
| `periodo_por_ajuste(t_eventos, eventos_por_periodo)` | período por ajuste de $t_n$ vs. $n$, con el parámetro que resuelve la trampa de las dos interrupciones por período (Colabs 04 y 09) |
| `propagar(expr, ...)` | propagación con tabla de contribuciones sobre expresiones de SymPy (Colab 02) |

Además, `redondear_con_error()` pasa a conservar **dos cifras significativas en la incerteza cuando
la primera es 1 o 2** (antes, sólo cuando era 1), que es la regla que usan todos los cuadernos nuevos.

### Datos de ejemplo que conviene reemplazar por los reales

Los cuadernos generan sus propios datos con semilla fija, pero tres conjuntos dependen de parámetros
del laboratorio que hay que confirmar:

- **Colab 03 — el faro.** Se supone $T_0 = 1{,}437$ s y $\sigma_{\text{ev}} = 0{,}17$ s para el
  cronometrado manual. Si el faro tiene período seteable, conviene fijarlo y anotarlo.
- **Colab 04 — el photogate.** `f_s = 1000.0` Hz. Todo el balance de incertezas de la sección 2.4
  depende de ese número. Si el sistema entrega directamente tiempos de flanco en lugar de la señal
  muestreada, la sección 2 hay que adaptarla.
- **Colab 05 — el resorte.** Se supone una tensión inicial de bobinado $F_0 = 0{,}15$ N. Si los
  resortes de la cátedra no la tienen, las secciones 2 a 5 se quedan sin fenómeno. Se verifica en dos
  minutos: colgar 5, 10 y 20 g y mirar si elonga desde el primer gramo o si hay un umbral.

---

## Cómo publicarlos y compartirlos

Ésta es la parte que conviene decidir bien de entrada, porque cambiar de esquema a mitad de
cuatrimestre significa reemplazar todos los enlaces de la página de la materia.

### La recomendación: repositorio público de GitHub

**Es el esquema que usaría.** El notebook vive en GitHub, no en tu Drive, y Colab lo abre desde ahí:

```
https://colab.research.google.com/github/USUARIO/labo1-colabs/blob/main/03_Estadistica_gaussiana_y_compatibilidad.ipynb
```

Ese enlace abre el notebook en Colab **en modo playground**: el estudiante puede ejecutarlo y
modificarlo, pero los cambios no van a ningún lado salvo que él haga `Guardar una copia en Drive`,
que la guarda en *su* Drive.

Por qué es la mejor opción:

- **Desacople total de tu cuenta.** Tu Drive no interviene. No hay permisos que revisar, ni riesgo
  de que un enlace mal configurado exponga otra cosa, ni tu nombre de cuenta asociado al archivo.
  Ésta es la respuesta directa a tu pregunta sobre cómo compartir sin dar acceso a nada más:
  el problema desaparece porque no hay nada tuyo del otro lado del enlace.
- **Enlaces estables.** La URL no cambia aunque edites el notebook. Podés corregir un error un
  martes a la noche y el enlace de la página de la materia ya apunta a la versión corregida.
- **Control de versiones.** Queda registro de qué cambiaste y cuándo. Si algo se rompe, volvés atrás.
  Con archivos sueltos en Drive esto no existe.
- **Los estudiantes pueden leerlos sin ejecutarlos.** GitHub renderiza `.ipynb` directamente en el
  navegador, así que quien solo quiera consultar un tema no necesita abrir Colab.
- **Reutilizable entre cuatrimestres y entre docentes.** Otra cátedra puede hacer un *fork* y
  adaptarlo, con atribución automática.

**Montaje, una sola vez (unos 20 minutos):**

1. Crear cuenta en github.com si no tenés, y un repositorio **público** llamado, por ejemplo,
   `labo1-colabs`.
2. Subir los `.ipynb`, el `lab1_utils.py` y este README (se puede desde el navegador: *Add file →
   Upload files*, arrastrando todo junto).
3. En cada notebook, reemplazar `TU-USUARIO/labo1-colabs` del badge de la primera celda por tu
   usuario y repositorio reales. (Un `sed` sobre los diez archivos: 
   `sed -i 's|TU-USUARIO/labo1-colabs|tuusuario/labo1-colabs|g' *.ipynb`)
4. En la página de WordPress de la materia, enlazar las URLs de `colab.research.google.com/github/...`.

**Para editar después:** abrís el notebook desde el enlace de Colab y usás
`Archivo → Guardar una copia en GitHub`, eligiendo el mismo repositorio y archivo. Colab pide
autorización a GitHub la primera vez y después es un clic.

**Archivos de datos** (los `.txt` de ejemplo, los exports de Tracker): también van al repositorio, y
los notebooks los leen con

```python
!wget -q https://raw.githubusercontent.com/USUARIO/labo1-colabs/main/datos/histograma_ejemplo.txt
```

### La alternativa: Google Drive

Funciona, y es lo que hacen la mayoría de las cátedras relevadas, pero tiene fricciones reales.

Si vas por este camino:

- **Compartí archivos, no carpetas.** Compartir una carpeta comparte todo lo que haya adentro
  **y todo lo que pongas adentro en el futuro**. Es la manera más común de exponer algo sin querer.
- **Siempre "Cualquier persona con el enlace → Lector"**, nunca *Editor*. Con permiso de editor, un
  estudiante que escriba en el notebook original lo rompe para todos los demás, sin mala intención.
- Compartir un archivo de Drive **no** da acceso a ningún otro archivo tuyo — eso es cierto y es una
  preocupación común pero infundada. Lo que sí queda expuesto es el nombre de tu cuenta como
  propietario y, para quien tenga el enlace, el archivo completo con su historial de revisiones.
- Si te incomoda mezclar la cuenta institucional con el material del curso, la opción prolija es una
  **cuenta de Google dedicada a la materia** (`labo1.acha@gmail.com` o similar), que además se puede
  transferir a quien tome la cátedra el año que viene.

Lo que hace incómodo a Drive para este uso no es la seguridad sino el mantenimiento: los enlaces son
opacos (`drive.google.com/file/d/1a2B3c...`), no hay historial legible, y actualizar un notebook
significa reemplazar el archivo o pisar el original sin registro de qué cambió.

### Lo que no recomiendo

- **Repositorio privado de GitHub.** El cargador de Colab desde GitHub no funciona para quien no
  tenga acceso, así que los estudiantes no podrían abrirlos. Para material docente, público es lo
  correcto.
- **Compartir la carpeta de Drive con permiso de edición** para que los estudiantes "entreguen ahí".
  Se desordena en dos semanas. Las entregas conviene canalizarlas por el campus o por correo.

### Qué hace el estudiante

Una sola instrucción, que conviene repetir en la Clase 1 y poner en la página:

> Abrí el enlace, y antes de escribir nada hacé `Archivo → Guardar una copia en Drive`. Vas a
> trabajar sobre tu copia. El original queda intacto para el resto del curso.

Cada notebook lo dice en su primera celda.

---

## Notas para el docente

- **Ejecutar antes de la clase.** Los diez notebooks fueron probados de punta a punta, pero Colab
  actualiza versiones de librerías sin aviso. Correr el notebook de la semana el día anterior evita
  sorpresas.
- **La Clase 2 conviene entregarla con celdas incompletas.** Es la de mayor densidad conceptual y
  cae cuando los estudiantes recién conocen el entorno. Si pelean con la sintaxis, no les queda
  atención para el concepto de estimador. Una versión "con huecos" del notebook 02, con los
  `TODO` marcados, es probablemente la mejor inversión de tiempo de preparación de todo el curso.
- **El notebook 10 es el cuello de botella.** El ajuste no lineal falla en convergencia si las
  semillas son malas. Conviene tener a mano un conjunto de datos de respaldo con solución conocida
  por si el experimento del día no da señal limpia.
- **Auditoría pendiente sobre material reutilizado.** Cualquier Colab de otras cátedras que se
  adapte debe revisarse buscando `np.std(` sin `ddof=1`. Es un error silencioso: no falla, solo
  devuelve un número un poco chico.

---

## Bibliografía citada en los notebooks

- Anscombe, F. J. "Graphs in Statistical Analysis", *The American Statistician* **27**(1), 17–21 (1973).
- Scott, D. W. "On optimal and data-based histograms", *Biometrika* **66**(3), 605–610 (1979).
- Spiess, A.-N. & Neumeyer, N. "An evaluation of R² as an inadequate measure for nonlinear models in
  pharmacological and biochemical research", *BMC Pharmacology* **10**:6 (2010).
- Cushing, J. T. "The spring-mass system revisited", *Am. J. Phys.* **52**, 925 (1984).
- Taylor, J. R. *An Introduction to Error Analysis*, 2.ª ed., University Science Books (1997).
- Bevington, P. R. & Robinson, D. K. *Data Reduction and Error Analysis for the Physical Sciences*,
  3.ª ed., McGraw-Hill (2002).
- Lyons, L. *A Practical Guide to Data Analysis for Physical Science Students*, Cambridge U. P. (1991).
- Kong, Q., Siauw, T. & Bayen, A. *Python Programming and Numerical Methods*, Academic Press (2020),
  de acceso abierto.
