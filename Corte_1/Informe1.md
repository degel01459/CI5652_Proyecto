---
lang: es-ES
header-includes:
   - \renewcommand{\tablename}{Tabla}
---

# Informe 1

Elaborado por:

Alejandro Zambrano (17-10684), Ángel Rodríguez (15-11669),
Francisco Márquez (12-11163), Kevin Briceño (15-11661), Sergio Carrillo
(14-11315).

Fecha: 01/02/2026

## Introducción

El problema de satisficibilidad booleana (SAT) es el problema de saber si una
fórmula en forma normal conjuntiva tiene una asignación de variables que la haga
verdadera. En casos más complejos o directamente en casos de contradicciones, no
se puede satisfacer la fórmula pero se busca minimizar las cláusulas que están
insatisfechas: este es el problema MaxSAT. Tal problema cuenta con varias
especializaciones sobre las cláusulas (que pueden tener ponderación o algún
tipo de condición de obligatoriedad). En este estudio se trabajó con el problema
MaxSAT puro.

Para el estudio de este problema se usó el _benchmark Causal Discovery_ el cual
es una transcripción de problemas de estadística para encontrar relaciones de
causalidad entre sucesos correlacionados. Este _benchmark_ explota al algoritmo
de un solucionador por dar lugar a cláusulas densas (con muchas variables) por
lo que es impráctico plantear un solucionador exacto en este estudio y se usará
uno disponible en la red para comparar la solución heurística y metaheurísitca
de búsqueda local propuesta.

## Descripción del problema

El problema de satisficibilidad booleana (SAT) es el problema de saber si dada
una fórmula booleana en forma normal conjuntiva (conjunción de disyunciones),
existe alguna asignación de variables tal que la fórmula sea verdadera (se
satisfaga). El problema MaxSAT es un problema de optimización donde dada una
fórmula SAT, que muy probablemente es insatisfacible (una contradicción) o donde
muy pocas asignaciones de variables la satisfacen, se quiere una asignación de
variables tal que se minimice la cantidad de cláusulas insatisfechas o,
equivalentemente, se maximice las cláusulas satisfechas. \[1\]

En la literatura, existen varias especializaciones de MaxSAT \[1\]:

- MaxSAT ponderado: existen cláusulas que tienen mayor preferencia para ser
  satisfechas (tienen mayor peso o ponderación).
- MaxSAT parcial: existen cláusulas que obligatoriamente deben satisfacerse
  (_hard_) y otras que se permite que no se satisfagan (_soft_).
- MaxSAT ponderado parcial: Una combinación de los casos anteriores.

En este estudio se trabajará con el problema simple MAXSAT donde puede decirse
que todas las cláusulas son _soft_ y con la misma ponderación con el enfoque de
minimización de cláusulas insatisfechas.

## _Benchmark_: Descubrimiento causal

El Descubrimiento Causal (_causal discovery_) es el área de la inteligencia
artificial y la estadística cuyo objetivo es inferir las relaciones de causa y
efecto a partir de datos, en lugar de simplemente encontrar correlaciones. \[2\]
Esto permite definir relaciones lógicas entre sucesos y puede traducirse al
lenguaje de la lógica proposicional (dominio de SAT) con el concepto de la
D-separación. \[3\] Cuando se trabaja con datos de muestras reales, las pruebas
estadísticas de independencia suelen producir resultados contradictorios debido
a la variabilidad estadística. Los métodos previos con frecuencia fallaban o se
volvían ineficientes ante estas inconsistencias. \[4\] Mediante optimización de
restricciones el problema de descubrimiento causal puede ser llevado al problema
de MaxSAT. \[4\] Se propuso el uso de MaxSAT en \[5\] como el motor de búsqueda
para el descubrimiento causal basado en restricciones.

Cada archivo del _benchmark_ plantea una cantidad de casos para los cuales se
aplica el descubrimiento causal y eso genera una cantidad de variables y
cláusulas que llevan a los resolvedores de MaxSAT a su límite. El resumen de los
archivos se presenta en la tabla 1.

Tabla 1. Descripción de los archivos del _benchmark_

| Archivo | Variables | Cláusulas |     Casos      | Lógica |
| :-----: | :-------: | :-------: | :------------: | :----: |
|   n5    |   61600   |  221790   | 500 - 1k - 10k | UAI13  |
|   n6    |  328107   |  1206162  |   500 - 10k    | UAI13  |
|   n6    |   12764   |   46236   | 500 - 1k - 10k | UAI14  |
|   n7    |   40290   |  145910   | 500 - 1k - 10k | UAI14  |

## Solución exacta: EvalMaxSAT

EvalMaxSAT \[6\] es un solucionador de MaxSAT implementado en C++ y basado en el
solucionador `CaDiCaL` \[7\] y el algortimo aprendizaje un nivel a la vez (OLL
por sus siglas en inglés: _One-Level-at-a-time Learning_) \[8\].

## Solución heurística

Se propuso una heurística de ramificación estática (_Static Branching
Heuristics_) basada en el conteo de literales. \[9,10\] La idea de la heurística
se centra en tomar la variable cuya frecuencia de aparición sea máxima (regla de
MOM \[11\]) para asignar el valor de verdad que maximice el número de cláusulas
que se satrisfacen en cada momento. La heurísitca es una simplificación de la
heurística de Jerislow-Wang \[12\] que trabaja con esta misma idea pero con una
función de costo afectada por el peso de cada cláusula.

## Búsqueda local: vecindad 1-_flip_

La vecindad propuesta para cada solución consiste en alternar el valor de verdad
asociado a una variable en una cláusula insatisfecha lo que implica que cada
vecino de la solución difiere en exactamente una variable.

## Búsqueda local iterada: perturbación y criterio de aceptación

La perturbación propuesta para esta metaheurística es la generalización _k-flip_
sobre variables en cláusulas satisfechas. Esto incrementa la probabilidad de
escape del cuenco de atracción que define la solución encontrada por la solución
heurística propuesta al cambiar drásticamente la forma de la solución inicial.
La búsqueda local se sigue haciendo con 1-_flip_.

## Resultados de corridas

En la tabla 2 se muestran los resultados en tiempo del promedio de treinta
ejecuciones de cada archivo del _benchmark_. El número entre paréntesis
corresponde a la desviación estándar sobre la última cifra significativa.

Tabla 2. Resultados de las corridas

| Archivo | Casos | S. exacto (s) | S. heurístico (s) | B. L. (s) | B. L. I. (s)
|:-------:|:-----:|:-------------:|:-----------------:|:---------:|:-----------:
|  n5 i2  |  500  |    0.25(5)    |     110(3)        |   0.0003  |  0.0095
|  n5 i4  |  500  |    0.24(2)    |     1.2(1)×10²    |   0.0003  |  0.0094
|  n5 i5  |  10k  |    0.25(2)    |     113(8)        |   0.0003  |  0.0091
|  n5 i7  |  1k   |    0.116(8)   |     108(2)        |   0.0003  |  0.0117
|  n5 i8  |  10k  |    0.089(7)   |     108.1(9)      |   0.0003  |  0.0094
|  n6 i1  |  500  |    1.03(5)    |     3.13(3)×10³   |   0.0019  |  0.0578
|  n6 i4  |  500  |    0.038(8)   |     3.80(4)       |   0.0021  |  0.0569
|  n6 i5  |  10k  |    1.27(4)    |     3.13(1)×10³   |   0.0022  |  0.0583
|  n6 i7  |  1k   |    0.027(6)   |     3.80(5)       |   0.0037  |  0.1110
|  n6 i8  |  1k   |    0.020(3)   |     3.79(2)       |   0.0037  |  0.1119
|  n6 i9  |  10k  |    0.028(6)   |     3.79(3)       |   0.0037  |  0.1208
|  n6 i9  |  1k   |    0.029(7)   |     3.80(7)       |   0.0037  |  0.1233
|  n7 i8  |  10k  |    0.062(9)   |     45.7(5)       |   0.0038  |  0.1153
|  n7 i8  |  1k   |    0.062(8)   |     45.8(5)       |   0.0200  |  0.6324
|  n7 i9  |  500  |    0.06(1)    |     45.9(7)       |   0.0204  |  0.6457

Se puede observar que la implementación de diferentes algoritmos lleva a
distintos tiempos de ejecución, también la diferencia con los tiempos se debe a
que la estructura de las cláusulas es muy densa \[5\] en estos problemas (cada
cláusula puede tener suficientes variables para actuar como un cuello de
botella)

## Referencias

\[1\] da Silva, P. F. M. (2010). _Max-SAT algorithms for real world instances_.
(Master Dissertation).

\[2\] Pearl, J. (2009). _Causality: Models, Reasoning, and Inference_ (2nd ed.).
Cambridge University Press.

\[3\] Hyttinen, A., Hoyer, P. O., Eberhardt, F., & Jarvisalo, M. (2013).
_Discovering cyclic causal models with latent variables: A general SAT-based
procedure_. arXiv preprint arXiv:1309.6836.

\[4\] Hyttinen, A., Eberhardt, F., & Järvisalo, M. (2014, July).
_Constraint-based Causal Discovery: Conflict Resolution with Answer Set
Programming_. In UAI (pp. 340-349).

\[5\] Berg, O. J., Hyttinen, A. J., & Järvisalo, M. J. (2019). _Applications of
MaxSAT in data analysis_. In International Conferences on Theory and
Applications of Satisfiability Testing (pp. 50-64). EasyChair Publications.

\[6\] Avellaneda, F. (2020). _A short description of the solver EvalMaxSAT_.
MaxSAT Evaluation, 8, 364.

\[7\] A. Biere, K. Fazekas, M. Fleury and M. Heisinger, _CaDiCaL, Kissat,
Paracooba, Plingeling and Treengeling Entering the SAT Competition 2020_.

\[8\] A. Morgado, C. Dodaro, and J. Marques-Silva, _Core-guided MaxSAT with soft
cardinality constraints_.

\[9\] Silva, J. M., & Sakallah, K. A. (1996, November). _GRASP-a new search
algorithm for satisfiability_. In Proceedings of International Conference on
Computer Aided Design (pp. 220-227). IEEE.

\[10\] Liang, J. H., Ganesh, V., Poupart, P., & Czarnecki, K. (2016, June).
_Learning rate based branching heuristic for SAT solvers_. In International
Conference on Theory and Applications of Satisfiability Testing (pp. 123-140).
Cham: Springer International Publishing.

\[11\] Dubois, O., André, P., Boufkhad, Y., & Carlier, J. (1996). _Sat versus
unsat_. Johnson and Trick, Second DIMACS Series in Discrete Mathematics and
Theoretical Computer Science, American Mathematical Society, 415-436.

\[12\] Jeroslow, R. G., & Wang, J. (1990). _Solving propositional satisfiability
problems_. Annals of mathematics and Artificial Intelligence, 1(1), 167-187.
