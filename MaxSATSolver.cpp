/**
 * @file MaxSATSolver.cpp
 * @brief Estructuras y funciones para el solver de MAXSAT. CI5652 EM26.
 * @author Alejandro, Ángel, Francisco, Kevin y Sergio
 * @date 2026/01 - 2026/03
 */

#include <algorithm>
#include <chrono>
#include <fstream>
#include <iostream>
#include <sstream>
#include <string>
#include <vector>
#include <cmath>
#include <numeric>
#include <iomanip>
#include <random> // Para generacion aleatoria segura en hilos
#include <omp.h>  // Para poder usar tu CPU al maximo
#include <mutex>  // Para que el texto no se mezcle en consola
#include <cmath>   // para sqrt, pow

using namespace std;

// Configuración
const int NUM_CORRIDAS = 1;

enum class TBool : int8_t
{
  Unknown = -1,
  False = 0,
  True = 1
};

struct Conteo
{
  int pos = 0;
  int neg = 0;
  void reset()
  {
    pos = 0;
    neg = 0;
  }
};

class Clausula
{
private:
  vector<int> variables;
  TBool satisfaccion;

public:
  Clausula() { satisfaccion = TBool::Unknown; }
  Clausula(vector<int> vars)
  {
    variables = vars;
    satisfaccion = TBool::Unknown;
  }

  const vector<int> &getVariables() const { return variables; }
  TBool getSatisfaccion() const { return satisfaccion; }
  void reset() { satisfaccion = TBool::Unknown; }

  void setSatisfaccion(const vector<TBool> &variablesGlobales, vector<Conteo> *ptrFrecuencias)
  {
    bool esperanza = false;
    for (const int &v : variables)
    {
      int indice = abs(v) - 1;
      TBool valorVar = variablesGlobales[indice];

      if ((v > 0 && valorVar == TBool::True) || (v < 0 && valorVar == TBool::False))
      {
        satisfaccion = TBool::True;
        actualizarFrecuencias(ptrFrecuencias);
        return;
      }
      if (valorVar == TBool::Unknown)
        esperanza = true;
    }

    if (!esperanza)
    {
      satisfaccion = TBool::False;
      actualizarFrecuencias(ptrFrecuencias);
    }
    else
    {
      satisfaccion = TBool::Unknown;
    }
  }

  bool esSatisfecha(const vector<TBool> &variablesGlobales) const
  {
    for (int v : variables)
    {
      int idx = abs(v) - 1;
      if (variablesGlobales[idx] == TBool::Unknown)
        continue;
      if ((v > 0 && variablesGlobales[idx] == TBool::True) ||
          (v < 0 && variablesGlobales[idx] == TBool::False))
      {
        return true;
      }
    }
    return false;
  }

  bool aparece(int variable)
  {
    for (int v : variables)
      if (abs(v) == variable + 1)
        return true;
    return false;
  }

  void actualizarFrecuencias(vector<Conteo> *ptrFrecuencias)
  {
    if (!ptrFrecuencias)
      return;
    for (int &variable : variables)
    {
      if (variable > 0)
        (*ptrFrecuencias)[abs(variable) - 1].pos--;
      else
        (*ptrFrecuencias)[abs(variable) - 1].neg--;
    }
  }
};

class Formula
{
private:
  vector<Clausula> clausulas;
  vector<vector<int>> apariciones;

public:
  Formula() {}
  Formula(const vector<Clausula> &clauses, int numVars)
  {
    clausulas = clauses;
    apariciones.resize(numVars);

    for (size_t i = 0; i < clausulas.size(); i++)
    {
      for (int v : clausulas[i].getVariables())
      {
        apariciones[abs(v) - 1].push_back(i);
      }
    }
  }

  int calcularCosto(const vector<TBool> &vars)
  {
    int costo = 0;
    for (const auto &c : clausulas)
    {
      if (!c.esSatisfecha(vars))
        costo++;
    }
    return costo;
  }

  int evaluarFlip(vector<TBool> &vars, int varIdx)
  {
    int costoAntes = 0;
    int costoDespues = 0;

    for (int idxClausula : apariciones[varIdx]) {
      if (!clausulas[idxClausula].esSatisfecha(vars)) costoAntes++;
    }

    vars[varIdx] = (vars[varIdx] == TBool::True) ? TBool::False : TBool::True;

    for (int idxClausula : apariciones[varIdx]) {
      if (!clausulas[idxClausula].esSatisfecha(vars)) costoDespues++;
    }

    vars[varIdx] = (vars[varIdx] == TBool::True) ? TBool::False : TBool::True;

    return costoDespues - costoAntes;
  }

  void solverConstructivo(vector<TBool> &variablesGlobales, vector<Conteo> frecs)
  {
    int variablesPendientes = variablesGlobales.size();
    while (variablesPendientes > 0)
    {
      auto moda = max_element(frecs.begin(), frecs.end(), [](const Conteo &a, const Conteo &b)
                              { return (a.pos + a.neg) < (b.pos + b.neg); });

      if (moda->pos <= 0 && moda->neg <= 0)
        break;

      int idModa = distance(frecs.begin(), moda);
      bool valor = moda->pos >= moda->neg;

      variablesGlobales[idModa] = valor ? TBool::True : TBool::False;
      moda->pos = -9999;
      moda->neg = -9999;

      for (int idxClausula : apariciones[idModa])
      {
        Clausula &clausula = clausulas[idxClausula];
        if (clausula.getSatisfaccion() == TBool::Unknown)
        {
          clausula.setSatisfaccion(variablesGlobales, &frecs);
        }
      }
      variablesPendientes--;
    }

    for (size_t i = 0; i < variablesGlobales.size(); i++)
    {
      if (variablesGlobales[i] == TBool::Unknown)
        variablesGlobales[i] = TBool::False;
    }
  }

  void busquedaLocal(vector<TBool> &vars)
  {
    bool mejora = true;
    while (mejora)
    {
      mejora = false;
      for (int i = 0; i < vars.size(); i++)
      {
        // Preguntamos: ¿Qué pasa si invierto la variable i?
        int delta = evaluarFlip(vars, i);
        if (delta < 0) // Si el delta es negativo, el costo baja (mejora!)
        {
          // Aplicamos el cambio permanentemente
          vars[i] = (vars[i] == TBool::True) ? TBool::False : TBool::True;
          mejora = true;
          break; // First-Improvement: aplicamos la primera mejora que encontramos
        }
      }
    }
  }
  /*void busquedaLocal(vector<TBool> &vars)
  {
    bool mejora = true;
    int costoActual = calcularCosto(vars);

    while (mejora)
    {
      mejora = false;
      vector<int> candidatos;
      for (const auto &c : clausulas)
      {
        if (!c.esSatisfecha(vars))
        {
          for (int v : c.getVariables())
            candidatos.push_back(abs(v) - 1);
        }
      }
      sort(candidatos.begin(), candidatos.end());
      candidatos.erase(unique(candidatos.begin(), candidatos.end()), candidatos.end());

      for (int idx : candidatos)
      {
        vars[idx] = (vars[idx] == TBool::True) ? TBool::False : TBool::True;
        int nuevoCosto = calcularCosto(vars);
        if (nuevoCosto < costoActual)
        {
          costoActual = nuevoCosto;
          mejora = true;
          break;
        }
        else
        {
          vars[idx] = (vars[idx] == TBool::True) ? TBool::False : TBool::True; // Revertir
        }
      }
    }
  }*/

  void busquedaLocalIterada(vector<TBool> &vars, int maxIteraciones, mt19937 &gen)
  {
    int mejorCosto = calcularCosto(vars);
    vector<TBool> mejorSolucion = vars;

    // Distribución uniforme para elegir variables al azar
    uniform_int_distribution<> dis(0, vars.size() - 1);

    for (int i = 0; i < maxIteraciones; i++)
    {
      vector<TBool> actual = mejorSolucion;

      // 1. Perturbación (Random k-flip 5%)
      int k = max(1, (int)(vars.size() * 0.05));
      for (int j = 0; j < k; j++)
      {
        int idx = dis(gen); // Usamos el generador seguro
        actual[idx] = (actual[idx] == TBool::True) ? TBool::False : TBool::True;
      }

      // 2. Búsqueda Local
      busquedaLocal(actual);

      // 3. Aceptación
      int costoActual = calcularCosto(actual);
      if (costoActual < mejorCosto)
      {
        mejorCosto = costoActual;
        mejorSolucion = actual;
      }
    }
    vars = mejorSolucion;
  }

  /**
   * Estructuras y funciones adicionales para Búsqueda Tabú
   */

  // Dentro de la clase Formula:
  void busquedaTabu(vector<TBool> &vars, int maxIteraciones, int tenureBase)
  {
      int n = vars.size();
      // 1. Estructura de Lista Tabú: almacena la iteración hasta la cual la variable está prohibida
      vector<int> tabuUntil(n, 0);

      vector<TBool> mejorSolucionGlobal = vars;
      int mejorCostoGlobal = calcularCosto(vars);
      int costoActual = mejorCostoGlobal;

      // Generador para tenure variable (opcional pero recomendado)
      random_device rd;
      mt19937 gen(rd());
      uniform_int_distribution<> disTenure(0, 5); // Variación de tenure

      for (int iter = 1; iter <= maxIteraciones; iter++)
      {
          int mejorVarIdx = -1;
          int mejorDelta = 1e9; // Buscamos el menor delta (incluso si es positivo/peor)

          // 2. Exploración de la vecindad 1-flip
          for (int i = 0; i < n; i++)
          {
              // Simular el flip
              int delta = evaluarFlip(vars, i);
              int nuevoCosto = costoActual + delta;

              // 3. Lógica de aceptación con Criterio de Aspiración
              bool esTabu = (iter < tabuUntil[i]);
              bool aspira = (nuevoCosto < mejorCostoGlobal);

              if (!esTabu || aspira)
              {
                  if (delta < mejorDelta)
                  {
                      mejorDelta = delta;
                      mejorVarIdx = i;
                  }
              }
          }

          // 4. Ejecutar el mejor movimiento encontrado (aunque sea peor que el actual)
          if (mejorVarIdx != -1)
          {
              vars[mejorVarIdx] = (vars[mejorVarIdx] == TBool::True) ? TBool::False : TBool::True;
              costoActual += mejorDelta;

              // Actualizar lista tabú con tenure variable
              tabuUntil[mejorVarIdx] = iter + tenureBase + disTenure(gen);

              // Actualizar mejor global si aplica
              if (costoActual < mejorCostoGlobal)
              {
                  mejorCostoGlobal = costoActual;
                  mejorSolucionGlobal = vars;
              }
          }
      }
      // Retornar la mejor solución encontrada en todo el proceso
      vars = mejorSolucionGlobal;
  }

  void recocidoSimulado(vector<TBool> &vars, mt19937 &gen, double tempInicial = 10.0, double alpha = 0.95, int iterPorTemp = 100) 
  {
      int n = vars.size();
      vector<TBool> actual = vars;
      vector<TBool> mejorSolucionGlobal = vars;

      int costoActual = calcularCosto(actual);
      int mejorCostoGlobal = costoActual;

      double T = tempInicial;
      double T_min = 0.01; // Temperatura de parada

      // Distribuciones para aleatoriedad
      uniform_int_distribution<> varDist(0, n - 1);
      uniform_real_distribution<> probDist(0.0, 1.0);

      while (T > T_min)
      {
          for (int i = 0; i < iterPorTemp; i++)
          {
              // 1. Elegir un vecino aleatorio (1-flip)
              int idx = varDist(gen);
              int delta = evaluarFlip(actual, idx);
              int nuevoCosto = costoActual + delta;

              // 2. Criterio de aceptación (Metrópolis)
              if (delta < 0)
              {
                  // Mejora directa
                  actual[idx] = (actual[idx] == TBool::True) ? TBool::False : TBool::True;
                  costoActual = nuevoCosto;
                  if (costoActual < mejorCostoGlobal)
                  {
                      mejorCostoGlobal = costoActual;
                      mejorSolucionGlobal = actual;
                  }
              }
              else
              {
                  // Movimiento peor: se acepta con probabilidad e^(-delta / T)
                  double probabilidad = exp(-delta / T);
                  if (probDist(gen) < probabilidad)
                  {
                      costoActual = nuevoCosto;
                  }
              }
          }
          // 3. Enfriamiento
          T *= alpha;
      }
      vars = mejorSolucionGlobal;
  }

  /**
   * Fase de Construcción de GRASP: Greedy Randomized
   * @param alpha Parámetro entre 0 y 1.
   * 0 = Totalmente Greedy, 1 = Totalmente Aleatorio.
   */
  void construccionGRASP(vector<TBool> &vars, vector<Conteo> frecs, double alpha, mt19937 &gen)
  {
    int n = vars.size();
    int variablesPendientes = n;

    // OPTIMIZACIÓN 1: Reutilizar el vector en lugar de crearlo en cada ciclo
    vector<int> rcl;
    rcl.reserve(n);

    // OPTIMIZACIÓN 2: Hacemos una copia local del estado de satisfacción
    // para no corromper el objeto Formula en las 20 iteraciones de GRASP.
    vector<TBool> estadoClausulas(clausulas.size(), TBool::Unknown);

    while (variablesPendientes > 0)
    {
      int s_min = numeric_limits<int>::max();
      int s_max = numeric_limits<int>::min();

      // 1. Encontrar S_min y S_max (sin asignar memoria)
      for (int j = 0; j < n; j++) {
        if (vars[j] == TBool::Unknown) {
          int beneficio = max(frecs[j].pos, frecs[j].neg);
          if (beneficio < s_min) s_min = beneficio;
          if (beneficio > s_max) s_max = beneficio;
        }
      }

      double umbral = s_max - alpha * (s_max - s_min);
      rcl.clear(); // Limpiamos el vector, pero no reasignamos memoria

      // 2. Llenar la RCL
      for (int j = 0; j < n; j++) {
        if (vars[j] == TBool::Unknown) {
          int beneficio = max(frecs[j].pos, frecs[j].neg);
          if (beneficio >= umbral) {
            rcl.push_back(j);
          }
        }
      }

      // 3. Selección aleatoria
      uniform_int_distribution<> dis(0, rcl.size() - 1);
      int idElegido = rcl[dis(gen)];

      bool valor = frecs[idElegido].pos >= frecs[idElegido].neg;
      vars[idElegido] = valor ? TBool::True : TBool::False;

      frecs[idElegido].pos = -99999;
      frecs[idElegido].neg = -99999;

      // 4. EL PASO ADAPTATIVO (Crucial para la calidad de GRASP)
      for (size_t c = 0; c < clausulas.size(); c++) {
        if (estadoClausulas[c] == TBool::Unknown && clausulas[c].aparece(idElegido)) {
          // Verificamos si la cláusula se acaba de satisfacer
          const auto& variablesC = clausulas[c].getVariables();
          for (int v : variablesC) {
            int idx = abs(v) - 1;
            if (idx == idElegido) {
              if ((v > 0 && vars[idx] == TBool::True) || (v < 0 && vars[idx] == TBool::False)) {
                estadoClausulas[c] = TBool::True;
                // Restamos frecuencias de las variables restantes en esta cláusula
                for (int varClausula : variablesC) {
                  if (vars[abs(varClausula) - 1] == TBool::Unknown) {
                    if (varClausula > 0) frecs[abs(varClausula) - 1].pos--;
                    else frecs[abs(varClausula) - 1].neg--;
                  }
                }
                break;
              }
            }
          }
        }
      }
      variablesPendientes--;
    }
  }
  /*void construccionGRASP(vector<TBool> &vars, vector<Conteo> frecs, double alpha, mt19937 &gen)
  {
      int n = vars.size();
      vector<int> asignados;

      // Trabajamos con una copia de las cláusulas para simular la construcción
      for (int i = 0; i < n; i++)
      {
          // 1. Encontrar el rango de beneficio (S_min y S_max)
          int s_min = numeric_limits<int>::max();
          int s_max = numeric_limits<int>::min();

          vector<pair<int, int>> candidatos; // {indice_var, beneficio}
          for (int j = 0; j < n; j++) {
              if (vars[j] == TBool::Unknown) {
                  int beneficio = max(frecs[j].pos, frecs[j].neg);
                  s_min = min(s_min, beneficio);
                  s_max = max(s_max, beneficio);
                  candidatos.push_back({j, beneficio});
              }
          }

          // 2. Definir el umbral para la RCL (Lista Restringida de Candidatos)
          // Umbral = S_max - alpha * (S_max - S_min)
          double umbral = s_max - alpha * (s_max - s_min);

          vector<int> rcl;
          for (auto& c : candidatos) {
              if (c.second >= umbral) {
                  rcl.push_back(c.first);
              }
          }

          // 3. Selección aleatoria de la RCL
          uniform_int_distribution<> dis(0, rcl.size() - 1);
          int idElegido = rcl[dis(gen)];

          // Asignar y actualizar (similar a tu solverConstructivo)
          bool valor = frecs[idElegido].pos >= frecs[idElegido].neg;
          vars[idElegido] = valor ? TBool::True : TBool::False;

          // Marcamos como procesada
          frecs[idElegido].pos = -99999;
          frecs[idElegido].neg = -99999;
      }
  }*/

  void busquedaGRASP(vector<TBool> &vars, int maxIteraciones, double alpha, mt19937 &gen, const vector<Conteo>& frecsOriginales) 
  {
      int mejorCostoGlobal = numeric_limits<int>::max();
      vector<TBool> mejorSolucionGlobal;

      for (int i = 0; i < maxIteraciones; i++)
      {
          vector<TBool> actual(vars.size(), TBool::Unknown);
          // Enviamos copia de las frecuencias ya que la fase constructiva las modifica
          construccionGRASP(actual, frecsOriginales, alpha, gen);
          busquedaLocal(actual);

          int costoFinal = calcularCosto(actual);
          if (costoFinal < mejorCostoGlobal) {
              mejorCostoGlobal = costoFinal;
              mejorSolucionGlobal = actual;
          }
      }
      vars = mejorSolucionGlobal;
  }

  void algoritmoGenetico(vector<TBool> &vars, int tamanoPoblacion, int maxGeneraciones, double probCruce, double probMutacion, mt19937 &gen) 
  {
      int n = vars.size();

      // Inicializar población
      vector<vector<TBool>> poblacion(tamanoPoblacion, vector<TBool>(n));
      uniform_int_distribution<> disBool(0, 1);

      // El individuo 0 es la solución inicial proporcionada (ej. de la heurística)
      poblacion[0] = vars;
      for (int i = 1; i < tamanoPoblacion; i++)
      {
          for (int j = 0; j < n; j++)
          {
              poblacion[i][j] = disBool(gen) ? TBool::True : TBool::False;
          }
      }

      vector<TBool> mejorSolucionGlobal = vars;
      int mejorCostoGlobal = calcularCosto(mejorSolucionGlobal);

      // Distribuciones para cruce, mutación y torneo
      uniform_real_distribution<> disProb(0.0, 1.0);
      uniform_int_distribution<> disTorneo(0, tamanoPoblacion - 1);

      // Lambda para Selección por Torneo (minimización de costo)
      auto seleccionTorneo = [&](const vector<int>& costos, int k) {
          int mejorIdx = disTorneo(gen);
          for (int i = 1; i < k; ++i) {
              int candidato = disTorneo(gen);
              if (costos[candidato] < costos[mejorIdx]) {
                  mejorIdx = candidato;
              }
          }
          return mejorIdx;
      };

      // Lambda para Mutación (Bit-flip)
      auto mutar = [&](vector<TBool>& ind) {
          for (int j = 0; j < n; j++) {
              if (disProb(gen) < probMutacion) {
                  ind[j] = (ind[j] == TBool::True) ? TBool::False : TBool::True;
              }
          }
      };

      // Ciclo Evolutivo
      for (int genIdx = 0; genIdx < maxGeneraciones; genIdx++)
      {
          vector<int> costos(tamanoPoblacion);
          for (int i = 0; i < tamanoPoblacion; i++)
          {
              costos[i] = calcularCosto(poblacion[i]);
              if (costos[i] < mejorCostoGlobal)
              {
                  mejorCostoGlobal = costos[i];
                  mejorSolucionGlobal = poblacion[i];
              }
          }

          vector<vector<TBool>> nuevaPoblacion(tamanoPoblacion, vector<TBool>(n));

          // Elitismo: Pasar al mejor de la generación actual directamente a la siguiente
          int idxMejorActual = distance(costos.begin(), min_element(costos.begin(), costos.end()));
          nuevaPoblacion[0] = poblacion[idxMejorActual];

          // Generar el resto de la población
          for (int i = 1; i < tamanoPoblacion; i += 2)
          {
              // Selección
              int p1 = seleccionTorneo(costos, 3);
              int p2 = seleccionTorneo(costos, 3);

              vector<TBool> hijo1 = poblacion[p1];
              vector<TBool> hijo2 = poblacion[p2];

              // Cruce Uniforme
              if (disProb(gen) < probCruce)
              {
                  for (int j = 0; j < n; j++)
                  {
                      if (disProb(gen) < 0.5)
                      {
                          swap(hijo1[j], hijo2[j]);
                      }
                  }
              }

              // Mutación
              mutar(hijo1);
              mutar(hijo2);

              nuevaPoblacion[i] = hijo1;
              if (i + 1 < tamanoPoblacion)
              {
                  nuevaPoblacion[i + 1] = hijo2;
              }
          }
          poblacion = move(nuevaPoblacion); // move para optimizar memoria
      }

      // Retornar la mejor solución encontrada
      vars = mejorSolucionGlobal;
  }

  // Función auxiliar: Selección por Torneo
  int seleccionTorneo(const vector<int>& costosPoblacion, int k, mt19937 &gen) {
    uniform_int_distribution<> distPoblacion(0, costosPoblacion.size() - 1);

    int mejorIdx = distPoblacion(gen);
    int mejorCosto = costosPoblacion[mejorIdx];

    for (int i = 1; i < k; i++) {
      int idx = distPoblacion(gen);
      if (costosPoblacion[idx] < mejorCosto) {
        mejorCosto = costosPoblacion[idx];
        mejorIdx = idx;
      }
    }
    return mejorIdx;
  }

  // Algoritmo Memético Principal
  void algoritmoMemetico(vector<TBool> &vars, int tamPoblacion, int maxGeneraciones, double tiempoLimiteSegundos, int maxGeneracionesSinMejora, mt19937 &gen) 
  {
      int n = vars.size();

      // Parámetros de excelencia basados en la literatura
      double probMutacion = 1.0 / n; // Tasa de mutación ideal (1/N)
      int k_torneo = 2;              // Torneo binario

      uniform_real_distribution<> probDist(0.0, 1.0);
      
      // Iniciar el reloj para el "Time-Boxing" (Comparación justa vs EvalMaxSAT)
      auto tiempoInicio = chrono::high_resolution_clock::now();

      vector<vector<TBool>> poblacion(tamPoblacion, vector<TBool>(n));
      vector<int> costosPoblacion(tamPoblacion);
      
      int mejorCostoGlobal = 1e9; // Empezamos en "infinito"
      vector<TBool> mejorSolucionGlobal = vars;
      int generacionesSinMejora = 0; // Contador para Early Stopping

      // ==========================================
      // FASE A: Inicialización de la Población
      // ==========================================
      for (int i = 0; i < tamPoblacion; i++) {
          // Crear individuo aleatorio
          for (int j = 0; j < n; j++) {
              poblacion[i][j] = (probDist(gen) < 0.5) ? TBool::True : TBool::False;
          }
          
          // ¡El toque memético! Bajamos el individuo al mínimo local antes de evaluar
          busquedaLocal(poblacion[i]);
          
          costosPoblacion[i] = calcularCosto(poblacion[i]);
          
          // Registrar el mejor global
          if (costosPoblacion[i] < mejorCostoGlobal) {
              mejorCostoGlobal = costosPoblacion[i];
              mejorSolucionGlobal = poblacion[i];
          }
      }

      // ==========================================
      // FASE B: Ciclo Evolutivo
      // ==========================================
      for (int genIdx = 0; genIdx < maxGeneraciones; genIdx++) {
          
          // 1. Criterio de Parada por TIEMPO (Corte de guillotina)
          auto tiempoActual = chrono::high_resolution_clock::now();
          chrono::duration<double> tiempoTranscurrido = tiempoActual - tiempoInicio;
          if (tiempoTranscurrido.count() >= tiempoLimiteSegundos) {
              // Descomenta la siguiente línea si quieres ver el aviso en la consola
              // cout << "c [Memetico] Limite de tiempo de " << tiempoLimiteSegundos << "s alcanzado." << endl;
              break; 
          }

          vector<vector<TBool>> nuevaPoblacion(tamPoblacion, vector<TBool>(n));
          vector<int> nuevosCostos(tamPoblacion);

          // 2. Elitismo Estricto: El mejor padre sobrevive siempre intacto
          auto itMejor = min_element(costosPoblacion.begin(), costosPoblacion.end());
          int idxMejorActual = distance(costosPoblacion.begin(), itMejor);
          
          nuevaPoblacion[0] = poblacion[idxMejorActual];
          nuevosCostos[0] = costosPoblacion[idxMejorActual];

          bool mejoraEnEstaGeneracion = false;

          // 3. Reproducción (Llenar el resto de la nueva población)
          for (int i = 1; i < tamPoblacion; i++) {
              
              // A. Selección
              int p1 = seleccionTorneo(costosPoblacion, k_torneo, gen);
              int p2 = seleccionTorneo(costosPoblacion, k_torneo, gen);

              vector<TBool> hijo(n);
              
              // B. Cruce Uniforme y Mutación en un solo barrido (Optimizado para CPU)
              for (int j = 0; j < n; j++) {
                  // Cruce Uniforme (~50% de cada padre)
                  if (probDist(gen) < 0.5) {
                      hijo[j] = poblacion[p1][j];
                  } else {
                      hijo[j] = poblacion[p2][j];
                  }

                  // Mutación (Probabilidad muy baja: 1/N)
                  if (probDist(gen) < probMutacion) {
                      hijo[j] = (hijo[j] == TBool::True) ? TBool::False : TBool::True;
                  }
              }

              // C. Búsqueda Local (El paso Memético que define al algoritmo)
              busquedaLocal(hijo);
              
              // D. Evaluación
              int costoHijo = calcularCosto(hijo);
              nuevaPoblacion[i] = hijo;
              nuevosCostos[i] = costoHijo;

              // E. Revisar si hay nuevo récord histórico
              if (costoHijo < mejorCostoGlobal) {
                  mejorCostoGlobal = costoHijo;
                  mejorSolucionGlobal = hijo;
                  mejoraEnEstaGeneracion = true;
              }
          }

          // 4. Reemplazo Generacional usando move() (Rápido y sin fugas de memoria)
          poblacion = move(nuevaPoblacion);
          costosPoblacion = move(nuevosCostos);

          // 5. Criterio de Parada por EARLY STOPPING (Convergencia)
          if (mejoraEnEstaGeneracion) {
              generacionesSinMejora = 0; // Reiniciamos el contador si hubo éxito
          } else {
              generacionesSinMejora++;
              if (generacionesSinMejora >= maxGeneracionesSinMejora) {
                  // cout << "c [Memetico] Convergencia alcanzada. Sin mejoras por " << maxGeneracionesSinMejora << " generaciones." << endl;
                  break; // Abortamos, ya estamos en un óptimo local del que no podemos salir
              }
          }
      }

      // ==========================================
      // FASE C: Entrega de Resultados
      // ==========================================
      vars = mejorSolucionGlobal;
  }
};

Clausula crearClausula(string linea, vector<Conteo> &frecuencias)
{
  stringstream ss(linea);
  int variable = 0;
  vector<int> variables;
  while (ss >> variable && variable != 0)
  {
    if (variable > 0)
      frecuencias[abs(variable) - 1].pos++;
    else
      frecuencias[abs(variable) - 1].neg++;
    variables.push_back(variable);
  }
  return Clausula(variables);
}

pair<int, int> leerPreambulo(string linea)
{
  stringstream ss(linea);
  string temp;
  ss >> temp >> temp;
  int vars, clausulas;
  ss >> vars >> clausulas;
  return {vars, clausulas};
}

// Funciones estadísticas
double promedio(const vector<double> &v)
{
  if (v.empty())
    return 0.0;
  return accumulate(v.begin(), v.end(), 0.0) / v.size();
}

double desviacionEstandar(const vector<double> &v, const double media)
{
  if (v.size() < 2)
    return 0.0;
  double sumCuadrados = accumulate(v.begin(), v.end(), 0.0,
      [media](double acc, double val) {
          return acc + pow(val - media, 2);
      });
  return sqrt(sumCuadrados / (v.size() - 1));
}

string formatearMedida(double media, double desv_est) {
    if (desv_est <= 0) return to_string(media) + "(0)";

    // 1. Encontrar la posición de la primera cifra significativa de la desv. estándar
    // Log10 nos dice la magnitud. Floor nos da el exponente entero.
    int exponente = floor(log10(desv_est));

    // 2. Extraer esa primera cifra
    // Multiplicamos para que la cifra sea la unidad, y redondeamos.
    int cifra = round(desv_est / pow(10, exponente));

    // Caso especial: si el redondeo nos lleva a 10 (ej: 0.099 -> 0.1)
    if (cifra == 10) {
        cifra = 1;
        exponente++;
    }

    // 3. Formatear la media con la precisión adecuada
    stringstream ss;
    if (exponente < 0) {
        // Si la cifra está en los decimales
        ss << fixed << setprecision(abs(exponente)) << media;
    } else {
        // Si la cifra es mayor o igual a 1, redondeamos a la unidad/decena pertinente
        double factor = pow(10, exponente);
        ss << fixed << setprecision(0) << round(media / factor) * factor;
    }

    // 4. Construir el string final
    return ss.str() + "(" + to_string(cifra) + ")";
}

int main(int argc, char const *argv[])
{
  // Optimizacion de I/O
  ios_base::sync_with_stdio(false);
  cin.tie(NULL);

  if (argc < 2)
  {
    cout << "Uso: ./solver archivo1.cnf [archivo2.cnf ...]" << endl;
    return 1;
  }

  cout << "==========================================================================================================" << endl;
  cout << " REPORTE COMPARATIVO 30 REPETICIONES: Heurística vs Búsqueda Local (LS) vs Búsqueda Local Iterada (ILS)" << endl;
  cout << "==========================================================================================================" << endl;

  // Encabezados ajustados para 3 metodos
  cout << left << setw(18) << "Archivo"
       << "| " << setw(9) << "Exacto"
       << "| " << setw(10) << "Costo H"
       << "| " << setw(10) << "T. H(s)"
       << "| " << setw(10) << "Costo LS"
       << "| " << setw(10) << "T. LS(s)"
       << "| " << setw(11) << "Costo ILS"
       << "| " << setw(11) << "T. ILS(s)"
       << "| " << setw(11) << "Costo TS"
       << "| " << setw(11) << "T. TS(s)"
       << "| " << setw(11) << "Costo SA"
       << "| " << setw(11) << "T. SA(s)"
       /*<< "| " << setw(11) << "Costo GRASP"
       << "| " << setw(11) << "T. GRASP(s)"*/
       << "| " << setw(11) << "Costo AG"
       << "| " << setw(11) << "T. AG(s)"
       /*<< "| " << setw(6) << "Gap H-I%"*/ << endl;
  cout << "----------------------------------------------------------------------------------------------------------" << endl;

#pragma omp parallel for schedule(dynamic)
  for (int f = 1; f < argc; f++)
  {
    string nombreArchivo = argv[f];

    // Semilla unica por hilo
    size_t seed = hash<string>{}(nombreArchivo) + omp_get_thread_num();
    mt19937 gen(seed);

    ifstream archivo(nombreArchivo);
    if (!archivo.is_open())
      continue;

    string linea;
    pair<int, int> datosFormula = {0, 0};
    vector<Clausula> clausulasBase;
    vector<Conteo> frecuenciasBase;

    while (getline(archivo, linea))
    {
      // 1. Buscar el primer carácter que no sea espacio o tabulación
      size_t inicio = linea.find_first_not_of(" \t\r\n");

      // Si la línea está vacía o solo tiene espacios, la saltamos
      if (inicio == string::npos)
        continue;

      // 2. Tomar el primer carácter real de la línea
      char primerCaracter = linea[inicio];

      if (primerCaracter == 'c')
      {
        continue;
      }
      else if (primerCaracter == 'p')
      {
        // Pasamos la línea completa a tu función
        datosFormula = leerPreambulo(linea);
        frecuenciasBase.resize(datosFormula.first);
        clausulasBase.reserve(datosFormula.second);
      }
      else if (isdigit(primerCaracter) || primerCaracter == '-')
      {
        clausulasBase.push_back(crearClausula(linea, frecuenciasBase));
      }
    }
    archivo.close();

    // Vectores para guardar promedios de los metodos
    vector<double> tH, tLS, tILS, tTS, tSA, /*tGRASP,*/ tAG, tAM;
    int cH, cLS, cILS, cTS, cSA, /*cGRASP,*/ cAG, cAM;

    for (int iter = 0; iter < NUM_CORRIDAS; iter++)
    {
      // 1. HEURISTICA CONSTRUCTIVA (Base)
      vector<TBool> vars = vector<TBool>(datosFormula.first, TBool::Unknown);
      vector<Conteo> frecs = frecuenciasBase;
      Formula problema(clausulasBase, datosFormula.first);
      auto start = chrono::high_resolution_clock::now();
      problema.solverConstructivo(vars, frecs); // Construimos solucion inicial
      auto end = chrono::high_resolution_clock::now();
      cH = problema.calcularCosto(vars);
      tH.push_back(chrono::duration<double>(end - start).count());

      // Copiamos la solucion de la heuristica para usarla en LS y en ILS por separado
      vector<TBool> varsParaLS = vars;
      vector<TBool> varsParaILS = vars;
      vector<TBool> varsParaTS = vars;
      vector<TBool> varsParaSA = vars;
      // Usamos una copia limpia de las variables (GRASP construye su propia solución)
      vector<TBool> varsParaGRASP(datosFormula.first, TBool::Unknown);
      vector<TBool> varsParaAG = vars;
      vector<TBool> varsParaAM = vars;

      // 2. BUSQUEDA LOCAL
      start = chrono::high_resolution_clock::now();
      problema.busquedaLocal(varsParaLS);
      end = chrono::high_resolution_clock::now();
      tLS.push_back(chrono::duration<double>(end - start).count());
      cLS = problema.calcularCosto(varsParaLS);

      // 3. BUSQUEDA LOCAL ITERADA
      start = chrono::high_resolution_clock::now();
      problema.busquedaLocalIterada(varsParaILS, 20, gen);
      end = chrono::high_resolution_clock::now();

      tILS.push_back(chrono::duration<double>(end - start).count());
      cILS = problema.calcularCosto(varsParaILS);

      // 4. BUSQUEDA TABU
      start = chrono::high_resolution_clock::now();
      int tenure = 7 + (datosFormula.first / 10); // Ejemplo de tenure proporcional
      problema.busquedaTabu(varsParaTS, 100, tenure);
      end = chrono::high_resolution_clock::now();
      tTS.push_back(chrono::duration<double>(end - start).count());
      cTS = problema.calcularCosto(varsParaTS);

      // 5. RECOCIDO SIMULADO
      start = chrono::high_resolution_clock::now();
      // Parámetros sugeridos: temp inicial 10, enfriamiento 0.98, 10 iter por nivel
      problema.recocidoSimulado(varsParaSA, gen, 10.0, 0.98, 10);
      end = chrono::high_resolution_clock::now();
      tSA.push_back(chrono::duration<double>(end - start).count());
      cSA = problema.calcularCosto(varsParaSA);

      // 6. GRASP
      /*start = chrono::high_resolution_clock::now();
      // Parámetros: 5 iteraciones, alpha = 0.2 (20% de aleatoriedad en RCL)
      problema.busquedaGRASP(varsParaGRASP, 5, 0.2, gen, frecuenciasBase);
      end = chrono::high_resolution_clock::now();

      tGRASP.push_back(chrono::duration<double>(end - start).count());
      cGRASP = problema.calcularCosto(varsParaGRASP);*/

      // 7. ALGORITMO GENÉTICO
      // Parámetros: Población=50, Generaciones=100, Cruce=0.85, Mutación=1.0/vars
      double probMutacion = 1.0 / datosFormula.first;

      start = chrono::high_resolution_clock::now();
      problema.algoritmoGenetico(varsParaAG, 50, 100, 0.85, probMutacion, gen);
      end = chrono::high_resolution_clock::now();

      tAG.push_back(chrono::duration<double>(end - start).count());
      cAG = problema.calcularCosto(varsParaAG);

      // 8. ALGORITMO MEMÉTICO
      // Parámetros: Población=20, Torneo=2, Mutación=1.0/vars Parada=60s o 30 generaciones
      start = chrono::high_resolution_clock::now();
      problema.algoritmoMemetico(varsParaAM, 50, 1000, 60.0, 25, gen);
      end = chrono::high_resolution_clock::now();

      tAM.push_back(chrono::duration<double>(end - start).count());
      cAM = problema.calcularCosto(varsParaAM);
    }

    // Promedios
    //double mCH = promedio(cH);
    double mTH = promedio(tH);
    //double mCLS = promedio(cLS);
    double mTLS = promedio(tLS);
    //double mCILS = promedio(cILS);
    double mTILS = promedio(tILS);
    //double mCTS = promedio(cTS);
    double mTTS = promedio(tTS);
    //double mCSA = promedio(cSA);
    double mTSA = promedio(tSA);
    /*//double mCGRASP = promedio(cGRASP);
    double mTGRASP = promedio(tGRASP);*/
    //double mCAG = promedio(cAG);
    double mTAG = promedio(tAG);
    double mTAM = promedio(tAM);

    // DesviacionesEstandar
    //double sdCH = desviacionEstandar(cH, mCH);
    double sdTH = desviacionEstandar(tH, mTH);
    //double sdCLS = desviacionEstandar(cLS, mCLS);
    double sdTLS = desviacionEstandar(tLS, mTLS);
    //double sdCILS = desviacionEstandar(cILS, mCILS);
    double sdTILS = desviacionEstandar(tILS, mTILS);
    //double sdCTS = desviacionEstandar(cTS, mCTS);
    double sdTTS = desviacionEstandar(tTS, mTTS);
    //double sdCSA = desviacionEstandar(cSA, mCSA);
    double sdTSA = desviacionEstandar(tSA, mTSA);
    /*//double sdCGRASP = desviacionEstandar(cGRASP, mCGRASP);
    double sdTGRASP = desviacionEstandar(tGRASP, mTGRASP);*/
    //double sdCAG = desviacionEstandar(cAG, mCAG);
    double sdTAG = desviacionEstandar(tAG, mTAG);
    double sdTAM = desviacionEstandar(tAM, mTAM);

    // Mejora total (Heuristica vs ILS)
    //double mejora = (mCH > 0) ? ((mCH - mCILS) / mCH) * 100.0 : 0.0;*/

#pragma omp critical
    {
      // cout << fixed << setprecision(2);
      string nombreCorto = (nombreArchivo.length() > 16) ? "..." + nombreArchivo.substr(nombreArchivo.length() - 30) : nombreArchivo;

      cout << left << setw(18) << nombreCorto
           << "| " << setw(9) << "-------- "
           << "| " << setw(10) << cH //formatearMedida(int(mCH), int(sdCH))
           << "| " << setw(10) << formatearMedida(mTH, sdTH)
           << "| " << setw(10) << cLS //formatearMedida(mCLS, sdCLS)
           << "| " << setw(10) << formatearMedida(mTLS, sdTLS)
           << "| " << setw(11) << cILS //formatearMedida(mCILS, sdCILS)
           << "| " << setw(11) << formatearMedida(mTILS, sdTILS)
           << "| " << setw(11) << cTS //formatearMedida(mCTS, sdCTS)
           << "| " << setw(11) << formatearMedida(mTTS, sdTTS)
           << "| " << setw(11) << cSA //formatearMedida(mCSA, sdCSA)
           << "| " << setw(11) << formatearMedida(mTSA, sdTSA)
           /*<< "| " << setw(11) << cGRASP //formatearMedida(mCGRASP, sdCGRASP)
           << "| " << setw(11) << formatearMedida(mTGRASP, sdTGRASP)*/
           << "| " << setw(11) << cAG //formatearMedida(mCAG, sdCAG)
           << "| " << setw(11) << formatearMedida(mTAG, sdTAG)
           << "| " << setw(11) << cAM //formatearMedida(mCAM, sdCAM)
           << "| " << setw(11) << formatearMedida(mTAM, sdTAM)
           /*<< "| " << setw(5) << mejora << "%"*/ << endl;
    }
  }

  cout << "==========================================================================================================" << endl;
  return 0;
}
