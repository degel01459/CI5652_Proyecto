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

using namespace std;

// Configuración
const int NUM_CORRIDAS = 30; // Requerimiento estadístico

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
  // Para resetear entre corridas
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

  // Reinicia la cláusula para una nueva corrida
  void reset() { satisfaccion = TBool::Unknown; }

  // Actualiza satisfacción y modifica frecuencias dinámicamente
  void setSatisfaccion(const vector<TBool> &variablesGlobales, vector<Conteo> *ptrFrecuencias)
  {
    bool esperanza = false;

    for (const int &v : variables)
    {
      int indice = abs(v) - 1;
      TBool valorVar = variablesGlobales[indice];

      // Si la cláusula se satisface
      if ((v > 0 && valorVar == TBool::True) || (v < 0 && valorVar == TBool::False))
      {
        satisfaccion = TBool::True;
        // Actualizar frecuencias (reducirlas porque la cláusula ya se cumplió)
        actualizarFrecuencias(ptrFrecuencias);
        return;
      }
      if (valorVar == TBool::Unknown)
        esperanza = true;
    }

    if (!esperanza)
    {
      satisfaccion = TBool::False;
      // También actualizamos frecuencias si ya es imposible satisfacerla
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

public:
  Formula() {}
  Formula(const vector<Clausula> &clauses)
  { // Copia profunda para resetear
    clausulas = clauses;
  }

  // Calcula costo (Cláusulas insatisfechas)
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

  // Heurística Constructiva
  void solverConstructivo(vector<TBool> &variablesGlobales, vector<Conteo> frecs)
  {
    // Nota: Recibimos 'frecs' por copia para poder modificarla sin dañar la original del main
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

      moda->pos = -9999; // Sacar de competencia
      moda->neg = -9999;

      for (Clausula &clausula : clausulas)
      {
        if (clausula.getSatisfaccion() == TBool::Unknown && clausula.aparece(idModa))
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

  // Búsqueda Local: Hill Climbing
  void busquedaLocal(vector<TBool> &vars)
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
          break; // First improvement
        }
        else
        {
          vars[idx] = (vars[idx] == TBool::True) ? TBool::False : TBool::True; // Revertir
        }
      }
    }
  }

  // ILS
  void busquedaLocalIterada(vector<TBool> &vars, int maxIteraciones)
  {
    int mejorCosto = calcularCosto(vars);
    vector<TBool> mejorSolucion = vars;

    for (int i = 0; i < maxIteraciones; i++)
    {
      vector<TBool> actual = mejorSolucion; // Empezamos desde la mejor conocida (Elitista)

      // 1. Perturbación (Random k-flip 5%)
      int k = max(1, (int)(vars.size() * 0.05));
      for (int j = 0; j < k; j++)
      {
        int idx = rand() % vars.size();
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
  double sum = accumulate(v.begin(), v.end(), 0.0);
  return sum / v.size();
}
double desviacion(const vector<double> &v, double mean)
{
  if (v.size() <= 1)
    return 0.0;
  double sq_sum = inner_product(v.begin(), v.end(), v.begin(), 0.0);
  double stdev = sqrt(sq_sum / v.size() - mean * mean);
  return stdev;
}

int main(int argc, char const *argv[])
{
  srand(time(0));
  if (argc < 2)
  {
    cout << "Uso: ./solver archivo.cnf" << endl;
    return 1;
  }

  ifstream archivo(argv[1]);
  if (!archivo.is_open())
  {
    cout << "Error abriendo archivo" << endl;
    return 1;
  }

  // 1. Lectura del archivo (solo una vez)
  string linea;
  pair<int, int> datosFormula = {0, 0};
  vector<Clausula> clausulasBase;
  vector<Conteo> frecuenciasBase;

  while (getline(archivo, linea))
  {
    if (linea.empty() || linea[0] == 'c')
      continue;
    if (linea[0] == 'p')
    {
      datosFormula = leerPreambulo(linea);
      frecuenciasBase.resize(datosFormula.first);
      clausulasBase.reserve(datosFormula.second);
    }
    else if (isdigit(linea[0]) || linea[0] == '-')
    {
      clausulasBase.push_back(crearClausula(linea, frecuenciasBase));
    }
  }

  // Vectores para estadísticas
  vector<double> tHeur, tBL, tILS;
  vector<double> cHeur, cBL, cILS;

  cout << "Iniciando " << NUM_CORRIDAS << " corridas. Paciencia..." << endl;
  cout << "Progreso: [";

  // 2. Bucle de 30 corridas
  for (int iter = 0; iter < NUM_CORRIDAS; iter++)
  {
    if (iter % 3 == 0)
    {
      cout << "#" << flush;
    } // Barra de progreso

    // Resetear estado para esta iteración
    vector<TBool> vars = vector<TBool>(datosFormula.first, TBool::Unknown);
    vector<Conteo> frecsIter = frecuenciasBase; // Copia fresca
    Formula problema(clausulasBase);            // Copia fresca de cláusulas

    // --- Heurística ---
    auto start = chrono::high_resolution_clock::now();
    problema.solverConstructivo(vars, frecsIter);
    auto end = chrono::high_resolution_clock::now();
    tHeur.push_back(chrono::duration<double>(end - start).count());
    cHeur.push_back(problema.calcularCosto(vars));

    // --- Búsqueda Local ---
    start = chrono::high_resolution_clock::now();
    problema.busquedaLocal(vars);
    end = chrono::high_resolution_clock::now();
    tBL.push_back(chrono::duration<double>(end - start).count()); // BL es acumulativo al tiempo de heurística
    cBL.push_back(problema.calcularCosto(vars));

    // --- ILS ---
    start = chrono::high_resolution_clock::now();
    problema.busquedaLocalIterada(vars, 20); // 20 iteraciones ILS
    end = chrono::high_resolution_clock::now();
    tILS.push_back(chrono::duration<double>(end - start).count());
    cILS.push_back(problema.calcularCosto(vars));
  }
  cout << "] Fin!" << endl
       << endl;

  // 3. Resultados Estadísticos
  cout << fixed << setprecision(4);
  cout << "===========================================" << endl;
  cout << " RESULTADOS PROMEDIO (" << NUM_CORRIDAS << " corridas)" << endl;
  cout << "===========================================" << endl;
  cout << "Algoritmo      | Tiempo (s) (Desv) | Costo (Insat) (Desv)" << endl;
  cout << "-------------------------------------------" << endl;

  double mTH = promedio(tHeur), dTH = desviacion(tHeur, mTH);
  double mCH = promedio(cHeur), dCH = desviacion(cHeur, mCH);
  cout << "Heuristica     | " << mTH << " (" << dTH << ") | " << mCH << " (" << dCH << ")" << endl;

  double mTB = promedio(tBL), dTB = desviacion(tBL, mTB);
  double mCB = promedio(cBL), dCB = desviacion(cBL, mCB);
  cout << "Busqueda Local | " << mTB << " (" << dTB << ") | " << mCB << " (" << dCB << ")" << endl;

  double mTILS = promedio(tILS), dTILS = desviacion(tILS, mTILS);
  double mCILS = promedio(cILS), dCILS = desviacion(cILS, mCILS);
  cout << "ILS            | " << mTILS << " (" << dTILS << ") | " << mCILS << " (" << dCILS << ")" << endl;
  cout << "===========================================" << endl;

  return 0;
}