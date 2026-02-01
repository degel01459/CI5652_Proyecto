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

using namespace std;

// Configuración
const int NUM_CORRIDAS = 30;

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

public:
  Formula() {}
  Formula(const vector<Clausula> &clauses) { clausulas = clauses; }

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
          break;
        }
        else
        {
          vars[idx] = (vars[idx] == TBool::True) ? TBool::False : TBool::True; // Revertir
        }
      }
    }
  }

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
  cout << left << setw(35) << "Archivo"
       << "| " << setw(8) << "Costo H"
       << "| " << setw(8) << "T. H(s)"
       << "| " << setw(8) << "Costo LS"
       << "| " << setw(8) << "T. LS(s)"
       << "| " << setw(9) << "Costo ILS"
       << "| " << setw(9) << "T. ILS(s)"
       << "| " << setw(6) << "Gap H-I%" << endl;
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
    archivo.close();

    // Vectores para guardar promedios de los 3 metodos
    vector<double> tH, tLS, tILS;
    vector<double> cH, cLS, cILS;

    for (int iter = 0; iter < NUM_CORRIDAS; iter++)
    {
      // 1. HEURISTICA CONSTRUCTIVA (Base)
      vector<TBool> vars = vector<TBool>(datosFormula.first, TBool::Unknown);
      vector<Conteo> frecs = frecuenciasBase;
      Formula problema(clausulasBase);

      auto start = chrono::high_resolution_clock::now();
      problema.solverConstructivo(vars, frecs); // Construimos solucion inicial
      auto end = chrono::high_resolution_clock::now();

      double costoH = problema.calcularCosto(vars);
      tH.push_back(chrono::duration<double>(end - start).count());
      cH.push_back(costoH);

      // Copiamos la solucion de la heuristica para usarla en LS y en ILS por separado
      vector<TBool> varsParaLS = vars;
      vector<TBool> varsParaILS = vars;

      // 2. BUSQUEDA LOCAL
      start = chrono::high_resolution_clock::now();
      problema.busquedaLocal(varsParaLS);
      end = chrono::high_resolution_clock::now();

      tLS.push_back(chrono::duration<double>(end - start).count());
      cLS.push_back(problema.calcularCosto(varsParaLS));

      // 3. BUSQUEDA LOCAL ITERADA
      start = chrono::high_resolution_clock::now();
      problema.busquedaLocalIterada(varsParaILS, 20, gen);
      end = chrono::high_resolution_clock::now();

      tILS.push_back(chrono::duration<double>(end - start).count());
      cILS.push_back(problema.calcularCosto(varsParaILS));
    }

    // Promedios
    double mCH = promedio(cH);
    double mTH = promedio(tH);
    double mCLS = promedio(cLS);
    double mTLS = promedio(tLS);
    double mCILS = promedio(cILS);
    double mTILS = promedio(tILS);

    // Mejora total (Heuristica vs ILS)
    double mejora = (mCH > 0) ? ((mCH - mCILS) / mCH) * 100.0 : 0.0;

#pragma omp critical
    {
      cout << fixed << setprecision(2);
      string nombreCorto = (nombreArchivo.length() > 33) ? "..." + nombreArchivo.substr(nombreArchivo.length() - 30) : nombreArchivo;

      cout << left << setw(35) << nombreCorto
           << "| " << setw(8) << mCH
           << "| " << setw(8) << setprecision(4) << mTH
           << "| " << setw(8) << setprecision(2) << mCLS
           << "| " << setw(8) << setprecision(4) << mTLS
           << "| " << setw(9) << setprecision(2) << mCILS
           << "| " << setw(9) << setprecision(4) << mTILS
           << "| " << setw(5) << setprecision(1) << mejora << "%" << endl;
    }
  }

  cout << "==========================================================================================================" << endl;
  return 0;
}