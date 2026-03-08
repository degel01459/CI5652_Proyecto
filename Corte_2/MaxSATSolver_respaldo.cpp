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

    for (int idxClausula : apariciones[varIdx])
    {
      if (!clausulas[idxClausula].esSatisfecha(vars))
        costoAntes++;
    }

    vars[varIdx] = (vars[varIdx] == TBool::True) ? TBool::False : TBool::True;

    for (int idxClausula : apariciones[varIdx])
    {
      if (!clausulas[idxClausula].esSatisfecha(vars))
        costoDespues++;
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
        int delta = evaluarFlip(vars, i);
        if (delta < 0)
        {
          vars[i] = (vars[i] == TBool::True) ? TBool::False : TBool::True;
          mejora = true;
          break; // First-Improvement
        }
      }
    }
  }

  void busquedaLocalIterada(vector<TBool> &vars, int maxIteraciones, mt19937 &gen)
  {
    int mejorCosto = calcularCosto(vars);
    vector<TBool> mejorSolucion = vars;
    uniform_int_distribution<> dis(0, vars.size() - 1);

    for (int i = 0; i < maxIteraciones; i++)
    {
      vector<TBool> actual = mejorSolucion;

      // 1. Perturbación (Random k-flip 5%)
      int k = max(1, (int)(vars.size() * 0.05));
      for (int j = 0; j < k; j++)
      {
        int idx = dis(gen);
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

  void busquedaTabu(vector<TBool> &vars, int maxIteraciones, int tenureBase)
  {
    int n = vars.size();
    vector<int> tabuUntil(n, 0);

    vector<TBool> mejorSolucionGlobal = vars;
    int mejorCostoGlobal = calcularCosto(vars);
    int costoActual = mejorCostoGlobal;

    random_device rd;
    mt19937 gen(rd());
    uniform_int_distribution<> disTenure(0, 5);

    for (int iter = 1; iter <= maxIteraciones; iter++)
    {
      int mejorVarIdx = -1;
      int mejorDelta = 1e9;

      for (int i = 0; i < n; i++)
      {
        int delta = evaluarFlip(vars, i);
        int nuevoCosto = costoActual + delta;

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

      if (mejorVarIdx != -1)
      {
        vars[mejorVarIdx] = (vars[mejorVarIdx] == TBool::True) ? TBool::False : TBool::True;
        costoActual += mejorDelta;

        tabuUntil[mejorVarIdx] = iter + tenureBase + disTenure(gen);

        if (costoActual < mejorCostoGlobal)
        {
          mejorCostoGlobal = costoActual;
          mejorSolucionGlobal = vars;
        }
      }
    }
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
    double T_min = 0.01;

    uniform_int_distribution<> varDist(0, n - 1);
    uniform_real_distribution<> probDist(0.0, 1.0);

    while (T > T_min)
    {
      for (int i = 0; i < iterPorTemp; i++)
      {
        int idx = varDist(gen);
        int delta = evaluarFlip(actual, idx);
        int nuevoCosto = costoActual + delta;

        if (delta < 0)
        {
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
          double probabilidad = exp(-delta / T);
          if (probDist(gen) < probabilidad)
          {
            actual[idx] = (actual[idx] == TBool::True) ? TBool::False : TBool::True;
            costoActual = nuevoCosto;
          }
        }
      }
      T *= alpha;
    }
    vars = mejorSolucionGlobal;
  }

  void construccionGRASP(vector<TBool> &vars, vector<Conteo> frecs, double alpha, mt19937 &gen)
  {
    int n = vars.size();
    int variablesPendientes = n;

    vector<int> rcl;
    rcl.reserve(n);

    vector<TBool> estadoClausulas(clausulas.size(), TBool::Unknown);

    while (variablesPendientes > 0)
    {
      int s_min = numeric_limits<int>::max();
      int s_max = numeric_limits<int>::min();

      for (int j = 0; j < n; j++)
      {
        if (vars[j] == TBool::Unknown)
        {
          int beneficio = max(frecs[j].pos, frecs[j].neg);
          if (beneficio < s_min)
            s_min = beneficio;
          if (beneficio > s_max)
            s_max = beneficio;
        }
      }

      double umbral = s_max - alpha * (s_max - s_min);
      rcl.clear();

      for (int j = 0; j < n; j++)
      {
        if (vars[j] == TBool::Unknown)
        {
          int beneficio = max(frecs[j].pos, frecs[j].neg);
          if (beneficio >= umbral)
          {
            rcl.push_back(j);
          }
        }
      }

      uniform_int_distribution<> dis(0, rcl.size() - 1);
      int idElegido = rcl[dis(gen)];

      bool valor = frecs[idElegido].pos >= frecs[idElegido].neg;
      vars[idElegido] = valor ? TBool::True : TBool::False;

      frecs[idElegido].pos = -99999;
      frecs[idElegido].neg = -99999;

      for (size_t c = 0; c < clausulas.size(); c++)
      {
        if (estadoClausulas[c] == TBool::Unknown && clausulas[c].aparece(idElegido))
        {
          const auto &variablesC = clausulas[c].getVariables();
          for (int v : variablesC)
          {
            int idx = abs(v) - 1;
            if (idx == idElegido)
            {
              if ((v > 0 && vars[idx] == TBool::True) || (v < 0 && vars[idx] == TBool::False))
              {
                estadoClausulas[c] = TBool::True;
                for (int varClausula : variablesC)
                {
                  if (vars[abs(varClausula) - 1] == TBool::Unknown)
                  {
                    if (varClausula > 0)
                      frecs[abs(varClausula) - 1].pos--;
                    else
                      frecs[abs(varClausula) - 1].neg--;
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

  void busquedaGRASP(vector<TBool> &vars, int maxIteraciones, double alpha, mt19937 &gen, const vector<Conteo> &frecsOriginales)
  {
    int mejorCostoGlobal = numeric_limits<int>::max();
    vector<TBool> mejorSolucionGlobal;

    for (int i = 0; i < maxIteraciones; i++)
    {
      vector<TBool> actual(vars.size(), TBool::Unknown);
      construccionGRASP(actual, frecsOriginales, alpha, gen);
      busquedaLocal(actual);

      int costoFinal = calcularCosto(actual);
      if (costoFinal < mejorCostoGlobal)
      {
        mejorCostoGlobal = costoFinal;
        mejorSolucionGlobal = actual;
      }
    }
    vars = mejorSolucionGlobal;
  }

  void algoritmoGenetico(vector<TBool> &vars, int tamanoPoblacion, int maxGeneraciones, double probCruce, double probMutacion, mt19937 &gen)
  {
    int n = vars.size();
    vector<vector<TBool>> poblacion(tamanoPoblacion, vector<TBool>(n));
    uniform_int_distribution<> disBool(0, 1);

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

    uniform_real_distribution<> disProb(0.0, 1.0);
    uniform_int_distribution<> disTorneo(0, tamanoPoblacion - 1);

    auto seleccionTorneo = [&](const vector<int> &costos, int k)
    {
      int mejorIdx = disTorneo(gen);
      for (int i = 1; i < k; ++i)
      {
        int candidato = disTorneo(gen);
        if (costos[candidato] < costos[mejorIdx])
        {
          mejorIdx = candidato;
        }
      }
      return mejorIdx;
    };

    auto mutar = [&](vector<TBool> &ind)
    {
      for (int j = 0; j < n; j++)
      {
        if (disProb(gen) < probMutacion)
        {
          ind[j] = (ind[j] == TBool::True) ? TBool::False : TBool::True;
        }
      }
    };

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

      int idxMejorActual = distance(costos.begin(), min_element(costos.begin(), costos.end()));
      nuevaPoblacion[0] = poblacion[idxMejorActual];

      for (int i = 1; i < tamanoPoblacion; i += 2)
      {
        int p1 = seleccionTorneo(costos, 3);
        int p2 = seleccionTorneo(costos, 3);

        vector<TBool> hijo1 = poblacion[p1];
        vector<TBool> hijo2 = poblacion[p2];

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

        mutar(hijo1);
        mutar(hijo2);

        nuevaPoblacion[i] = hijo1;
        if (i + 1 < tamanoPoblacion)
        {
          nuevaPoblacion[i + 1] = hijo2;
        }
      }
      poblacion = move(nuevaPoblacion);
    }
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
                                   [media](double acc, double val)
                                   {
                                     return acc + pow(val - media, 2);
                                   });
  return sqrt(sumCuadrados / (v.size() - 1));
}

string formatearMedida(double media, double desv_est)
{
  if (desv_est <= 0)
    return to_string(media) + "(0)";

  int exponente = floor(log10(desv_est));
  int cifra = round(desv_est / pow(10, exponente));

  if (cifra == 10)
  {
    cifra = 1;
    exponente++;
  }

  stringstream ss;
  if (exponente < 0)
  {
    ss << fixed << setprecision(abs(exponente)) << media;
  }
  else
  {
    double factor = pow(10, exponente);
    ss << fixed << setprecision(0) << round(media / factor) * factor;
  }

  return ss.str() + "(" + to_string(cifra) + ")";
}

int main(int argc, char const *argv[])
{
  ios_base::sync_with_stdio(false);
  cin.tie(NULL);

  if (argc < 2)
  {
    cout << "Uso: ./solver archivo1.cnf [archivo2.cnf ...]" << endl;
    return 1;
  }

  cout << "=================================================================================================================================" << endl;
  cout << "                           REPORTE COMPARATIVO: Heuristica, LS, ILS, Tabu, SA, GRASP, Genetico                                   " << endl;
  cout << "=================================================================================================================================" << endl;

  cout << left << setw(18) << "Archivo"
       << "| " << setw(9) << "Exacto"
       << "| " << setw(10) << "Costo H"
       << "| " << setw(10) << "T. H(s)"
       << "| " << setw(10) << "Costo LS"
       << "| " << setw(10) << "T. LS(s)"
       << "| " << setw(11) << "Costo ILS"
       << "| " << setw(11) << "T. ILS(s)"
       << "| " << setw(10) << "Costo TS"
       << "| " << setw(10) << "T. TS(s)"
       << "| " << setw(10) << "Costo SA"
       << "| " << setw(10) << "T. SA(s)"
       << "| " << setw(11) << "Costo GRASP"
       << "| " << setw(11) << "T. GRASP(s)"
       << "| " << setw(10) << "Costo AG"
       << "| " << setw(10) << "T. AG(s)" << endl;
  cout << "---------------------------------------------------------------------------------------------------------------------------------" << endl;

#pragma omp parallel for schedule(dynamic)
  for (int f = 1; f < argc; f++)
  {
    string nombreArchivo = argv[f];

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
      size_t inicio = linea.find_first_not_of(" \t\r\n");
      if (inicio == string::npos)
        continue;

      char primerCaracter = linea[inicio];

      if (primerCaracter == 'c')
        continue;
      else if (primerCaracter == 'p')
      {
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

    vector<double> tH, tLS, tILS, tTS, tSA, tGRASP, tAG;
    int cH, cLS, cILS, cTS, cSA, cGRASP, cAG;

    for (int iter = 0; iter < NUM_CORRIDAS; iter++)
    {
      // 1. HEURISTICA CONSTRUCTIVA (Base)
      vector<TBool> vars = vector<TBool>(datosFormula.first, TBool::Unknown);
      vector<Conteo> frecs = frecuenciasBase;
      Formula problema(clausulasBase, datosFormula.first);
      auto start = chrono::high_resolution_clock::now();
      problema.solverConstructivo(vars, frecs);
      auto end = chrono::high_resolution_clock::now();
      cH = problema.calcularCosto(vars);
      tH.push_back(chrono::duration<double>(end - start).count());

      vector<TBool> varsParaLS = vars;
      vector<TBool> varsParaILS = vars;
      vector<TBool> varsParaTS = vars;
      vector<TBool> varsParaSA = vars;
      vector<TBool> varsParaGRASP(datosFormula.first, TBool::Unknown);
      vector<TBool> varsParaAG = vars;

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
      int tenure = 7 + (datosFormula.first / 10);
      problema.busquedaTabu(varsParaTS, 100, tenure);
      end = chrono::high_resolution_clock::now();
      tTS.push_back(chrono::duration<double>(end - start).count());
      cTS = problema.calcularCosto(varsParaTS);

      // 5. RECOCIDO SIMULADO
      start = chrono::high_resolution_clock::now();
      problema.recocidoSimulado(varsParaSA, gen, 10.0, 0.98, 10);
      end = chrono::high_resolution_clock::now();
      tSA.push_back(chrono::duration<double>(end - start).count());
      cSA = problema.calcularCosto(varsParaSA);

      // 6. GRASP
      start = chrono::high_resolution_clock::now();
      problema.busquedaGRASP(varsParaGRASP, 5, 0.2, gen, frecuenciasBase);
      end = chrono::high_resolution_clock::now();
      tGRASP.push_back(chrono::duration<double>(end - start).count());
      cGRASP = problema.calcularCosto(varsParaGRASP);

      // 7. ALGORITMO GENÉTICO
      double probMutacion = 1.0 / datosFormula.first;
      start = chrono::high_resolution_clock::now();
      problema.algoritmoGenetico(varsParaAG, 50, 100, 0.85, probMutacion, gen);
      end = chrono::high_resolution_clock::now();
      tAG.push_back(chrono::duration<double>(end - start).count());
      cAG = problema.calcularCosto(varsParaAG);
    }

    double mTH = promedio(tH);
    double mTLS = promedio(tLS);
    double mTILS = promedio(tILS);
    double mTTS = promedio(tTS);
    double mTSA = promedio(tSA);
    double mTGRASP = promedio(tGRASP);
    double mTAG = promedio(tAG);

    double sdTH = desviacionEstandar(tH, mTH);
    double sdTLS = desviacionEstandar(tLS, mTLS);
    double sdTILS = desviacionEstandar(tILS, mTILS);
    double sdTTS = desviacionEstandar(tTS, mTTS);
    double sdTSA = desviacionEstandar(tSA, mTSA);
    double sdTGRASP = desviacionEstandar(tGRASP, mTGRASP);
    double sdTAG = desviacionEstandar(tAG, mTAG);

#pragma omp critical
    {
      string nombreCorto = (nombreArchivo.length() > 16) ? "..." + nombreArchivo.substr(nombreArchivo.length() - 15) : nombreArchivo;

      cout << left << setw(18) << nombreCorto
           << "| " << setw(9) << "-------- " // <-- Aquí colocarás tus resultados de EvalMaxSAT luego
           << "| " << setw(10) << cH
           << "| " << setw(10) << formatearMedida(mTH, sdTH)
           << "| " << setw(10) << cLS
           << "| " << setw(10) << formatearMedida(mTLS, sdTLS)
           << "| " << setw(11) << cILS
           << "| " << setw(11) << formatearMedida(mTILS, sdTILS)
           << "| " << setw(10) << cTS
           << "| " << setw(10) << formatearMedida(mTTS, sdTTS)
           << "| " << setw(10) << cSA
           << "| " << setw(10) << formatearMedida(mTSA, sdTSA)
           << "| " << setw(11) << cGRASP
           << "| " << setw(11) << formatearMedida(mTGRASP, sdTGRASP)
           << "| " << setw(10) << cAG
           << "| " << setw(10) << formatearMedida(mTAG, sdTAG) << endl;
    }
  }

  cout << "=================================================================================================================================" << endl;
  return 0;
}