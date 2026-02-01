/**
 * @file MaxSATSolver.cpp
 * @brief Estructuras y funciones para el solver de MAXSAT. CI5652 EM26.
 * @author Alejandro, Ángel, Francisco, Kevin y Sergio
 * @date 2026/01 - 2026/03
 */

#include <algorithm>
#include <cctype>
#include <chrono>
// #include <charconv>
#include <fstream>
#include <iostream>
#include <iterator>
// #include <map>
#include <ranges>
#include <sstream>
#include <string>
#include <utility>
#include <vector>

using namespace std;

/**
 * @brief Estados válidos para un literal.
 * 
 * Este enum define los estados válidos para un literal de una fórmula de SAT.
 */
enum class TBool : int8_t {
  Unknown = -1, ///< No se le ha asignado ningún valor.
  False,        ///< Se le asignó el valor False.
  True          ///< Se le asignó el valor True.
};

/**
 * @brief Contador de apariciones de un literal dentro de una fórmula.
 * 
 * Estructura de datos pura para llevar la frecuencia de un literal.
 */
struct Conteo {
  int pos = 0; ///< Cuenta las veces que el literal aparece sin negar.
  int neg = 0; ///< Cuenta las veces que el literal aparece negado.

  bool operator<(const Conteo& otro) const {
    return (this->pos + this->neg) < (otro.pos + otro.neg);
  }
};

/**
 * @brief Clase principal para el control de los literales en una fórmula.
 * 
 * Esta clase mantiene un vector de variables y su valor de satisfacción actual.
 */
class Clausula {
  private:
    vector<int> variables; ///< Lleva las variables presentes en la cláusula.
    TBool satisfaccion;    ///< Mantiene el estado de verdad de la cláusula.
    vector<Conteo>* ptrFrecuencias; ///< Apuntador al vector de frecuencias.

  public:
    /**
     * @brief Constructor por defecto.
     * Inicializa la Clausula sin variables y con valor de verdad sin asignar.
     */
    Clausula() {satisfaccion = TBool::Unknown;}

    /**
     * @brief Constructor que recibe la lista de variables de la cláusula y un
     * apuntador al vector global de frecuencias.
     * @param vars Variables que aparecen en la cláusula.
     * @param frecs Apuntados al vector global de frecuencias.
     * Inicializa la Clausula con estas variables y no asigna valor de verdad.
     * Almacena la dirección del vector de frecuencias para actualizarlo en su
     * momento.
     */
    Clausula(vector<int> vars, vector<Conteo>* frecs) {
      variables = vars;
      satisfaccion = TBool::Unknown;
      ptrFrecuencias = frecs;
    }


    /**
     * @brief Obtiene el estado de satisfacción actual.
     * @return El valor actual del @ref satisfaccion.
     */
    TBool getSatisfaccion() {return satisfaccion;}

    /**
     * @brief Calcula el estado de satisfacción actual.
     * @param variablesGlobales Vector global de variables.
     * Consulta el valor global de cada variable en la cláusula y calcula la
     * la satisfacción. Despliega el método que actualiza las frecuencias de las
     * variables.
     */
    void setSatisfaccion(const vector<TBool>& variablesGlobales) {
      bool esperanza = false;
      switch(satisfaccion) {
        case TBool::Unknown:
          for (const int &v : variables) {
            TBool valorV = variablesGlobales[abs(v) - 1];
            if (v > 0) {
              switch(valorV) {
                case TBool::True:
                  satisfaccion = TBool::True;
                  actualizarFrecuencias();
                  return;
                default:
                  break;
              }
            } else {
              switch(valorV) {
                case TBool::False:
                  satisfaccion = TBool::True;
                  actualizarFrecuencias();
                  return;
                default:
                  break;
              }
            }
            if (valorV == TBool::Unknown) {
              esperanza = true;
            }
          }
          if (!esperanza) {
            satisfaccion = TBool::False;
            actualizarFrecuencias();
          }
          return;
        default: return;
      }
    }

    /**
     * @brief Función auxiliar que indica si una variable ocurre en la cláusula.
     * 
     * @param variable Identificador de la variable buscada.
     * @return true si y solo si la variable se encuentra en la cláusula.
     */
    bool aparece(int variable) {
      auto pos = find(variables.begin(), variables.end(), variable + 1);
      return pos != variables.end();
    }

    /**
     * @brief Función auxiliar que actualiza las frecuencias de las variables
     * en el vector global de frecuencias.
     */
    void actualizarFrecuencias() {
      if (!ptrFrecuencias) {
        return;
      }
      for (int &variable : variables) {
        if (variable > 0) {
          (*ptrFrecuencias)[abs(variable) - 1].pos--;  
        } else {
          (*ptrFrecuencias)[abs(variable) - 1].neg--;  
        }
      }
    }
};

/**
 * @brief Clase principal para el control de las cláusulas en una fórmula.
 * 
 * Esta clase mantiene un vector de cláusulas y el número de cláusulas que no
 * se pueden satisfacer.
 */
class Formula {
  private:
    vector<Clausula> clausulas; ///< Lleva las cláusulas de la fórmula.
    int insatisfechas;          ///< Cuenta las cláusulas insatisfechas.
    int pendientes;          ///< Cuenta las cláusulas insatisfechas.

  public:
    /**
     * @brief Constructor por defecto.
     * Inicializa la Formula sin cláusulas y sin cláusulas insatisfechas.
     */
    Formula() {insatisfechas = pendientes = 0;}

    /**
     * @brief Constructor que recibe la lista de cláusulas.
     * Inicializa la Formula con estas cláusulas y con cero cláusulas sin
     * satisfacer.
     */
    Formula(vector<Clausula>&& clauses, int claus) {
      clausulas = move(clauses);
      insatisfechas = 0;
      pendientes = claus;
    }

    /**
     * @brief Calcula el número de cláusulas insatisfechas de la fórmula.
     * Revisa cada cláusula y cuenta cada vez que alguna es falsa.
     */
    void actualizarClausulas() {
      insatisfechas = 0;
      pendientes = clausulas.size();
      for (Clausula &c : clausulas) {
        TBool estado = c.getSatisfaccion();
        if (estado != TBool::Unknown) {
          pendientes--;
        }
        if (estado == TBool::False) {
          insatisfechas++;
        }
      }
    }

    /**
     * @brief Obtiene la lista de cláusulas de la fórmula.
     * @return El valor actual de @ref clausulas.
     */
    const vector<Clausula>& getClausulas() const {return clausulas;}

    /**
     * @brief Obtiene el conteo de cláusulas insatisfechas.
     * @return El valor actual del @ref insatisfechas.
     */
    int getInsatisfechas() const {return insatisfechas;}

    /**
     * @brief Función principal que asigna los valores de verdad de las
     * variables y el número de cláusulas insatisfechas.
     * 
     * @param variablesGlobales vector de variables globales.
     * @param variablesGlobales vector global de frecuencias.
     * @return true si y solo si existe un cláusula que sea desconocida.
     */
    void solver(vector<TBool>& variablesGlobales, vector<Conteo>& frecs) {
      int variablesPendientes = variablesGlobales.size();
      do {
        auto moda = max_element(frecs.begin(), frecs.end());
        int idModa = distance(frecs.begin(), moda); // Ya es el índice
        bool valor = moda->pos >= moda->neg;
        variablesGlobales[idModa] = valor ? TBool::True : TBool::False;
        moda->pos = moda->neg = 0;
        for (Clausula &clausula : clausulas) {
          if (clausula.aparece(idModa)) {
            clausula.setSatisfaccion(variablesGlobales);
          }
        }
        actualizarClausulas();
        variablesPendientes--;
      } while (variablesPendientes > 0);
    }

    vector<Clausula> getClausInsatisfechas() {
      vector<Clausula> noSatisfechas;
      copy_if(
        clausulas.begin(), clausulas.end(), back_inserter(noSatisfechas),
        [](Clausula c) {return c.getSatisfaccion() == TBool::False;}
      );
      return noSatisfechas;
    }

  void flip(vector<TBool>& variablesGlobales, int var) {
    int indice = abs(var) - 1;
    switch(variablesGlobales[indice]) {
      case TBool::False;
        test(indice, TBool::True);
        break;
      case TBool::True;
        test(indice, TBool::False);
        break;
      default:
        break;
    }
  }

  void test(int var, TBool valor) {
    for (Clausula &c : clausulas) {
      if (c.getSatisfaccion() == TBool::False && c.aparece(var)) {
        
      }
    }
  }

};

/**
 * @brief Función auxiliar que crea el vector de variables globales.
 * 
 * @param total Número de variables a crear.
 * @return El vector de variables globales incializadas en TBool::Unknown.
 */
vector<TBool> crearVariablesGlobales(int total) {
  vector<TBool> variables;
  variables.assign(total, TBool::Unknown);
  return variables;
}

/**
 * @brief Función auxiliar que crea una cláusula y actualiza el vector de
 * frecuencias a partir de un string en formato DIMACS.
 * 
 * @param linea Texto en formato en DIMACS a traducir.
 * @param frecuencias Vector global de frecuencias.
 * @return Cláusula creada con las variables.
 */
Clausula crearClausula(string linea, vector<Conteo>& frecuencias) {
  stringstream ss(linea);
  int variable = 0;
  vector<int> variables;
  ss >> variable;
  while (variable) {
    if (variable > 0) {
      frecuencias[abs(variable) - 1].pos++;  
    } else {
      frecuencias[abs(variable) - 1].neg++;  
    }
    variables.push_back(variable);
    ss >> variable;
  }
  return Clausula(move(variables), &frecuencias);
}

/**
 * @brief Función auxiliar que lee el preámbulo en formato DIMACS y devuelve
 * como un par de enteros con el número de variables y el número de cláusulas.
 * 
 * @param linea Texto en formato en DIMACS a traducir.
 * @return Par de enteros {variables, cláusulas}.
 */
pair<int, int> leerPreambulo(string linea) {
  stringstream ss(linea.substr(6));
  int vars, clausulas;
  ss >> vars;
  ss >> clausulas;
  return {vars, clausulas};
}

int main(int argc, char const *argv[]) {
  // Se crean las variables del programa
  ifstream archivo(argv[1]);
  int insatisfechas;
  pair<int,int> formula;
  string linea, solucion;
  vector<Clausula> clausulas;
  vector<TBool> variablesGlobales;
  vector<Conteo> frecuencias;

  while (getline(archivo, linea)) {
    if (linea.empty()) {
      continue;
    }
    switch(linea[0]) {
      case 'c':
        break;
      case 'p':
        formula = leerPreambulo(linea);
        variablesGlobales = crearVariablesGlobales(formula.first);
        frecuencias.assign(formula.first, Conteo());
        solucion.reserve(formula.first);
        clausulas.reserve(formula.second);
        break;
      default:
        if (isdigit(linea[0]) || linea[0] == '-') {
          clausulas.push_back(crearClausula(linea, frecuencias));
        } else {
          break;
        }
    }
  }

  auto inicio = chrono::high_resolution_clock::now();
  Formula problema = Formula(move(clausulas), formula.second);
  problema.solver(variablesGlobales, frecuencias);
  auto fin = std::chrono::high_resolution_clock::now();

  chrono::duration<double> duracion = fin - inicio;


  for (int i = 0; i < formula.first; i++) {
    switch(variablesGlobales[i]) {
      case TBool::True:
        solucion[i] = '1';
        break;
      case TBool::False:
        solucion[i] = '0';
        break;
      default:
        solucion[i] = '?';
        break;
    }
  }

  // cout << "Sol.: ";
  // for (int i = 0; i < formula.first; i++) {
  //   cout << solucion[i];
  // }
  insatisfechas = problema.getInsatisfechas();
  cout << "Insatisfechas: " << insatisfechas << endl;
  cout << "Tiempo de ejecución: " << duracion.count() << " s" << endl;

  /* Búsqueda local
  
  [0-1-?]+ --> cambiar 1 *Guardar el valor original que ha cambiado*
  contar insatisfacción
  Mejora? Terminar
  No mejora? Devuelve el cambio y prueba otro.
  Repetir hasta que se completen variablesGlobales.size() / 2 cambios.
  Decir que no se encontró algo mejor.
  */

  clausulas.clear();
  clausulas = move(problema.getClausInsatisfechas());

  if (!problema.getInsatisfechas()) {
    inicio = chrono::high_resolution_clock::now();
    for (Clausula &c : clausulas) {
      for (const int var : const c.getVariables()) {
        problema.flip(var);
        if problema.getInsatisfechas() <= insatisfechas {
          insatisfechas = problema.getInsatisfechas();
        } else {
          switch(variablesGlobales[var]) {
            case TBool::True:
              variablesGlobales[var] = TBool::False;
              break;
            case TBool::False;
              variablesGlobales[var] = TBool::True;
              break;
            default:
              break;
          }
        }
      }
    }
    fin = std::chrono::high_resolution_clock::now();

    duracion = fin - inicio;

    cout << "Búsqueda local" << endl;
    cout << "Insatisfechas: " << insatisfechas << endl;
    cout << "Tiempo de ejecución: " << duracion.count() << " s" << endl;
  }

  /* Búsqueda local iterada vs guiada */
  /*if (!problema.getInsatisfechas()) {
    inicio = chrono::high_resolution_clock::now();
    for (Clausula &c : clausulas) {
      for (const int var : const c.getVariables()) {
        problema.perturbar(var);
        if problema.getInsatisfechas() <= insatisfechas {
          insatisfechas = problema.getInsatisfechas();
        } else {
          switch(variablesGlobales[var]) {
            case TBool::True:
              variablesGlobales[var] = TBool::False;
              break;
            case TBool::False;
              variablesGlobales[var] = TBool::True;
              break;
            default:
              break;
          }
        }
      }
    }
    fin = std::chrono::high_resolution_clock::now();

    duracion = fin - inicio;
    
    cout << "Búsqueda local iterada" << endl;
    cout << "Insatisfechas: " << insatisfechas << endl;
    cout << "Tiempo de ejecución: " << duracion.count() << " s" << endl;
  }*/
  return 0;
}