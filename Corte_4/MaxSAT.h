// MaxSAT.h
#pragma once
#include <vector>
#include <random>
#include <string>

using namespace std;

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
    Clausula();
    Clausula(vector<int> vars);
    const vector<int> &getVariables() const;
    TBool getSatisfaccion() const;
    void reset();
    void setSatisfaccion(const vector<TBool> &variablesGlobales, vector<Conteo> *ptrFrecuencias);
    bool esSatisfecha(const vector<TBool> &variablesGlobales) const;
    bool aparece(int variable);
    void actualizarFrecuencias(vector<Conteo> *ptrFrecuencias);
};

class Formula
{
private:
    vector<Clausula> clausulas;
    vector<vector<int>> apariciones;

public:
    Formula() {}
    Formula(const vector<Clausula> &clauses, int numVars);

    int calcularCosto(const vector<TBool> &vars);
    int evaluarFlip(vector<TBool> &vars, int varIdx);

    // Declaración de todas tus metaheurísticas:
    void solverConstructivo(vector<TBool> &variablesGlobales, vector<Conteo> frecs);
    void busquedaLocal(vector<TBool> &vars);
    void busquedaLocalIterada(vector<TBool> &vars, int maxIteraciones, mt19937 &gen);
    void busquedaTabu(vector<TBool> &vars, int maxIteraciones, int tenureBase);
    void recocidoSimulado(vector<TBool> &vars, mt19937 &gen, double tempInicial = 10.0, double alpha = 0.95, int iterPorTemp = 100);
    void construccionGRASP(vector<TBool> &vars, vector<Conteo> frecs, double alpha, mt19937 &gen);
    void busquedaGRASP(vector<TBool> &vars, int maxIteraciones, double alpha, mt19937 &gen, const vector<Conteo> &frecsOriginales);
    void algoritmoGenetico(vector<TBool> &vars, int tamanoPoblacion, int maxGeneraciones, double probCruce, double probMutacion, mt19937 &gen);
    void algoritmoMemetico(vector<TBool> &vars, int tamPoblacion, int maxGeneraciones, double tiempoLimiteSegundos, int maxGeneracionesSinMejora, mt19937 &gen);
    void busquedaDispersa(vector<TBool> &vars, int tamRefSet, double tiempoLimiteSegundos, mt19937 &gen);
    void optimizacionColoniaHormigas(vector<TBool> &vars, int numHormigas, double tiempoLimiteSegundos, mt19937 &gen);
    void optimizacionReaccionQuimica(vector<TBool> &mejorSolucionGlobal, int popSizeInicial, double tiempoLimite, mt19937 &gen);
};

// Funciones Auxiliares Globales
Clausula crearClausula(string linea, vector<Conteo> &frecuencias);
pair<int, int> leerPreambulo(string linea);
double promedio(const vector<double> &v);
double desviacionEstandar(const vector<double> &v, const double media);
string formatearMedida(double media, double desv_est);