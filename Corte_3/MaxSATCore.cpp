// MaxSATCore.cpp
#include "MaxSAT.h"
#include <cmath>
#include <numeric>
#include <sstream>
#include <iomanip>

using namespace std;

// ======================================================================
// IMPLEMENTACIÓN DE LA CLASE CLAUSULA
// ======================================================================

Clausula::Clausula() { satisfaccion = TBool::Unknown; }

Clausula::Clausula(vector<int> vars)
{
    variables = vars;
    satisfaccion = TBool::Unknown;
}

const vector<int> &Clausula::getVariables() const { return variables; }
TBool Clausula::getSatisfaccion() const { return satisfaccion; }
void Clausula::reset() { satisfaccion = TBool::Unknown; }

void Clausula::setSatisfaccion(const vector<TBool> &variablesGlobales, vector<Conteo> *ptrFrecuencias)
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

bool Clausula::esSatisfecha(const vector<TBool> &variablesGlobales) const
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

bool Clausula::aparece(int variable)
{
    for (int v : variables)
        if (abs(v) == variable + 1)
            return true;
    return false;
}

void Clausula::actualizarFrecuencias(vector<Conteo> *ptrFrecuencias)
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

// ======================================================================
// IMPLEMENTACIÓN BÁSICA DE LA CLASE FORMULA
// ======================================================================

Formula::Formula(const vector<Clausula> &clauses, int numVars)
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

int Formula::calcularCosto(const vector<TBool> &vars)
{
    int costo = 0;
    for (const auto &c : clausulas)
    {
        if (!c.esSatisfecha(vars))
            costo++;
    }
    return costo;
}

int Formula::evaluarFlip(vector<TBool> &vars, int varIdx)
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

// ======================================================================
// FUNCIONES AUXILIARES (Lectura, Estadísticas y Formateo)
// ======================================================================

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