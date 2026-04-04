#include "MaxSAT.h"
#include <algorithm>
#include <limits>
#include <random>

using namespace std;

// 1. FASE CONSTRUCTIVA DE GRASP
void Formula::construccionGRASP(vector<TBool> &vars, vector<Conteo> frecs, double alpha, mt19937 &gen)
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

        // Encontrar los valores máximo y mínimo de beneficio (apariciones)
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

        // Calcular el umbral para la RCL basado en alpha
        double umbral = s_max - alpha * (s_max - s_min);
        rcl.clear();

        // Llenar la RCL con los candidatos que superen o igualen el umbral
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

        // Elegir un candidato al azar de la RCL
        uniform_int_distribution<> dis(0, rcl.size() - 1);
        int idElegido = rcl[dis(gen)];

        // Asignar el valor de verdad que más le convenga
        bool valor = frecs[idElegido].pos >= frecs[idElegido].neg;
        vars[idElegido] = valor ? TBool::True : TBool::False;

        frecs[idElegido].pos = -99999;
        frecs[idElegido].neg = -99999;

        // Actualizar frecuencias de las variables restantes basadas en las cláusulas satisfechas
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

// 2. FASE PRINCIPAL DE GRASP (Iteración)
void Formula::busquedaGRASP(vector<TBool> &vars, int maxIteraciones, double alpha, mt19937 &gen, const vector<Conteo> &frecsOriginales)
{
    int mejorCostoGlobal = numeric_limits<int>::max();
    vector<TBool> mejorSolucionGlobal;

    for (int i = 0; i < maxIteraciones; i++)
    {
        vector<TBool> actual(vars.size(), TBool::Unknown);

        // Paso 1: Construcción codiciosa aleatorizada
        construccionGRASP(actual, frecsOriginales, alpha, gen);

        // Paso 2: Búsqueda Local para alcanzar el óptimo local (llama a tu módulo BL)
        busquedaLocal(actual);

        // Paso 3: Actualización de la mejor solución encontrada
        int costoFinal = calcularCosto(actual);
        if (costoFinal < mejorCostoGlobal)
        {
            mejorCostoGlobal = costoFinal;
            mejorSolucionGlobal = actual;
        }
    }

    // Retornamos la mejor solución después de todas las iteraciones
    vars = mejorSolucionGlobal;
}