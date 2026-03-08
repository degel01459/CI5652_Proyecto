#include "MaxSAT.h"
#include <random>

using namespace std;

void Formula::busquedaTabu(vector<TBool> &vars, int maxIteraciones, int tenureBase)
{
    int n = vars.size();

    // Lista Tabú: Guarda hasta qué iteración una variable está "bloqueada"
    vector<int> tabuUntil(n, 0);

    vector<TBool> mejorSolucionGlobal = vars;
    int mejorCostoGlobal = calcularCosto(vars);
    int costoActual = mejorCostoGlobal;

    // Generador local para el tenure dinámico
    random_device rd;
    mt19937 gen(rd());
    uniform_int_distribution<> disTenure(0, 5);

    for (int iter = 1; iter <= maxIteraciones; iter++)
    {
        int mejorVarIdx = -1;
        int mejorDelta = 1e9; // Inicializamos con un valor muy grande

        // Explorar todo el vecindario (todos los flips posibles)
        for (int i = 0; i < n; i++)
        {
            int delta = evaluarFlip(vars, i);
            int nuevoCosto = costoActual + delta;

            bool esTabu = (iter < tabuUntil[i]);
            bool aspira = (nuevoCosto < mejorCostoGlobal); // Criterio de aspiración

            // Consideramos la variable si NO es tabú, o si siendo tabú mejora el óptimo global
            if (!esTabu || aspira)
            {
                if (delta < mejorDelta)
                {
                    mejorDelta = delta;
                    mejorVarIdx = i;
                }
            }
        }

        // Aplicar el mejor movimiento encontrado (incluso si empeora la solución actual)
        if (mejorVarIdx != -1)
        {
            vars[mejorVarIdx] = (vars[mejorVarIdx] == TBool::True) ? TBool::False : TBool::True;
            costoActual += mejorDelta;

            // Bloquear la variable que acabamos de tocar (Tenure dinámico)
            tabuUntil[mejorVarIdx] = iter + tenureBase + disTenure(gen);

            // Actualizar el óptimo global si es necesario
            if (costoActual < mejorCostoGlobal)
            {
                mejorCostoGlobal = costoActual;
                mejorSolucionGlobal = vars;
            }
        }
    }

    // Guardamos la mejor solución histórica encontrada
    vars = mejorSolucionGlobal;
}