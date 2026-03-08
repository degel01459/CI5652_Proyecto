#include "MaxSAT.h"
#include <cmath>
#include <random>

using namespace std;

void Formula::recocidoSimulado(vector<TBool> &vars, mt19937 &gen, double tempInicial, double alpha, int iterPorTemp)
{
    int n = vars.size();
    vector<TBool> actual = vars;
    vector<TBool> mejorSolucionGlobal = vars;

    int costoActual = calcularCosto(actual);
    int mejorCostoGlobal = costoActual;

    double T = tempInicial;
    double T_min = 0.01;

    // Distribuciones para elegir variable aleatoria y calcular probabilidad
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
                // Si mejora el costo, lo aceptamos siempre
                actual[idx] = (actual[idx] == TBool::True) ? TBool::False : TBool::True;
                costoActual = nuevoCosto;

                // Actualizamos el óptimo global si es necesario
                if (costoActual < mejorCostoGlobal)
                {
                    mejorCostoGlobal = costoActual;
                    mejorSolucionGlobal = actual;
                }
            }
            else
            {
                // Si empeora, calculamos la probabilidad de aceptarlo (Criterio de Metropolis)
                double probabilidad = exp(-delta / T);
                if (probDist(gen) < probabilidad)
                {
                    actual[idx] = (actual[idx] == TBool::True) ? TBool::False : TBool::True;
                    costoActual = nuevoCosto;
                }
            }
        }
        // Enfriamiento del sistema
        T *= alpha;
    }

    // Retornamos la mejor solución histórica encontrada
    vars = mejorSolucionGlobal;
}