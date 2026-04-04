#include "MaxSAT.h"
#include <algorithm>

using namespace std;

void Formula::busquedaLocalIterada(vector<TBool> &vars, int maxIteraciones, mt19937 &gen)
{
    int mejorCosto = calcularCosto(vars);
    vector<TBool> mejorSolucion = vars;

    // Distribución para elegir una variable al azar
    uniform_int_distribution<> dis(0, vars.size() - 1);

    for (int i = 0; i < maxIteraciones; i++)
    {
        vector<TBool> actual = mejorSolucion;

        // 1. Perturbación (Random k-flip del 5% de las variables)
        int k = max(1, (int)(vars.size() * 0.05));
        for (int j = 0; j < k; j++)
        {
            int idx = dis(gen);
            actual[idx] = (actual[idx] == TBool::True) ? TBool::False : TBool::True;
        }

        // 2. Búsqueda Local (Llamamos a tu otro módulo)
        busquedaLocal(actual);

        // 3. Criterio de Aceptación (Aceptamos solo si mejora estrictamente)
        int costoActual = calcularCosto(actual);
        if (costoActual < mejorCosto)
        {
            mejorCosto = costoActual;
            mejorSolucion = actual;
        }
    }

    // Guardamos la mejor solución encontrada en las variables originales
    vars = mejorSolucion;
}