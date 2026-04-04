#include "MaxSAT.h"
#include <algorithm>
#include <random>

using namespace std;

void Formula::algoritmoGenetico(vector<TBool> &vars, int tamanoPoblacion, int maxGeneraciones, double probCruce, double probMutacion, mt19937 &gen)
{
    int n = vars.size();

    // Matriz para guardar la población actual
    vector<vector<TBool>> poblacion(tamanoPoblacion, vector<TBool>(n));
    uniform_int_distribution<> disBool(0, 1);

    // Inicialización de la población: el individuo 0 es la solución dada (probablemente heurística)
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

    // Lambda: Selección por Torneo (k individuos pelean y sobrevive el que tiene menor costo)
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

    // Lambda: Mutación (invierte un gen/variable con una probabilidad baja)
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

    // Bucle principal evolutivo
    for (int genIdx = 0; genIdx < maxGeneraciones; genIdx++)
    {
        // 1. Evaluar población actual
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

        // 2. Elitismo: El mejor individuo pasa intacto a la siguiente generación
        int idxMejorActual = distance(costos.begin(), min_element(costos.begin(), costos.end()));
        nuevaPoblacion[0] = poblacion[idxMejorActual];

        // 3. Selección, Cruce y Mutación (se generan de 2 en 2)
        for (int i = 1; i < tamanoPoblacion; i += 2)
        {
            int p1 = seleccionTorneo(costos, 3);
            int p2 = seleccionTorneo(costos, 3);

            vector<TBool> hijo1 = poblacion[p1];
            vector<TBool> hijo2 = poblacion[p2];

            // Operador de cruce uniforme
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

            // Aplicar mutación a los hijos
            mutar(hijo1);
            mutar(hijo2);

            // Agregar hijos a la nueva generación
            nuevaPoblacion[i] = hijo1;
            if (i + 1 < tamanoPoblacion)
            {
                nuevaPoblacion[i + 1] = hijo2;
            }
        }

        // Reemplazo poblacional
        poblacion = move(nuevaPoblacion);
    }

    // Devolvemos el mejor global de todas las generaciones
    vars = mejorSolucionGlobal;
}