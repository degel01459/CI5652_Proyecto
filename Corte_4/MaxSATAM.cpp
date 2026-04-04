#include "MaxSAT.h"
#include <algorithm>
#include <chrono>
#include <random>

using namespace std;

void Formula::algoritmoMemetico(vector<TBool> &vars, int tamPoblacion, int maxGeneraciones, double tiempoLimiteSegundos, int maxGeneracionesSinMejora, mt19937 &gen)
{
    int n = vars.size();
    auto tiempoInicio = chrono::high_resolution_clock::now();

    // 1. Inicialización de la Población
    vector<vector<TBool>> poblacion(tamPoblacion, vector<TBool>(n));
    uniform_int_distribution<> disBool(0, 1);

    // El individuo 0 es nuestra solución inicial (ej. Heurística)
    poblacion[0] = vars;
    busquedaLocal(poblacion[0]); // Optimizamos el individuo inicial con BL

    for (int i = 1; i < tamPoblacion; i++)
    {
        for (int j = 0; j < n; j++)
        {
            poblacion[i][j] = disBool(gen) ? TBool::True : TBool::False;
        }
    }

    vector<TBool> mejorSolucionGlobal = poblacion[0];
    int mejorCostoGlobal = calcularCosto(mejorSolucionGlobal);
    int genSinMejora = 0;

    uniform_real_distribution<> disProb(0.0, 1.0);
    uniform_int_distribution<> disTorneo(0, tamPoblacion - 1);

    // Parámetros internos heredados (puedes ajustarlos si lo necesitas)
    double probCruce = 0.85;
    double probMutacion = 1.0 / n;

    // Lambda: Selección por Torneo
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

    // Lambda: Mutación Simple
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

    // 2. Ciclo Evolutivo Memético
    for (int genIdx = 0; genIdx < maxGeneraciones; genIdx++)
    {
        // Control de tiempo máximo de ejecución
        auto tiempoActual = chrono::high_resolution_clock::now();
        chrono::duration<double> tiempoTranscurrido = tiempoActual - tiempoInicio;
        if (tiempoTranscurrido.count() > tiempoLimiteSegundos)
        {
            break;
        }

        // Control de estancamiento (Early Stopping)
        if (genSinMejora >= maxGeneracionesSinMejora)
        {
            break;
        }

        vector<int> costos(tamPoblacion);
        for (int i = 0; i < tamPoblacion; i++)
        {
            costos[i] = calcularCosto(poblacion[i]);
        }

        vector<vector<TBool>> nuevaPoblacion(tamPoblacion, vector<TBool>(n));

        // Elitismo: Pasar al mejor a la siguiente generación
        int idxMejorActual = distance(costos.begin(), min_element(costos.begin(), costos.end()));
        nuevaPoblacion[0] = poblacion[idxMejorActual];

        bool mejoraEnGeneracion = false;

        // Reproducción
        for (int i = 1; i < tamPoblacion; i += 2)
        {
            int p1 = seleccionTorneo(costos, 3);
            int p2 = seleccionTorneo(costos, 3);

            vector<TBool> hijo1 = poblacion[p1];
            vector<TBool> hijo2 = poblacion[p2];

            // Cruce Uniforme
            if (disProb(gen) < probCruce)
            {
                for (int j = 0; j < n; j++)
                {
                    if (disProb(gen) < 0.5)
                        swap(hijo1[j], hijo2[j]);
                }
            }

            // Mutación
            mutar(hijo1);
            mutar(hijo2);

            // ------------------------------------------------------------------
            // EL FACTOR MEMÉTICO: Búsqueda Local intensiva a los nuevos hijos
            // ------------------------------------------------------------------
            busquedaLocal(hijo1);
            busquedaLocal(hijo2);

            // Evaluar los nuevos hijos para registrar mejoras globales inmediatas
            int costoH1 = calcularCosto(hijo1);
            if (costoH1 < mejorCostoGlobal)
            {
                mejorCostoGlobal = costoH1;
                mejorSolucionGlobal = hijo1;
                mejoraEnGeneracion = true;
            }

            int costoH2 = calcularCosto(hijo2);
            if (costoH2 < mejorCostoGlobal)
            {
                mejorCostoGlobal = costoH2;
                mejorSolucionGlobal = hijo2;
                mejoraEnGeneracion = true;
            }

            nuevaPoblacion[i] = hijo1;
            if (i + 1 < tamPoblacion)
                nuevaPoblacion[i + 1] = hijo2;
        }

        poblacion = move(nuevaPoblacion);

        // Actualizar contador de estancamiento
        if (mejoraEnGeneracion)
        {
            genSinMejora = 0;
        }
        else
        {
            genSinMejora++;
        }
    }

    // Retornamos la mejor solución encontrada tras las evoluciones y mejoras locales
    vars = mejorSolucionGlobal;
}