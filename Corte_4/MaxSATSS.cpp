#include "MaxSAT.h"
#include <vector>
#include <algorithm>
#include <chrono>
#include <random>

using namespace std;

// Implementación de Búsqueda Dispersa (Scatter Search)
void Formula::busquedaDispersa(vector<TBool> &mejorSolucion, int tamRefSet, double tiempoLimite, mt19937 &gen)
{
    auto start_time = chrono::high_resolution_clock::now();
    int numVariables = mejorSolucion.size();

    // 1. Inicializar el Conjunto de Referencia (RefSet)
    // El RefSet guardará pares de <Costo, Solucion>
    vector<pair<int, vector<TBool>>> refSet;

    // Evaluamos la solución inicial que nos pasaron por parámetro
    int costoMejor = calcularCosto(mejorSolucion);
    refSet.push_back({costoMejor, mejorSolucion});

    // Generamos el resto del RefSet con variaciones aleatorias de la solución inicial
    // para asegurar "Diversidad" (un pilar del Scatter Search)
    uniform_int_distribution<> disBool(0, 1);
    uniform_int_distribution<> disVar(0, numVariables - 1);

    while (refSet.size() < (size_t)tamRefSet)
    {
        vector<TBool> nuevaSol = mejorSolucion;
        // Mutamos aleatoriamente un 20% de las variables para crear diversidad
        int numMutaciones = max(1, numVariables / 5);
        for (int i = 0; i < numMutaciones; i++)
        {
            int varIdx = disVar(gen);
            // Invertir el valor de la variable
            // Asumimos que TBool tiene estados equivalentes a true/false o 1/0
            if (nuevaSol[varIdx] == mejorSolucion[varIdx])
            {
                nuevaSol[varIdx] = disBool(gen) ? mejorSolucion[varIdx] : (TBool)(1 - (int)mejorSolucion[varIdx]);
            }
        }

        int costoNueva = calcularCosto(nuevaSol);
        refSet.push_back({costoNueva, nuevaSol});
    }

    // Ordenamos el RefSet de menor a mayor costo (asumiendo que menor costo es mejor)
    // Si en tu problema "mayor costo es mejor", cambia el '<' por '>'
    sort(refSet.begin(), refSet.end(), [](const pair<int, vector<TBool>> &a, const pair<int, vector<TBool>> &b)
         { return a.first < b.first; });

    // 2. Ciclo Principal del Scatter Search
    bool mejora = true;
    while (mejora)
    {
        mejora = false;
        vector<pair<int, vector<TBool>>> nuevosCandidatos;

        // Comprobamos el tiempo límite
        auto current_time = chrono::high_resolution_clock::now();
        double elapsed = chrono::duration<double>(current_time - start_time).count();
        if (elapsed >= tiempoLimite)
            break;

        // Combinación de Subconjuntos (Pares del RefSet)
        for (int i = 0; i < tamRefSet; i++)
        {
            for (int j = i + 1; j < tamRefSet; j++)
            {

                vector<TBool> hijo = refSet[i].second;

                // Operador de Combinación simple (Crossover Uniforme)
                for (int k = 0; k < numVariables; k++)
                {
                    if (disBool(gen))
                    {
                        hijo[k] = refSet[j].second[k];
                    }
                }

                // Aquí idealmente aplicarías un Método de Mejora (Búsqueda Local) al 'hijo'
                // busquedaLocal(hijo); // Descomenta si tienes esta función disponible y quieres más potencia

                int costoHijo = calcularCosto(hijo);
                nuevosCandidatos.push_back({costoHijo, hijo});
            }
        }

        // 3. Actualización del RefSet
        for (auto &candidato : nuevosCandidatos)
        {
            // Si el candidato es mejor que el peor del RefSet
            if (candidato.first < refSet.back().first)
            {
                // Verificamos que no sea idéntico a uno que ya está para mantener diversidad
                bool duplicado = false;
                for (const auto &ref : refSet)
                {
                    if (ref.first == candidato.first)
                    { // Simplificación rápida de igualdad
                        duplicado = true;
                        break;
                    }
                }

                if (!duplicado)
                {
                    refSet.back() = candidato; // Reemplazamos al peor
                    // Reordenamos
                    sort(refSet.begin(), refSet.end(), [](const pair<int, vector<TBool>> &a, const pair<int, vector<TBool>> &b)
                         { return a.first < b.first; });
                    mejora = true; // Hubo un cambio, seguimos iterando
                }
            }
        }
    }

    // 4. Retornar la mejor solución encontrada
    mejorSolucion = refSet[0].second;
}