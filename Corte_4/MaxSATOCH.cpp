#include "MaxSAT.h"
#include <algorithm>
#include <chrono>
#include <random>

using namespace std;

void Formula::optimizacionColoniaHormigas(vector<TBool> &vars, int numHormigas, double tiempoLimiteSegundos, mt19937 &gen)
{
    int n = vars.size();
    auto tiempoInicio = chrono::high_resolution_clock::now();

    // 1. Inicialización de feromonas
    // feromonas[i] indica la "conveniencia" de asignar True a la variable i.
    // Iniciamos en 0.5 (neutralidad total).
    vector<double> feromonas(n, 0.5);

    // Parámetros de ACO
    double evaporacion = 0.1; // Tasa de evaporación (rho)
    double q0 = 0.8;          // Probabilidad de explotación vs exploración

    vector<TBool> mejorSolucionGlobal = vars;
    int mejorCostoGlobal = calcularCosto(mejorSolucionGlobal);

    uniform_real_distribution<> disProb(0.0, 1.0);

    // Bucle principal de generaciones de hormigas
    while (true)
    {
        // Control de tiempo máximo de ejecución
        auto tiempoActual = chrono::high_resolution_clock::now();
        chrono::duration<double> tiempoTranscurrido = tiempoActual - tiempoInicio;
        if (tiempoTranscurrido.count() > tiempoLimiteSegundos)
        {
            break;
        }

        vector<TBool> mejorSolucionIteracion;
        int mejorCostoIteracion = 1e9;

        // 2. Construcción de soluciones por cada hormiga
        for (int h = 0; h < numHormigas; h++)
        {
            vector<TBool> solucionHormiga(n, TBool::Unknown);

            for (int i = 0; i < n; i++)
            {
                // Transición de estado pseudo-aleatoria proporcional
                if (disProb(gen) < q0)
                {
                    // Explotación: elegir determinísticamente lo que dicta la feromona
                    solucionHormiga[i] = (feromonas[i] >= 0.5) ? TBool::True : TBool::False;
                }
                else
                {
                    // Exploración: ruleta probabilística basada en el nivel de feromona
                    solucionHormiga[i] = (disProb(gen) < feromonas[i]) ? TBool::True : TBool::False;
                }
            }

            // Opcional pero recomendado: Aplicar Búsqueda Local a la hormiga
            // Esto convierte el algoritmo en un ACO Híbrido, muchísimo más potente
            busquedaLocal(solucionHormiga);

            int costoHormiga = calcularCosto(solucionHormiga);

            // Guardar la mejor hormiga de esta iteración
            if (costoHormiga < mejorCostoIteracion)
            {
                mejorCostoIteracion = costoHormiga;
                mejorSolucionIteracion = solucionHormiga;
            }
        }

        // Actualizar el óptimo global histórico
        if (mejorCostoIteracion < mejorCostoGlobal)
        {
            mejorCostoGlobal = mejorCostoIteracion;
            mejorSolucionGlobal = mejorSolucionIteracion;
        }

        // 3. Actualización de Feromonas (Evaporación + Aporte)
        for (int i = 0; i < n; i++)
        {
            // Evaporación constante
            feromonas[i] = (1.0 - evaporacion) * feromonas[i];

            // Aporte de feromonas (Depositado por la mejor hormiga global)
            // Si en la mejor solución la variable es True, reforzamos hacia 1.0
            if (mejorSolucionGlobal[i] == TBool::True)
            {
                feromonas[i] += evaporacion;
            }

            // Límites de feromonas al estilo MAX-MIN Ant System (MMAS)
            // Mantiene las feromonas entre [0.01 y 0.99] para nunca perder del todo la exploración
            feromonas[i] = max(0.01, min(0.99, feromonas[i]));
        }
    }

    // Retornamos la mejor solución de la colonia
    vars = mejorSolucionGlobal;
}