#include "MaxSAT.h"
#include <algorithm>

using namespace std;

void Formula::solverConstructivo(vector<TBool> &variablesGlobales, vector<Conteo> frecs)
{
    int variablesPendientes = variablesGlobales.size();

    while (variablesPendientes > 0)
    {
        // Buscar la variable con mayor cantidad de apariciones (moda)
        auto moda = max_element(frecs.begin(), frecs.end(), [](const Conteo &a, const Conteo &b)
                                { return (a.pos + a.neg) < (b.pos + b.neg); });

        // Criterio de parada: si ya no hay variables con apariciones positivas/negativas
        if (moda->pos <= 0 && moda->neg <= 0)
            break;

        int idModa = distance(frecs.begin(), moda);
        bool valor = moda->pos >= moda->neg;

        // Asignar el valor a la variable global
        variablesGlobales[idModa] = valor ? TBool::True : TBool::False;

        // Marcar la variable como procesada
        moda->pos = -9999;
        moda->neg = -9999;

        // Actualizar el estado de las cláusulas donde aparece esta variable
        for (int idxClausula : apariciones[idModa])
        {
            Clausula &clausula = clausulas[idxClausula];
            if (clausula.getSatisfaccion() == TBool::Unknown)
            {
                clausula.setSatisfaccion(variablesGlobales, &frecs);
            }
        }
        variablesPendientes--;
    }

    // Rellenar las variables que no fueron asignadas (por no tener impacto restante) con False
    for (size_t i = 0; i < variablesGlobales.size(); i++)
    {
        if (variablesGlobales[i] == TBool::Unknown)
            variablesGlobales[i] = TBool::False;
    }
}