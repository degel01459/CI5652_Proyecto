#include "MaxSAT.h"

using namespace std;

void Formula::busquedaLocal(vector<TBool> &vars)
{
    bool mejora = true;
    while (mejora)
    {
        mejora = false;
        for (size_t i = 0; i < vars.size(); i++)
        {
            // Preguntamos: ¿Qué pasa si invierto la variable i?
            int delta = evaluarFlip(vars, i);

            if (delta < 0) // Si el delta es negativo, el costo baja (¡mejora!)
            {
                // Aplicamos el cambio permanentemente
                vars[i] = (vars[i] == TBool::True) ? TBool::False : TBool::True;
                mejora = true;
                break; // First-Improvement: aplicamos la primera mejora que encontramos y reiniciamos
            }
        }
    }
}