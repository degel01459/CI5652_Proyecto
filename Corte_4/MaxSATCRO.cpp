#include "MaxSAT.h"
#include <algorithm>
#include <chrono>
#include <random>

using namespace std;

// Representación química de una solución
struct Molecula {
    vector<TBool> vars;
    int PE;        // Energía Potencial (Costo)
    double KE;     // Energía Cinética (Temperatura local)
    int numHits;   // Número de colisiones (Estancamiento)
};

void Formula::optimizacionReaccionQuimica(vector<TBool> &varsGlobal, int popSizeInicial, double tiempoLimite, mt19937 &gen)
{
    int n = varsGlobal.size();
    auto tiempoInicio = chrono::high_resolution_clock::now();

    // Parámetros termodinámicos de CRO
    double KE_inicial = 1000.0;   
    double bufferEnergia = 0.0;   
    double KELossRate = 0.2;      
    double moleCollRatio = 0.2;   
    int maxHits = 10;             

    // Inicializar población
    vector<Molecula> poblacion;
    poblacion.reserve(popSizeInicial * 2); // Pre-reservar memoria para evitar reasignaciones
    uniform_int_distribution<> disBool(0, 1);
    
    int mejorPEGlobal = calcularCosto(varsGlobal);
    vector<TBool> mejorSolucionGlobal = varsGlobal;

    for (int i = 0; i < popSizeInicial; i++)
    {
        Molecula m;
        m.vars = varsGlobal;
        for (int j = 0; j < n; j++) {
            if (disBool(gen)) m.vars[j] = (m.vars[j] == TBool::True) ? TBool::False : TBool::True;
        }
        m.PE = calcularCosto(m.vars);
        m.KE = KE_inicial;
        m.numHits = 0;
        poblacion.push_back(m);

        if (m.PE < mejorPEGlobal) {
            mejorPEGlobal = m.PE;
            mejorSolucionGlobal = m.vars;
        }
    }

    uniform_real_distribution<> probDist(0.0, 1.0);
    uniform_int_distribution<> varDist(0, n - 1);
    uniform_int_distribution<> molDist(0, 1);
    unsigned long long iteraciones = 0; // Optimización 3: Contador para no llamar al reloj siempre

    while (true)
    {
        // Control Termodinámico del tiempo (Solo consultamos al SO cada 100 iteraciones)
        if (iteraciones++ % 10000 == 0) {
            auto tiempoActual = chrono::high_resolution_clock::now();
            chrono::duration<double> tiempoTranscurrido = tiempoActual - tiempoInicio;
            if (tiempoTranscurrido.count() > tiempoLimite) break;
        }

        bool colisionBimolecular = (probDist(gen) < moleCollRatio) && (poblacion.size() > 1);

        if (!colisionBimolecular) 
        {
            // --- REACCIONES UNIMOLECULARES ---
            molDist.param(uniform_int_distribution<>::param_type(0, poblacion.size() - 1));
            int idx = molDist(gen);
            Molecula &m = poblacion[idx];

            if (m.numHits > maxHits) 
            {
                int flip1 = varDist(gen);
                int flip2 = varDist(gen);

                int nuevoPE1 = m.PE + evaluarFlip(m.vars, flip1);
                int nuevoPE2 = m.PE + evaluarFlip(m.vars, flip2);

                if (m.PE + m.KE + bufferEnergia >= nuevoPE1 + nuevoPE2) {
                    double energiaRestante = (m.PE + m.KE + bufferEnergia) - (nuevoPE1 + nuevoPE2);
                    
                    Molecula m2 = m;
                    m2.vars[flip2] = (m2.vars[flip2] == TBool::True) ? TBool::False : TBool::True;
                    m2.PE = nuevoPE2;
                    m2.KE = energiaRestante / 2.0;
                    m2.numHits = 0;
                    
                    m.vars[flip1] = (m.vars[flip1] == TBool::True) ? TBool::False : TBool::True;
                    m.PE = nuevoPE1;
                    m.KE = energiaRestante / 2.0;
                    m.numHits = 0;
                    
                    // --- CORRECCIÓN VITAL AQUÍ ---
                    // Actualizamos récords MIENTRAS 'm' SIGUE A SALVO
                    if (nuevoPE1 < mejorPEGlobal) { mejorPEGlobal = nuevoPE1; mejorSolucionGlobal = m.vars; }
                    if (nuevoPE2 < mejorPEGlobal) { mejorPEGlobal = nuevoPE2; mejorSolucionGlobal = m2.vars; }
                    
                    // Ahora sí añadimos m2. Si el vector colapsa y mueve la memoria,
                    // ya no importa porque no volveremos a usar 'm' en este bloque.
                    poblacion.push_back(std::move(m2));
                    bufferEnergia = 0.0;
                    
                } else {
                    m.numHits++;
                }
            } 
            else
            {
                // 2. COLISIÓN CON LA PARED (Se mantiene igual de rápida)
                int flipIdx = varDist(gen);
                int deltaPE = evaluarFlip(m.vars, flipIdx); 
                int nuevoPE = m.PE + deltaPE;

                if (m.PE + m.KE >= nuevoPE) {
                    double KE_nuevo = (m.PE - nuevoPE + m.KE);
                    double disipado = KE_nuevo * KELossRate;
                    
                    m.KE = KE_nuevo - disipado;
                    bufferEnergia += disipado;
                    
                    m.vars[flipIdx] = (m.vars[flipIdx] == TBool::True) ? TBool::False : TBool::True;
                    m.PE = nuevoPE;
                    m.numHits++;

                    if (m.PE < mejorPEGlobal) {
                        mejorPEGlobal = m.PE;
                        mejorSolucionGlobal = m.vars;
                    }
                } else {
                    m.numHits++;
                }
            }
        } 
        else 
        {
            // --- REACCIONES BIMOLECULARES ---
            molDist.param(uniform_int_distribution<>::param_type(0, poblacion.size() - 1));
            int idx1 = molDist(gen);
            int idx2 = molDist(gen);
            while (idx1 == idx2) idx2 = molDist(gen);

            Molecula &m1 = poblacion[idx1];
            Molecula &m2 = poblacion[idx2];

            if (m1.KE < 10.0 && m2.KE < 10.0) {
                // 3. SÍNTESIS OPTIMIZADA
                Molecula nuevaMol;
                nuevaMol.vars = m1.vars; 
                
                uint32_t bitsAleatorios = 0; // Guardará 32 decisiones aleatorias de golpe

                for (int j = 0; j < n; j++) {
                    // Solo le pedimos un nuevo número al generador cada 32 variables
                    if (j % 32 == 0) {
                        bitsAleatorios = gen(); 
                    }

                    // Usamos el último bit (Bitwise AND). Si es 1, copiamos el gen del padre 2.
                    if (bitsAleatorios & 1) {
                        nuevaMol.vars[j] = m2.vars[j];
                    }

                    // Desplazamos a la derecha para usar el siguiente bit en la próxima iteración
                    bitsAleatorios >>= 1; 
                }
                
                nuevaMol.PE = calcularCosto(nuevaMol.vars); // Aquí sí es obligatorio
                
                if (m1.PE + m2.PE + m1.KE + m2.KE >= nuevaMol.PE) {
                    nuevaMol.KE = (m1.PE + m2.PE + m1.KE + m2.KE) - nuevaMol.PE;
                    nuevaMol.numHits = 0;
                    
                    // 1. Evaluamos los récords MIENTRAS nuevaMol ESTÁ INTACTA
                    if (nuevaMol.PE < mejorPEGlobal) {
                        mejorPEGlobal = nuevaMol.PE;
                        mejorSolucionGlobal = nuevaMol.vars;
                    }
                    
                    // 2. AHORA SÍ le robamos la memoria (Swap and Pop)
                    if (idx2 == poblacion.size() - 1) {
                        poblacion[idx1] = std::move(nuevaMol); // idx2 se borra abajo
                    } else if (idx1 == poblacion.size() - 1) {
                        poblacion[idx2] = std::move(nuevaMol); // idx1 se borra abajo
                    } else {
                        poblacion[idx1] = std::move(nuevaMol);
                        poblacion[idx2] = std::move(poblacion.back()); // Llenamos el hueco de idx2
                    }
                    poblacion.pop_back(); // Siempre borramos el último

                } else {
                    m1.numHits++; m2.numHits++;
                }
            }
            else
            {
                // 4. COLISIÓN INTERMOLECULAR OPTIMIZADA
                int f1 = varDist(gen);
                int f2 = varDist(gen);
                
                // Evaluamos sin clonar memoria
                int nuevoPE1 = m1.PE + evaluarFlip(m1.vars, f1);
                int nuevoPE2 = m2.PE + evaluarFlip(m2.vars, f2);

                if (m1.PE + m2.PE + m1.KE + m2.KE >= nuevoPE1 + nuevoPE2) {
                    double energiaRestante = (m1.PE + m2.PE + m1.KE + m2.KE) - (nuevoPE1 + nuevoPE2);
                    
                    // Aplicamos cambios en el lugar
                    m1.vars[f1] = (m1.vars[f1] == TBool::True) ? TBool::False : TBool::True;
                    m2.vars[f2] = (m2.vars[f2] == TBool::True) ? TBool::False : TBool::True;
                    
                    m1.PE = nuevoPE1;
                    m2.PE = nuevoPE2;
                    m1.KE = energiaRestante / 2.0;
                    m2.KE = energiaRestante / 2.0;
                    m1.numHits++; m2.numHits++;

                    if (m1.PE < mejorPEGlobal) { mejorPEGlobal = m1.PE; mejorSolucionGlobal = m1.vars; }
                    if (m2.PE < mejorPEGlobal) { mejorPEGlobal = m2.PE; mejorSolucionGlobal = m2.vars; }
                } else {
                    m1.numHits++; m2.numHits++;
                }
            }
        }
    }

    varsGlobal = mejorSolucionGlobal;
}