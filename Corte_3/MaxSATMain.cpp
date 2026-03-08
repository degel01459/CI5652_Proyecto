#include "MaxSAT.h"
#include <iostream>
#include <iomanip>
#include <fstream>
#include <chrono>
#include <algorithm>
#include <omp.h>
#include <string>
#include <cmath>

using namespace std;

int main(int argc, char const *argv[])
{
    // Optimizacion de I/O
    ios_base::sync_with_stdio(false);
    cin.tie(NULL);

    if (argc < 2)
    {
        cout << "Uso Terminal: ./solver <repeticiones> <algoritmo> archivo1.cnf [archivo2.cnf ...]" << endl;
        cout << "Uso VS Code : ./solver archivo1.cnf [archivo2.cnf ...] (Usa valores por defecto)" << endl;
        cout << "Algoritmos disponibles: all, heuristico, ls, ils, tabu, sa, grasp, ag, memetico, aco, ss" << endl;
        return 1;
    }

    int numCorridas = 1;
    string algoritmo = "all";
    int startFileIdx = 1;

    // 1. Detección Inteligente de Argumentos
    string primerArg = argv[1];
    bool esNumero = !primerArg.empty() && std::all_of(primerArg.begin(), primerArg.end(), ::isdigit);

    if (esNumero)
    {
        if (argc < 4)
        {
            cout << "[ERROR] Si indicas repeticiones, debes indicar el algoritmo y al menos un archivo." << endl;
            return 1;
        }
        numCorridas = stoi(primerArg);
        algoritmo = argv[2];
        startFileIdx = 3;
        transform(algoritmo.begin(), algoritmo.end(), algoritmo.begin(), ::tolower);
    }
    else
    {
        cout << "[INFO] Modo VS Code detectado. Usando " << numCorridas << " repeticiones y algoritmo '" << algoritmo << "' por defecto." << endl;
    }

    // 2. Interruptores de ejecución
    bool runAll = (algoritmo.find("all") != string::npos);
    bool runH = runAll || (algoritmo.find("heuristico") != string::npos);
    bool runLS = runAll || (algoritmo.find("ls") != string::npos);
    bool runILS = runAll || (algoritmo.find("ils") != string::npos);
    bool runTS = runAll || (algoritmo.find("tabu") != string::npos);
    bool runSA = runAll || (algoritmo.find("sa") != string::npos);
    bool runGRASP = runAll || (algoritmo.find("grasp") != string::npos);
    bool runAG = runAll || (algoritmo.find("ag") != string::npos);
    bool runAM = runAll || (algoritmo.find("memetico") != string::npos);
    bool runACO = runAll || (algoritmo.find("aco") != string::npos);
    bool runSS = runAll || (algoritmo.find("ss") != string::npos);

    bool necesitaHeuristica = runH || runLS || runILS || runTS || runSA || runAG || runAM || runSS || runACO;

    // --- CONFIGURACIÓN DEL CSV ---
    string nombreArchivoLimpio = algoritmo;
    replace(nombreArchivoLimpio.begin(), nombreArchivoLimpio.end(), ',', '_');
    replace(nombreArchivoLimpio.begin(), nombreArchivoLimpio.end(), ' ', '_');

    string nombreCSV = "resultados_" + nombreArchivoLimpio + "_" + to_string(numCorridas) + "rep.csv";

    ofstream archivoCSVCabecera(nombreCSV);
    if (archivoCSVCabecera.is_open())
    {
        // Se quitó "Casos"
        archivoCSVCabecera << "Archivo,Exacto";
        if (runH)
            archivoCSVCabecera << ",Costo H,T. H(s)";
        if (runLS)
            archivoCSVCabecera << ",Costo LS,T. LS(s)";
        if (runILS)
            archivoCSVCabecera << ",Costo ILS,T. ILS(s)";
        if (runTS)
            archivoCSVCabecera << ",Costo TS,T. TS(s)";
        if (runSA)
            archivoCSVCabecera << ",Costo SA,T. SA(s)";
        if (runGRASP)
            archivoCSVCabecera << ",Costo GRASP,T. GRASP(s)";
        if (runAG)
            archivoCSVCabecera << ",Costo AG,T. AG(s)";
        if (runAM)
            archivoCSVCabecera << ",Costo AM,T. AM(s)";
        if (runSS)
            archivoCSVCabecera << ",Costo SS,T. SS(s)";
        if (runACO)
            archivoCSVCabecera << ",Costo ACO,T. ACO(s)";
        archivoCSVCabecera << "\n";
        archivoCSVCabecera.close();
    }

    // 3. Cabecera Dinámica (Consola)
    int wC = 11;
    int wT = 13;

    // Ajustamos la base de la línea (quitando los ~11 caracteres que ocupaba Casos)
    int lineLength = 29 + (runH + runLS + runILS + runTS + runSA + runGRASP + runAG + runAM + runSS + runACO) * (wC + wT + 6);

    cout << string(lineLength, '=') << "\n";
    cout << " REPORTE COMPARATIVO | Repeticiones: " << numCorridas << " | Algoritmo Evaluado: " << algoritmo << "\n";
    cout << " CSV Generado: " << nombreCSV << "\n";
    cout << string(lineLength, '=') << "\n";

    // Se quitó "Casos"
    cout << left << setw(18) << "Archivo" << " | " << setw(8) << "Exacto";
    if (runH)
        cout << " | " << setw(wC) << "Costo H" << " | " << setw(wT) << "T. H(s)";
    if (runLS)
        cout << " | " << setw(wC) << "Costo LS" << " | " << setw(wT) << "T. LS(s)";
    if (runILS)
        cout << " | " << setw(wC) << "Costo ILS" << " | " << setw(wT) << "T. ILS(s)";
    if (runTS)
        cout << " | " << setw(wC) << "Costo TS" << " | " << setw(wT) << "T. TS(s)";
    if (runSA)
        cout << " | " << setw(wC) << "Costo SA" << " | " << setw(wT) << "T. SA(s)";
    if (runGRASP)
        cout << " | " << setw(wC) << "Costo GRASP" << " | " << setw(wT) << "T. GRASP(s)";
    if (runAG)
        cout << " | " << setw(wC) << "Costo AG" << " | " << setw(wT) << "T. AG(s)";
    if (runAM)
        cout << " | " << setw(wC) << "Costo AM" << " | " << setw(wT) << "T. AM(s)";
    if (runSS)
        cout << " | " << setw(wC) << "Costo SS" << " | " << setw(wT) << "T. SS(s)";
    if (runACO)
        cout << " | " << setw(wC) << "Costo ACO" << " | " << setw(wT) << "T. ACO(s)";
    cout << "\n"
         << string(lineLength, '-') << "\n";

// 4. Ejecución Condicionada
#pragma omp parallel for schedule(dynamic)
    for (int f = startFileIdx; f < argc; f++)
    {
        string nombreArchivo = argv[f];

        size_t seed = hash<string>{}(nombreArchivo) + omp_get_thread_num();
        mt19937 gen(seed);

        ifstream archivo(nombreArchivo);
        if (!archivo.is_open())
            continue;

        string linea;
        pair<int, int> datosFormula = {0, 0};
        vector<Clausula> clausulasBase;
        vector<Conteo> frecuenciasBase;

        while (getline(archivo, linea))
        {
            size_t inicio = linea.find_first_not_of(" \t\r\n");
            if (inicio == string::npos)
                continue;

            char primerCaracter = linea[inicio];
            if (primerCaracter == 'c')
                continue;
            else if (primerCaracter == 'p')
            {
                datosFormula = leerPreambulo(linea);
                frecuenciasBase.resize(datosFormula.first);
                clausulasBase.reserve(datosFormula.second);
            }
            else if (isdigit(primerCaracter) || primerCaracter == '-')
            {
                clausulasBase.push_back(crearClausula(linea, frecuenciasBase));
            }
        }
        archivo.close();

        vector<double> tH, tLS, tILS, tTS, tSA, tGRASP, tAG, tAM, tSS, tACO;
        vector<double> cH, cLS, cILS, cTS, cSA, cGRASP, cAG, cAM, cSS, cACO;

        for (int iter = 0; iter < numCorridas; iter++)
        {
            Formula problema(clausulasBase, datosFormula.first);
            vector<TBool> varsInicial(datosFormula.first, TBool::Unknown);

            if (necesitaHeuristica)
            {
                vector<Conteo> frecs = frecuenciasBase;
                auto start = chrono::high_resolution_clock::now();
                problema.solverConstructivo(varsInicial, frecs);
                auto end = chrono::high_resolution_clock::now();
                if (runH)
                {
                    cH.push_back(problema.calcularCosto(varsInicial));
                    tH.push_back(chrono::duration<double>(end - start).count());
                }
            }

            if (runLS)
            {
                vector<TBool> varsParaLS = varsInicial;
                auto start = chrono::high_resolution_clock::now();
                problema.busquedaLocal(varsParaLS);
                auto end = chrono::high_resolution_clock::now();
                tLS.push_back(chrono::duration<double>(end - start).count());
                cLS.push_back(problema.calcularCosto(varsParaLS));
            }

            if (runILS)
            {
                vector<TBool> varsParaILS = varsInicial;
                auto start = chrono::high_resolution_clock::now();
                problema.busquedaLocalIterada(varsParaILS, 20, gen);
                auto end = chrono::high_resolution_clock::now();
                tILS.push_back(chrono::duration<double>(end - start).count());
                cILS.push_back(problema.calcularCosto(varsParaILS));
            }

            if (runTS)
            {
                vector<TBool> varsParaTS = varsInicial;
                auto start = chrono::high_resolution_clock::now();
                int tenure = 7 + (datosFormula.first / 10);
                problema.busquedaTabu(varsParaTS, 100, tenure);
                auto end = chrono::high_resolution_clock::now();
                tTS.push_back(chrono::duration<double>(end - start).count());
                cTS.push_back(problema.calcularCosto(varsParaTS));
            }

            if (runSA)
            {
                vector<TBool> varsParaSA = varsInicial;
                auto start = chrono::high_resolution_clock::now();
                problema.recocidoSimulado(varsParaSA, gen, 10.0, 0.98, 10);
                auto end = chrono::high_resolution_clock::now();
                tSA.push_back(chrono::duration<double>(end - start).count());
                cSA.push_back(problema.calcularCosto(varsParaSA));
            }

            if (runGRASP)
            {
                vector<TBool> varsParaGRASP(datosFormula.first, TBool::Unknown);
                auto start = chrono::high_resolution_clock::now();
                problema.busquedaGRASP(varsParaGRASP, 5, 0.2, gen, frecuenciasBase);
                auto end = chrono::high_resolution_clock::now();
                tGRASP.push_back(chrono::duration<double>(end - start).count());
                cGRASP.push_back(problema.calcularCosto(varsParaGRASP));
            }

            if (runAG)
            {
                vector<TBool> varsParaAG = varsInicial;
                double probMutacion = 1.0 / datosFormula.first;
                auto start = chrono::high_resolution_clock::now();
                problema.algoritmoGenetico(varsParaAG, 50, 100, 0.85, probMutacion, gen);
                auto end = chrono::high_resolution_clock::now();
                tAG.push_back(chrono::duration<double>(end - start).count());
                cAG.push_back(problema.calcularCosto(varsParaAG));
            }

            if (runAM)
            {
                vector<TBool> varsParaAM = varsInicial;
                auto start = chrono::high_resolution_clock::now();
                problema.algoritmoMemetico(varsParaAM, 30, 50, 10.0, 15, gen);
                auto end = chrono::high_resolution_clock::now();
                tAM.push_back(chrono::duration<double>(end - start).count());
                cAM.push_back(problema.calcularCosto(varsParaAM));
            }

            if (runSS)
            {
                vector<TBool> varsParaSS = varsInicial;
                auto start = chrono::high_resolution_clock::now();
                problema.busquedaDispersa(varsParaSS, 10, 10.0, gen);
                auto end = chrono::high_resolution_clock::now();
                tSS.push_back(chrono::duration<double>(end - start).count());
                cSS.push_back(problema.calcularCosto(varsParaSS));
            }

            if (runACO)
            {
                vector<TBool> varsParaACO = varsInicial;
                auto start = chrono::high_resolution_clock::now();
                problema.optimizacionColoniaHormigas(varsParaACO, 20, 10.0, gen);
                auto end = chrono::high_resolution_clock::now();
                tACO.push_back(chrono::duration<double>(end - start).count());
                cACO.push_back(problema.calcularCosto(varsParaACO));
            }
        }

// 5. Impresión de Resultados Modular (Protegida por Critical)
#pragma omp critical
        {
            string nombreCorto = (nombreArchivo.length() > 17) ? ".." + nombreArchivo.substr(nombreArchivo.length() - 15) : nombreArchivo;

            // Se quitó la impresión de datosFormula.second (Casos)
            cout << left << setw(18) << nombreCorto << " | " << setw(8) << "--------";

            if (runH)
                cout << " | " << setw(wC) << (int)round(promedio(cH)) << " | " << setw(wT) << formatearMedida(promedio(tH), desviacionEstandar(tH, promedio(tH)));
            if (runLS)
                cout << " | " << setw(wC) << (int)round(promedio(cLS)) << " | " << setw(wT) << formatearMedida(promedio(tLS), desviacionEstandar(tLS, promedio(tLS)));
            if (runILS)
                cout << " | " << setw(wC) << (int)round(promedio(cILS)) << " | " << setw(wT) << formatearMedida(promedio(tILS), desviacionEstandar(tILS, promedio(tILS)));
            if (runTS)
                cout << " | " << setw(wC) << (int)round(promedio(cTS)) << " | " << setw(wT) << formatearMedida(promedio(tTS), desviacionEstandar(tTS, promedio(tTS)));
            if (runSA)
                cout << " | " << setw(wC) << (int)round(promedio(cSA)) << " | " << setw(wT) << formatearMedida(promedio(tSA), desviacionEstandar(tSA, promedio(tSA)));
            if (runGRASP)
                cout << " | " << setw(wC) << (int)round(promedio(cGRASP)) << " | " << setw(wT) << formatearMedida(promedio(tGRASP), desviacionEstandar(tGRASP, promedio(tGRASP)));
            if (runAG)
                cout << " | " << setw(wC) << (int)round(promedio(cAG)) << " | " << setw(wT) << formatearMedida(promedio(tAG), desviacionEstandar(tAG, promedio(tAG)));
            if (runAM)
                cout << " | " << setw(wC) << (int)round(promedio(cAM)) << " | " << setw(wT) << formatearMedida(promedio(tAM), desviacionEstandar(tAM, promedio(tAM)));
            if (runSS)
                cout << " | " << setw(wC) << (int)round(promedio(cSS)) << " | " << setw(wT) << formatearMedida(promedio(tSS), desviacionEstandar(tSS, promedio(tSS)));
            if (runACO)
                cout << " | " << setw(wC) << (int)round(promedio(cACO)) << " | " << setw(wT) << formatearMedida(promedio(tACO), desviacionEstandar(tACO, promedio(tACO)));

            cout << "\n";
            cout.flush();

            // --- ESCRITURA EN CSV ---
            ofstream archivoCSVSync(nombreCSV, ios::app);
            if (archivoCSVSync.is_open())
            {
                // Se quitó la impresión de datosFormula.second (Casos)
                archivoCSVSync << nombreArchivo << ",--------";
                if (runH)
                    archivoCSVSync << "," << (int)round(promedio(cH)) << "," << formatearMedida(promedio(tH), desviacionEstandar(tH, promedio(tH)));
                if (runLS)
                    archivoCSVSync << "," << (int)round(promedio(cLS)) << "," << formatearMedida(promedio(tLS), desviacionEstandar(tLS, promedio(tLS)));
                if (runILS)
                    archivoCSVSync << "," << (int)round(promedio(cILS)) << "," << formatearMedida(promedio(tILS), desviacionEstandar(tILS, promedio(tILS)));
                if (runTS)
                    archivoCSVSync << "," << (int)round(promedio(cTS)) << "," << formatearMedida(promedio(tTS), desviacionEstandar(tTS, promedio(tTS)));
                if (runSA)
                    archivoCSVSync << "," << (int)round(promedio(cSA)) << "," << formatearMedida(promedio(tSA), desviacionEstandar(tSA, promedio(tSA)));
                if (runGRASP)
                    archivoCSVSync << "," << (int)round(promedio(cGRASP)) << "," << formatearMedida(promedio(tGRASP), desviacionEstandar(tGRASP, promedio(tGRASP)));
                if (runAG)
                    archivoCSVSync << "," << (int)round(promedio(cAG)) << "," << formatearMedida(promedio(tAG), desviacionEstandar(tAG, promedio(tAG)));
                if (runAM)
                    archivoCSVSync << "," << (int)round(promedio(cAM)) << "," << formatearMedida(promedio(tAM), desviacionEstandar(tAM, promedio(tAM)));
                if (runSS)
                    archivoCSVSync << "," << (int)round(promedio(cSS)) << "," << formatearMedida(promedio(tSS), desviacionEstandar(tSS, promedio(tSS)));
                if (runACO)
                    archivoCSVSync << "," << (int)round(promedio(cACO)) << "," << formatearMedida(promedio(tACO), desviacionEstandar(tACO, promedio(tACO)));
                archivoCSVSync << "\n";
                archivoCSVSync.close();
            }
        }
    }

    cout << string(lineLength, '=') << "\n";
    return 0;
}