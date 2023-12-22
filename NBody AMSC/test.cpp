#include <iostream>
#include <vector>
#include <chrono>
#include <omp.h>
#include <cmath>
#include "Constants.hpp"
#include "Particle.hpp"
#include "Arrows.hpp"
#include "simulationfunctions.hpp"
#include "simulationFunctionParallel.hpp"






int main() {
    

    // Misura il tempo per la versione seriale
    
    auto start_serial = std::chrono::high_resolution_clock::now();
    simulationfunctions::doSim();
    auto end_serial = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> elapsed_serial = end_serial - start_serial;

    // Misura il tempo per la versione parallela
    auto start_parallel = std::chrono::high_resolution_clock::now();
    simulationFunctionsParallel::doParallelSim();
    //stepSimParallel(particles);
    auto end_parallel = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> elapsed_parallel = end_parallel - start_parallel;

    // Calcola e stampa lo speedup
    double speedup = elapsed_serial.count() / elapsed_parallel.count();
    std::cout << "Tempo di esecuzione seriale: " << elapsed_serial.count() << " secondi\n";
    std::cout << "Tempo di esecuzione parallela: " << elapsed_parallel.count() << " secondi\n";
    std::cout << "Speedup: " << speedup << "\n";

    return 0;
}
