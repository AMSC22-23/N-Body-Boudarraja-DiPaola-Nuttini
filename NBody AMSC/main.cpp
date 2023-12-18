#include <cmath>
#include "Constants.hpp"
#include "Particle.hpp"
#include "simulationfunctions.hpp"
#include "simulationFunctionParallel.hpp"

/*
il main viene principalmente utilizzato per avviare il programma, al suo interno sono contenuti gli header files delle classi utilizzate per la simulazione
tra cui 
1)Costants.hpp -> header file utilizzato per definire delle costanti che verranno utilizzate nel programma
2)Particle.hpp ->header file in cui è contenuta la classe che rappresenta il corpo di una particella
3)simulationfunctions.hpp -> La classe in cui viene eseguita effettivamente la simulazione
il seguente blocco di programma di occupa di chiamare il metodo statico doSim() che  da il via alla simulazione
*/

int main(){
    
    bool parallel=true;
    
    if (parallel==true){
        simulationFunctionsParallel::doParallelSim();
    } else {
        simulationfunctions::doSim();
    }
    
    return 0;
}
