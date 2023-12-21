


#ifndef CONSTANTS_H
#define CONSTANTS_H

//Dimensione del problema 1D 2D 3D
const unsigned int dim=3;


//Costante Gravitazionale Universale G

const double G = -6.67408e-11;



//Delta temporale espresso in secondi
double dt=0.1;
      
//Tempo totale della simulazione espresso in secondi
 double totalTime =1000;


//Numeri di cicli
const int cycles = static_cast<int>(std::floor(totalTime/dt));

//Numero di particelle generate                 //modifica del 17 dicembre 2023, è stata creata una costante per definire 
const int numberOfParticles = 50;              //il numero di particelle da inizializzare per la simulazione.


#endif  // CONSTANTS_H
