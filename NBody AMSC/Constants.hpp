/*
ATTENZIONE!!!
Il seguente header serve a definire le costanti del problema e deve essere ancora implementato in una maniera coerente 
per poter essere utilizzato intelligentemente nel programma

(#ifndef CONSTANTS_H
#define CONSTANTS_H


//Dimensione del problema 1D 2D 3D
const unsigned int dim=3;


//Costante Gravitazionale Universale G




//Delta temporale espresso in secondi
const double dt=1;
      
//Tempo totale della simulazione espresso in secondi
const double totalTime =5;


//Numeri di cicli
const int cycles = static_cast<int>(std::floor(totalTime/dt));

//Numero di particelle generate                 //modifica del 17 dicembre 2023, è stata creata una costante per definire 
const int numberOfParticles = 100;              //il numero di particelle da inizializzare per la simulazione.


#endif  // CONSTANTS_H
