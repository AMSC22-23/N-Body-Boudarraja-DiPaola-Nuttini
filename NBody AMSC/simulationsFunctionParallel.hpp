#ifndef PARTICLE_UTIL_HPP_PARALLEL
#define PARTICLE_UTIL_HPP_PARALLEL

#include "Particle.hpp"
#include "Constants.hpp"
#include <vector>
#include <random>
#include <iostream>
#include <fstream>
#include <string>
#include <omp.h>
#include<chrono>


class simulationFunctionsParallel{
public:
    
    
    // Metodo statico per generare una collezione di oggetti Particle
    static std::vector<Particle<dim>> generateParticles() {
        std::vector<Particle<dim>> particles;
        particles.reserve(numberOfParticles);
        //PER DEBUGGING
        
        // Metodo che genera un numero prefissato di particelle di dimensione DIM con valori casuali
        std::random_device rd;
        std::mt19937 gen(rd());
        std::uniform_real_distribution<double> distribution(-10.0, 10.0); // Cambia il range se necessario
        std::uniform_real_distribution<double> massDistribution(0.1, 10.0); // Massa compresa tra 0.1 e 10.0 (valori arbitrari)
 // #pragma omp parallel num_threads(8) 
 //   {  
  //  #pragma omp for
        for (unsigned int i = 0; i < numberOfParticles; ++i) {
            // Genera valori casuali per posizione, velocità, massa, ecc.
            Arrows<dim> randomPosition;
            Arrows<dim> randomVelocity;
            //Arrows<dim> randomAcceleration;
            //Arrows<dim> randomCoefficients;
            
            for (unsigned int j = 0; j < dim; ++j) {
                randomPosition[j] = distribution(gen);
                randomVelocity[j] = distribution(gen);
                // Genera valori casuali per altre proprietà se necessario
            }

            double randomMass = massDistribution(gen);

            // Crea e aggiungi una particella con valori casuali alla collezione
            particles.push_back(Particle<dim>(i + 1, randomPosition, randomVelocity, Arrows<dim>(), Arrows<dim>(), randomMass));
     }

  // }
        


        return particles;
    }

   
     //Salvataggio delle posizioni su una vector. utile per generare il file txt che dovrà essere letto per l'animazione grafica
    static void memPos(const std::vector<Particle<dim>> particles, std::vector<double>& positions){
        for(Particle particle : particles){
            for(unsigned int i=0; i<dim; ++i){
                positions.emplace_back(particle.getPositionCoordinate(i));
            }
        }
    }

     //INSERIMENTO BLOCCHI PARALLELI:
     static void doParallelSim()
   {
    std::vector<double> positions = std::vector<double>();
    positions.reserve(numberOfParticles);

      
       std::vector<Particle<dim>> particles = std::vector<Particle<dim>>();
       particles.reserve(numberOfParticles);
       particles = generateParticles();


      
      
      /*
      il seguente ciclo for viene utilizzato per eseguire la simulazione dell'Nbody problem sulle N particelle appena generate dal metodo 
      generateParticles(), . Il numero ci cili è stabilito nel file Constants.hpp dove il tempo totale della simulazione viene suddiviso per il numero di intervalli 
      desiderati*/

       for(unsigned int i=0 ; i<cycles;i++){
         
           stepParallelSim(particles,positions);
                                                  
    
       }
       
         writePositionToTXT(positions);                            

       
   }
   /*
   il seguente metodo viene utilizzato per calcolare per ogni step dt l'iterazione tra le particles della nostra simulazione, 
   utilizzando la openMp vengono parallelizzati tre processi:
   1) il ciclo di inizializzazione dei coefficienti delle particelle ad ogni ciclo(viene impostato a 0)
   2) vengono ricavati i coefficienti per ciascuna particella attraverso i due cicli for annidati, il thread safe è garantito dal fatto che le particelle si limitano
   ad eseguire una lettura esclusivamente delle posizioni e delle masse ai fini dei calcolo dei coefficenti
   3)Aggiornamento di velocità e posizione */
   static void stepParallelSim(std::vector<Particle<dim>>& particles,std::vector<double>& positions)
   { 
       memPos(particles,positions);//Riempi le posizioni a ciclo                       

        
         #pragma omp parallel shared( particles)
         {
            omp_set_num_threads(omp_get_max_threads());

             #pragma omp for schedule(static,omp_get_num_threads())
             for(size_t i=0; i< particles.size(); ++i){
                particles[i].setToZero();

             }
                         
            
        #pragma omp barrier  


        Arrows<dim> temp1=Arrows <dim>();
        bool collision=false;
           #pragma omp for  schedule(static,numberOfParticles/omp_get_num_threads()) private (temp1)
            for(unsigned int i = 0; i<particles.size(); ++i){
                    temp1*=0;
                
                     for(unsigned int j=0; j<particles.size(); ++j){

                        if(i==j){continue;}
                        temp1 += particles[i].calcCoefficients(particles[j]);
                        if(particles[i].collision(particles[j])){collision=true;}
                    }
                            
                    particles[i].coefficientsSetter(temp1);
                        if(collision){
                            particles[i].calcAccelleration();
                        }
                        else
                        {
                            particles[i].calcAccellerationAfterCollision();
                        }
                        
                            
            } //NB: the pragma omp barrier could be avoided because it is implicit at the end of the parallel for
            //#pragma omp barrier
            #pragma  omp for schedule (static,numberOfParticles/omp_get_num_threads()) 
            for(size_t i=0; i<particles.size();++i)
            {
                 particles[i].updateVelocity(dt);
                 particles[i].updatePosition(dt);
            }
            

     }   
           
     } //steoParallelSim e doParallelSim sono del blocco parallelo
    
   
    void solveCollisions(std::vector<Particle<dim>> particles){
        
    }
    
    //Determinazione del numero di cicli. Deve essere fornito un lasso di tempo totale e un dt
    unsigned int numberOfCycles(double totalTime, double dt){
        double temp = totalTime/dt;
        return static_cast<unsigned int>(std::floor(temp));
    }
    
    //Scrivi le coordinate delle particelle su file TXT. Necessario per implementare la graficazione della simulazione
    static void writePositionToTXT(const std::vector<double> positions) {
        std::ofstream outputFile("positionsParallel.txt");
        if (!outputFile.is_open()) {
            std::cerr << "Impossibile aprire il file!" << std::endl;
            return;
        }

        int coordinates = dim;
        int rows = positions.size() / coordinates;
        int dif = 0;
        int cycle = 0;
        int counter = cycle + 1;

        for (int i = 0; i < rows; ++i) {
            outputFile << cycle << ","; // Scrivi il numero del ciclo
            if(i >= numberOfParticles)
            {
                outputFile << i - dif << ","; //Scrivi l'ID della particella dopo il primo ciclo
            } 
            else 
            {
                outputFile << i << ","; // Scrivi l'ID della particella per il primo ciclo
            }

            for (int j = 0; j < coordinates; ++j) {
                outputFile << positions[i * coordinates + j]; // Scrivi le coordinate
                if (j != coordinates - 1) {
                    outputFile << ",";
                }
            }
            outputFile << std::endl; // Vai a capo dopo ogni riga
            

            if(i == numberOfParticles*counter - 1)
            {
                dif = counter*numberOfParticles;
                ++counter;
                ++cycle;
            }
            
        }

        outputFile.close();
    }

};

#endif // PARTICLE_UTIL_HPP
  
