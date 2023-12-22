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
        std::uniform_real_distribution<double> distribution(-100.0, 100.0); // Cambia il range se necessario
        std::uniform_real_distribution<double> distributionVel(-1.0,1.0);
        std::uniform_real_distribution<double> massDistribution(10.0, 1000000.0); // Massa compresa tra 0.1 e 10.0 (valori arbitrari)
 
        for (unsigned int i = 0; i < numberOfParticles; ++i) {
            // Genera valori casuali per posizione, velocità, massa, ecc.
            Arrows<dim> randomPosition;
            Arrows<dim> randomVelocity;
           
            
            for (unsigned int j = 0; j < dim; ++j) {
                randomPosition[j] = distribution(gen);
                randomVelocity[j] = distributionVel(gen);
                // Genera valori casuali per altre proprietà se necessario
            }

            double randomMass = massDistribution(gen);

            // Crea e aggiungi una particella con valori casuali alla collezione
            particles.push_back(Particle<dim>(i + 1, randomPosition, randomVelocity, Arrows<dim>(), Arrows<dim>(), randomMass));
     }


        


        return particles;
    }

   
   

     //INSERIMENTO BLOCCHI PARALLELI:
     static void doParallelSim()
   {
    
      
       std::vector<Particle<dim>> particles = std::vector<Particle<dim>>();
       particles.reserve(numberOfParticles);
       particles = generateParticles();
       stepParallelSim(particles);          
         
   }
   
  
   
    
    
    static void stepParallelSim(std::vector<Particle<dim>>& particles)
   { 

        double local_dt = dt;
        
        

    std::vector<Arrows<dim>> tempCoefficients(particles.size());

    std::vector<bool> collisions(particles.size(), false);

    #pragma omp parallel num_threads(4) shared(particles, local_dt, tempCoefficients, collisions)
    {              
        std::vector<double> positions = std::vector<double>();
        positions.reserve(numberOfParticles);

        for(unsigned int c=0 ; c<cycles;c++){
            #pragma omp for schedule(static,numberOfParticles/omp_get_num_threads())
            for(size_t i=0; i< particles.size(); ++i){
                particles[i].setToZero();
            }
              

                    //FOR EQUIVALENTE DI MEMPOS
                #pragma omp for schedule(static,numberOfParticles/omp_get_num_threads()) 
                for(size_t i=0; i<particles.size();++i)
                    {                                                                                           
                        for(unsigned int j=0; j<dim; ++j){
                        positions.emplace_back(particles[i].getPositionCoordinate(j));
                         }
                    }

                #pragma omp for schedule(static)
                for (unsigned int i = 0; i < particles.size(); ++i) {
                    Arrows<dim> temp1 = Arrows<dim>();

                    for (unsigned int j = 0; j < particles.size(); ++j) {
                        if (i == j) { continue; }
                            temp1 += particles[i].calcCoefficients(particles[j]);
                        if (particles[i].collision(particles[j])) { collisions[i] = true; }
                    }
                        
                    tempCoefficients[i] = temp1;
                
                }

                #pragma omp for schedule(static)
                for (unsigned int i = 0; i < particles.size(); ++i) {
                    particles[i].coefficientsSetter(tempCoefficients[i]);
                    if (collisions[i]) {
                        particles[i].calcAccellerationAfterCollision();
                    }
                    else {
                        particles[i].calcAccelleration();
                    }
            
                }

                #pragma omp for schedule(static)
                for (unsigned int i = 0; i < particles.size(); ++i) {
                    particles[i].updatePosition(local_dt);
                }

                #pragma omp for schedule(static)
                 for (unsigned int i = 0; i < particles.size(); ++i) {
                    particles[i].updateVelocity(local_dt);
                }

            }//alla fine dei cicli calcolati in base al total time e dt viene scritto su file positions
        writePositionToTXT(positions); 
    }
        
 } 
           
         
    //Scrivi le coordinate delle particelle su file TXT. Necessario per implementare la graficazione della simulazione
    static void writePositionToTXT(const std::vector<double> positions) {
        std::ofstream outputFile("positions.txt");
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
