#ifndef PARTICLE_UTIL_HPP
#define PARTICLE_UTIL_HPP

#include "Particle.hpp"
#include "Constants.hpp"
#include <vector>
#include <random>
#include <iostream>
#include <fstream>
#include <string>



class simulationfunctions{
public:
    
    
    // Metodo statico per generare una collezione di oggetti Particle
    static std::vector<Particle<dim>> generateParticles() {
        std::vector<Particle<dim>> particles;
        //PER DEBUGGING
        /*particles.push_back(Particle<dim>(1, Arrows<dim>({1.0, 2.0, 3.0}), Arrows<dim>({0.0, 0.0, 0.0}), Arrows<dim>({0.0,0.0,0.0}), Arrows<dim>({0.0,0.0,0.0}), 5.0));
        particles.push_back(Particle<dim>(2, Arrows<dim>({3.0, 2.0, 1.0}), Arrows<dim>({0.0, 0.0, 0.0}), Arrows<dim>({0.0,0.0,0.0}), Arrows<dim>({0.0,0.0,0.0}), 7.0));
        particles.push_back(Particle<dim>(3, Arrows<dim>({4.0, 5.0, 6.0}), Arrows<dim>({0.0, 0.0, 0.0}), Arrows<dim>({0.0,0.0,0.0}), Arrows<dim>({0.0,0.0,0.0}), 10.0));
        particles.push_back(Particle<dim>(4, Arrows<dim>({6.0, 5.0, 4.0}), Arrows<dim>({0.0, 0.0, 0.0}), Arrows<dim>({0.0,0.0,0.0}), Arrows<dim>({0.0,0.0,0.0}), 1.0));*/
        // Metodo che genera un numero prefissato di particelle di dimensione DIM con valori casuali
        std::random_device rd;
        std::mt19937 gen(rd());
        std::uniform_real_distribution<double> distribution(-10.0, 10.0); // Cambia il range se necessario
        std::uniform_real_distribution<double> massDistribution(0.1, 10.0); // Massa compresa tra 0.1 e 10.0 (valori arbitrari)

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


        return particles;
    }

    //Metodo che esegue la simulazione. Tale metodo provvederà a chiamare il generatore di Particles. PER ORA SOLVER BASE!!!
    static void doSim(){
        
        std::vector<Particle<dim>> particles = std::vector<Particle<dim>>();

        std::vector<double> positions = std::vector<double>();

        particles = generateParticles();
        for(Particle particle : particles){
            particle.setToZero();
            particle.print();
            std::cout<<" "<<std::endl;
        }

        
        
        for(unsigned int i=0 ; i<cycles; i++){
            std::cout<<"!!!! INIZIO DELLA SIMULAZIONE !!!!"<<std::endl;
            std::cout<<" "<<std::endl;
            std::cout<<"CICLO "<< i << ":" << std::endl;
            std::cout<<" "<<std::endl;
            
            stepSim(particles,positions);
            updateSim(particles);
        }
        writePositionToTXT(positions);
    
    }
 
    
    //Salvataggio delle posizioni su una vector. utile per generare il file txt che dovrà essere letto per l'animazione grafica
    static void memPos(const std::vector<Particle<dim>> particles, std::vector<double>& positions){
        for(Particle particle : particles){
            for(unsigned int i=0; i<dim; ++i){
                positions.emplace_back(particle.getPositionCoordinate(i));
            }
        }
    }

    //Algoritmo Base
    static void stepSim(std::vector<Particle<dim>>& particles,std::vector<double>& positions){  
        memPos(particles,positions);//Riempi le posizioni a ciclo
        for(unsigned int i = 0; i<particles.size(); ++i){
            Arrows<dim> temp = Arrows<dim>();
            for(unsigned int j=0; j<particles.size(); ++j){
                if(i==j){continue;}
                temp += particles[i].calcCoefficients(particles[j]);
            }
            particles[i].coefficientsSetter(temp);
            particles[i].calcAccelleration();
        }
    } 

    //Algoritmo ridotto
    static void reducedSim(std::vector<Particle<dim>>& particles){
        unsigned int middleIndex = particles.size() / 2;
        for(unsigned int i = 0; i < middleIndex; ++i){
            Arrows<dim> temp = Arrows<dim>();
            Arrows<dim> temp2 = Arrows<dim>();
            for(unsigned int j = i + 1; j < particles.size(); ++j){
                temp += particles[i].calcCoefficients(particles[j]);
                particles[i].coefficientsSetter(temp);
                temp2 = particles[j].getCoefficients() - particles[i].calcCoefficients(particles[j]);
                particles[j].coefficientsSetter(temp2);
            }
            particles[i].calcAccelleration();
            particles[particles.size()-i-1].calcAccelleration();
        }
    }
 
        
    static void updateSim(std::vector<Particle<dim>>& particles){ 
        double local_dt = dt;
        for(Particle<dim>& particle : particles){
            particle.updateVelocity(local_dt);
            particle.updatePosition(local_dt);
            particle.print();
            std::cout<<" "<<std::endl;
            particle.setToZero();
        }
    } 
    
   
    void solveCollisions(std::vector<Particle<dim>> particles){
        
    }
    
    //Determinazione del numero di cicli. Deve essere fornito un lasso di tempo totale e un dt
    unsigned int numberOfCycles(double totalTime, double dt){
        double temp = totalTime/dt;
        return static_cast<unsigned int>(std::floor(temp));
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
