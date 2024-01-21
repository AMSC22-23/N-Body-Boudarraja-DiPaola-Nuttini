#include "Arrows.hpp"
#include "Constants.hpp"
#include <iostream>
#include <omp.h>
#ifndef PARTICLE_H
#define PARTICLE_H
template <unsigned int dim>
class Particle {
   private:
       unsigned int ID;
       Arrows<dim> position;
       Arrows<dim> velocity;
       Arrows<dim> accelleration;
       Arrows<dim> coefficients;
       double mass;
  
   public:
       //Costruttore di Default
       Particle() : ID(0), position(), velocity(), accelleration(), coefficients(), mass(0.0) {}


       //Costruttore con parametri
       Particle(unsigned int id, const Arrows<dim>& pos, const Arrows<dim>& vel, const Arrows<dim>& acc, const Arrows<dim>& coeff, double m): ID(id), position(pos), velocity(vel), accelleration(acc),coefficients(coeff), mass(m) {}


       //Metodo per il calcolo dei coefficienti che richiede come argomento un'altra particella e torna una Arrow di coefficienti
       Arrows<dim> calcCoefficients(const Particle& otherParticle){
           //Calcolo in primis una arrow che abbia come componenti le differenze tra le coordinate
           Arrows<dim> differences = Arrows<dim>();
             
           
           for(unsigned int i=0; i<dim; ++i){
               differences[i]=position[i]-otherParticle.position[i];
           }
           double distance = differences.normEclidea();
           return differences * otherParticle.mass /(distance*distance*distance);
       }
      
      
       //Meteodo per il calcolo delle accellerazioni
       //di fatto si occupa di andare a moltiplicare i coefficienti della particella per la costante G, definita in Constants.hpp
       void calcAccelleration(){
           accelleration = coefficients*G;
       }

        //Metodo per calcolo Accellerazione dopo collisione. AGGIUNTO DA KARIM
        void calcAccellerationAfterCollision(){
            accelleration = coefficients*G;
            accelleration = accelleration * (-1);
        }



       //Metodo per check collisioni FATTO DA KARIM
       bool collision(const Particle<dim>& otherParticle){
          //Compara posizioni tra particelle
          for(unsigned int i=0 ; i<dim ; ++i){
            if(position[i] != otherParticle.position[i]){
                return false;
            }
          }
          return true;
       }


       //Usato principalmente per motivi di Debugging
       void print(){
           std::cout<<"POSIZIONI"<<std::endl;
           std::cout<<position<<std::endl;
           std::cout<<"ACCELLERAZIONI"<<std::endl;
           std::cout<<accelleration<<std::endl;
           std::cout<<"COEFFICIENTI"<<std::endl;
           std::cout<<coefficients<<std::endl;
           std::cout<<"FORZA"<<std::endl;
           Arrows<dim> forces = Arrows<dim>();
           forces = accelleration*mass;
           std::cout<<forces<<std::endl;


       }




       //Metodo per calcolare future Velocità. Si basa sulle leggi di moto del m.u.a. La velocità V_t sarà data da V_t_0 + Accellerazione * dt
       Arrows<dim> calcNextVelocity(double& dt){
           
           return velocity + accelleration * dt;
       }


       //Metodo per calcolare future Posizioni. Si basa sulle leggi di moto del m.u.a. La posizione S_t sarà data da S_t_0 + V_t_0 * Accellerazione + 0.5 * Accellerazione * dt * dt
       Arrows<dim> calcNextPosition(double& dt){
           Arrows<dim> futurePos = Arrows<dim>();
           return futurePos = position + velocity * dt + accelleration * dt * dt * 0.5;
       }


       //Metodo per aggiornare Velocità.
       void updateVelocity(double& dt){
           velocity = calcNextVelocity(dt);
       }


       //Metodo per aggiornare Posizioni
       void updatePosition(double& dt){
           position = calcNextPosition(dt);
       }


       //Metodo per set valori coefficienti
       void setToZero(){
           if constexpr (dim == 1){
               /*coefficients = Arrows<dim>({0.0});
               accelleration = Arrows<dim>({0.0});*/                    //MODIFICHE APPORTATE AI FINI DI ALLEGGERIMENTO
               coefficients*= 0;
               accelleration*=0;
           }
           else if constexpr (dim==2)
           {
               //coefficients = Arrows<dim>({0.0, 0.0});
               //accelleration = Arrows<dim>({0.0, 0.0});
               coefficients*=0;
               accelleration*=0;
           }
           else if constexpr(dim==3){
              // coefficients = Arrows<dim>({0.0, 0.0, 0.0});
              // accelleration = Arrows<dim>({0.0, 0.0, 0.0});
              coefficients  *= 0;
              accelleration *=0;
           }


       }
      
       void coefficientsSetter (Arrows<dim> otherCoefficients){
            
           coefficients=otherCoefficients; 
       }
       

       Arrows<dim> getCoefficients() const {
           return coefficients;
        }
        
        //MODIFICHE 17/12
        double getPositionCoordinate(const unsigned int i) const{
            return position[i];
        }

        auto getID() const{
            return ID;
        }

         auto getMass() const{
            return mass;
        }  
      
};
#endif


