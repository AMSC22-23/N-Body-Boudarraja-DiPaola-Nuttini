#include "Particle.hpp"
#include "Arrows.hpp"
#include "Constants.hpp"
#include <iostream>


using namespace std;


// Function that calculates the multipole expansion of a particle
Arrows<dim> calcMultipoleExpansion(const Particle<dim>& particle) {
   Arrows<dim> coefficients=Arrows<dim>();


   for (int i = 0; i < 2; ++i) {
       coefficients[i] = particle.getMass() / (4.0 * M_PI);
   }


   return coefficients;
}


// Function that calculates the local expansion of a particle
Arrows<dim> calcLocalExpansion(const Particle<dim>& particle) {
   Arrows<dim> forces=Arrows<dim>();


   double distance = 1.0;


   for (int i = 0; i < 2; ++i) {
       forces[i] = particle.getMass() * particle.position[i] / distance;
   }


   return forces;
}


// Function that calculates the multipole-to-local translation
Arrows<dim> M2L(const Arrows<dim>& a) {
   Arrows<dim> b= Arrows<dim>


   b = a / (4.0 * M_PI);


   return b;
}


// Function that calculates the local-to-local translation
Arrows<dim> L2L(const Arrows<dim>& b) {
   Arrows<dim> a= Arrows<dim>();


   a = b * 4.0 * M_PI;


   return a;
}


// Function that calculates the multipole-to-multipole translation
Arrows<dim> M2M(const Arrows<dim>& a1, const Arrows<dim>& a2) {
   Arrows<dim> a=Arrows<dim>();


   a = a1 + a2;


   return a;
}


// Function that evaluates the potential at a point
double evaluate(Arrows<dim>& position, Arrows<dim>& multipoleCoefficients) {
   double potential = 0.0;


   for (int i = 0; i < 2; ++i) {
       potential += multipoleCoefficients[i] / (4.0 * M_PI * (position[i] - multipoleCoefficients[i]));
   }


   return potential;
}


void fmm(const std::vector<Particle<dim>>& particles, std::vector<double>& u) {
   unsigned int J = floor(log2(particles.size()));


   // Compute the weight
   for (int k = 0; k < 1 << J; ++k) {
        for(size_t i=0; i<particles.size();++i)
            {
                particles[i].calcCoefficients(particles[k]);
            }
       
   }


   for (int L = J - 1; L >= 0; --L) {
       for (int k = 0; k < 1 << L; ++k) {
           for (int s = 0; s < 2; ++s) {
                particles[1 << L * 2 + k].coefficientsSetter( M2L(particles[1 << L * 2 + s].getCoefficients()));

           }
       }
   }


   // Evaluation
   for (int L = 2; L <= J - 1; ++L) {
       for (int k = 0; k < 1 << L; ++k) {
           for (const Particle<dim>& particle : particles) {
               particles[1 << (L + 1) * 2 + k].coefficients = M2L(particle.coefficients);
           }
       }
   }


   for (int k = 0; k < 1 << J; ++k) {
       for (const Particle<dim>& particle : particles) {
           u[k] = evaluate(particle.position, particles[k].coefficients) + particle.mass;
       }
   }
}
