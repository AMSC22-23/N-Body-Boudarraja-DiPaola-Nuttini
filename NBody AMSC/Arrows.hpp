
//Committed by Boudarraja Karim

#ifndef ARROWS_HPP
#define ARROWS_HPP
#include <array>
#include <cmath>
#include <ostream>
#include <istream>
#include <iostream>

#include "Constants.hpp"


template <unsigned int dim>
class Arrows
{
private:
   //La classe è formata da un array di double di dimensione pari a DIM. il ragionamento si basa sul fatto che i vettori siano effettivamente composti
   //da "componenti" vettori chew sono le proiezioni del vettore lungo gli assi cartesiani x,y,z (a seconda di DIM)
   std::array<double,dim> components;
public:
   //Costruttore di default della classe Arrow che inizializza a 0 gli elementi.
   Arrows<dim>() : components({0}) {}


   //Costruttore che inizializza arrow sulla base dell'array passato come argomento
   Arrows<dim>(const std::array<double,dim>& arr) : components(arr) {}


   //Metodo per il ritorno dello specifico componente della Arrow. Qui torna un riferimento double& dato il riferimento int& index.
   double& operator[](const unsigned int& index){return components[index];}
   //Metodo per il ritorno dello specifico componente della Arrow. Qui torna un valore double dato il riferimento int& index.
   double operator[](const unsigned int& index) const {return components[index];}


   //Metodo Arrow per scalare. Richiede come argomento un riferimento a uno scalare double e successivamente moltiplica ogni singolo componente per scal.
   //Il metodo richiede implementazione friend al fine di ottenere come argomenti un Arrow e uno scal
   Arrows<dim>& operator*=(const double& scal){
       /*for(unsigned int i=0; i<dim; ++i)
           components[i] *= scal;
       return *this;*/
       components[0]*= scal;
       if(dim==2)
       {
        components[1]*= scal;
       }
       if (dim==3)
       {
        components[1]*= scal;

        components[2]*= scal;
       }
       return *this;
   }
   //Formato friend
   friend Arrows<dim> operator*(Arrows<dim> arrow, const double& scal){
       arrow *= scal;
       return arrow;
   }
   //Metodo Arrow divisione per scalare. Richiede come argomento un riferimento a uno scalare double e successivamente divide ogni singolo componente per scal.
   //Il metodo richiede implementazione friend al fine di ottenere come argomenti un Arrow e uno scal
   Arrows<dim>& operator/=(const double& scal){
       for(unsigned int i=0; i<dim; ++i)
           components[i] /= scal;
       return *this;
   }
   //Formato friend
   friend Arrows<dim> operator/(Arrows<dim> arrow, const double& scal){
       arrow /= scal;
       return arrow;
   }


   //Metodo per somma tra Arrow. Il metodo logicamente richiederà anche un metodo friend che richiede due arrows.
   //Il metodo operator += torna Arrow<DIM>& e richiede come argomento un riferimento costante ad un altra Arrow.
   Arrows<dim>& operator+=(const Arrows<dim>& secondArrow){
       for(unsigned int i=0; i<dim; ++i)
           components[i]+=secondArrow[i];
       return *this;
   }
   


   //Formato friend. Usato per sommare due Arrows quindi richiede due Arrows come argomenti e usa l'operator += definito in precedenza.
   //Gli argomenti richiesti sono il primo Arrow da sommare e il riferimento const al secondo Arrow da sommare (logicamente const & in quanto il secondo non varia).
   //Ritorno il primo Arrow quindi così sommato.
   friend Arrows<dim> operator+(Arrows<dim> firstArrow, const Arrows<dim>& secondArrow){
       firstArrow += secondArrow;
       return firstArrow;
   }


   //Metodo per differenza tra Arrow. Il metodo logicamente richiederà anche un metodo friend che richiede due arrows.
   //Il metodo operator += torna Arrow<DIM>& e richiede come argomento un riferimento costante ad un altra Arrow.
   Arrows<dim>& operator-=(const Arrows<dim>& secondArrow){
       for(unsigned int i=0; i<dim; ++i)
           components[i]-=secondArrow[i];
       return *this;
   }


   //Formato friend. Usato per eseguire differenza due Arrows quindi richiede due Arrows come argomenti e usa l'operator -= definito in precedenza.
   //Gli argomenti richiesti sono il primo Arrow dal quale sottraggo il riferimento const al secondo Arrow (logicamente const & in quanto il secondo non varia).
   //Ritorno il primo Arrow quindi così differernziato.
   friend Arrows<dim> operator-(Arrows<dim> firstArrow, const Arrows<dim>& secondArrow){
       firstArrow -= secondArrow;
       return firstArrow;
   }


   //Prodotto scalare tra due vettori che torna uno scalare. //controlla che sia stato modificato correttamente
   double operator*(const Arrows<dim>& secondArrow) const{
       double product = 0.;
       for (unsigned int i=0; i<dim; ++i)
           product += components[i]*secondArrow[i];
       return product;
   }


   //Prodotto vettoriale tra due vettori. Anche questo metodo richiederà l'uso di un metodo friend che segua la logica vista in precedenza
   Arrows<dim>& operator^=(const Arrows<dim>& secondArrow){
       double x=0,y=0,z=0;
       if constexpr(dim == 3)
       {
           x = components[1]*secondArrow[2] - components[2]*secondArrow[1];
           y = components[2]*secondArrow[0] - components[0]*secondArrow[2];
           z = components[0]*secondArrow[1] - components[1]*secondArrow[0];
          
       }
       else if constexpr (dim == 2)
       {
           x = components[0]*secondArrow[1] - components[1]*secondArrow[0];
       }
       else if constexpr(dim == 1)
       {
           x = 0;
       }
       components = {x, y, z};
       return *this;
   }


   //Forma Friend che richiede i due vettori
   friend Arrows<dim> operator^(Arrows<dim> firstArrow, Arrows<dim>& secondArrow){
       firstArrow ^= secondArrow;
       return firstArrow;
   }


   //Metodo di stampa che richiede operator<<
   friend std::ostream& operator<<(std::ostream& out, Arrows<dim>& arrow){
       out << "[";
       for(unsigned int i=0; i<dim; ++i)
           out << arrow[i] << " ";
       out << "]";
       return out;
   }
  
  //bisogna rivedere questo metodo 
   double normEclidea() const{
       double tot=0.;
       for (unsigned int i = 0; i < dim; ++i)
           tot += components[i] * components[i];
       return std::sqrt(tot);
   }
};
#endif // ARROWS_HPP
