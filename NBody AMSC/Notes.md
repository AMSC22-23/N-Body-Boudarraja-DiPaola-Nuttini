# A few notes #
The code is still work ongoing sine there are errors that need to be fixed.
For instance 
```c++
#include "simulationFunctionParallel.hpp"
```
should be
```c++
#include "simulationsFunctionParallel.hpp"
```
(or the name of the file should be changed).

As I told at lectuer `emplace_back()` is more powerful than `push_back()` since it avoids the copy of the object. So,
```c++
            particles.push_back(Particle<dim>(i + 1, randomPosition, randomVelocity, Arrows<dim>(), Arrows<dim>(), randomMass));
```
can be more efficiently stated as
```c++
            particles.emplace_back(i + 1, randomPosition, randomVelocity, Arrows<dim>(), Arrows<dim>(), randomMass);
```
avoing the creation and copy of a temporary object.

Nice having commented the code. I would suggest to use the Doxygen style for comments. It allows the automatic generation of documentation. See for instance [here](https://www.doxygen.nl/manual/docblocks.html).

You should use `constexpr` for immutable constants. For instance
```c++
    constexpr double G = 9.81;
```
I do not understand
```c++
    //Delta temporale espresso in secondi
const double dt=1;
```
since `dt` should be a variable that you may want to change.  

I am looking forward to your presentation of the final version of the code.

