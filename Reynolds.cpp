#include <iostream>
#include<cmath>

/*
Reynolds is a function that calculates the Reynolds number, a 
nondimensional number that describes the flow pattern along an 
object. Reynolds number determines if the boundary layer is 
laminar or turbulent.
*/

double Reynolds(double rho, double u, double l, double mu) {

    return (rho * u * l) / mu; // Return reynolds number
}