#include <iostream>
#include <cmath>
#include <chrono>

double Mach(double u, double gamma, double R, double T);
double Reynolds(double rho, double u, double l, double mu);

double AirTemp(double h);
double AirPressure(double h);
double AirDensity(double h);
double AirViscosity(double h);
double Gravity(double h);

int main() {

    // Start the timer
    auto start = std::chrono::high_resolution_clock::now();

    for(int h = 0; h <= 30000; h+=1) {
        double T = AirTemp(h);
        double rho = AirDensity(h);
        double mu = AirViscosity(h);

        double Re = Reynolds(rho, 100, 1, mu);
        double M = Mach(100, 1.4, 287.052874247, T);
    }

    // Stop the timer
    auto stop = std::chrono::high_resolution_clock::now();

    // Calculate the duration
    auto duration = std::chrono::duration_cast<std::chrono::milliseconds>(stop - start);

    // Print the duration
    std::cout << "Time taken: " << duration.count() << " milliseconds" << std::endl;

    return 0;
}
