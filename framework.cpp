#include <vector>
#include <iostream>
#include <fstream>
#include <cmath>

//#############################################################################
// Calculate vector dSigma/dt
//#############################################################################

void CalcFdens(std::vector<double> radii,
               std::vector<double> sdens,
               std::vector<double> fdens,
               double visc)
{
    size_t N = radii.size();        // Size of vectors
    std::vector<double> temp(N);    // Temporary vector

    // temp = sqrt(R)*d/dr(Signa*sqrt(R))
    temp[0] =
            sqrt(radii[0])*(sdens[1]*sqrt(radii[1]) - sdens[0]*sqrt(radii[0]))/
            (radii[1] - radii[0]);
    for (int i = 1; i < N - 1; i++)
        temp[i] =
                sqrt(radii[i])*(sdens[i + 1]*sqrt(radii[i + 1]) -
                                sdens[i - 1]*sqrt(radii[i - 1]))/
                (radii[i + 1] - radii[i - 1]);
    temp[N - 1] =
            sqrt(radii[N - 1])*(sdens[N - 1]*sqrt(radii[N - 1]) -
                                sdens[N - 2]*sqrt(radii[N - 2]))/
            (radii[N - 1] - radii[N - 2]);

    // fdens = (3*visc/R)*d/dr(temp)
    fdens[0] = (3.0*visc/radii[0])*(temp[1] - temp[0])/(radii[1] - radii[0]);
    for (int i = 1; i < N - 1; i++)
        fdens[i] = (3.0*visc/radii[i])*(temp[i + 1] - temp[i - 1])/
                   (radii[i + 1] - radii[i - 1]);
    fdens[N - 1] = (3.0*visc/radii[N - 1])*(temp[N - 1] - temp[N - 2])/
                   (radii[N - 1] - radii[N - 2]);
}

//#############################################################################
// Do one time step using Euler method
//#############################################################################

void DoTimeStep(std::vector<double> radii,
                std::vector<double> sdens,
                std::vector<double> fdens,
                double dt,
                double visc)
{
    // Calculate dSigma/dt
    CalcFdens(radii, sdens, fdens, visc);

    size_t N = radii.size();
    // Update using Euler's method
    for (int i = 0; i < N; i++)
        sdens[i] += dt*fdens[i];
}

//#############################################################################
// Main function
//#############################################################################

int main()
{
    size_t N = 128;       // Number of radial bins
    double rMin = 0.4;    // Inner radius
    double rMax = 2.5;    // Outer radius
    double visc = 1.0*pow(10,-5); // Viscosity
    double maxT = 100.0;  // Maximum simulation time
    double dt = 0.01;     // Time step

    std::vector<double> radii(N);   // Vector of radii
    std::vector<double> sdens(N);   // Vector of surface densities
    std::vector<double> fdens(N);   // Vector of dSigma/dt

    // Fill vector of radii
    for (int i = 0; i < N; i++)
        radii[i] = rMin + (double) i*(rMax - rMin)/(double) (N - 1);

    // Fill vector of surface densities
    for (int i = 0; i < N; i++)
        sdens[i] = 0.001/sqrt(radii[i]);

    // Take time steps until t = maxT
    double t = 0.0;
    while (t < maxT) {
        DoTimeStep(radii, sdens, fdens, dt, visc);
        t += dt;
    }

    // Output result to text file
    std::ofstream output_file("./output.txt");
    for (int i = 0; i < N; i++)
        output_file << radii[i] << " " << sdens[i] << std::endl;
    output_file.close();

    cout << dt << endl;
    cout << sdens << endl;
    cout << i << endl;
    cout << i << endl;

    return 0;
}
