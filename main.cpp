#include <iostream>
#include "constants.h"
#include <vector>
#include <fstream>
#include <cmath>
#include <cstdlib>
#include <iomanip>

using namespace std;



    void CalcFdens(std::vector<double> radii,//R
                   std::vector<double> sdens,//sigma
                   std::vector<double> fdens,//dSigma (???)
                   double visc,
                   std::vector<double> lambda,
                   std::vector<double> lambda2) //viscosity

    {
        size_t N = radii.size();        // Size of vectors
        std::vector<double> temp(N);    // Temporary vector

        //temp = sqrt(R)*d/dr(Sigma*sqrt(R))
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

        double Mj =  9.5376e-4; //mass of jupiter like planet in Msun
        int M = 1;//mass of sun in Msun
        double Ms = 2.858860e-4;
        double Rj = 1;//rad in AU
        double Rs = 1.3;
        double k = radii[0] - Rj;
        double m = radii[0] - Rs;
        double g = 1;
        double q = Mj/M;
        double r = Ms/M;

        /*lambda = rate of angular momentum per unit mass by tidal interaction
        * = ((R - Rj)/0.05R)*((0.23((Mj/Ms)^2)*G*Ms)/2R)*(R/|R-Rj|)
        * define lambda function here*/

        lambda [0] = ((k)/0.05)*((0.23*(pow(q,2)*g*M))/2*radii[0])*(pow((radii[0]/abs(k)),4));
        //for (int i = 1; i < N; i++)

        lambda2 [0] = ((m)/(0.5*Rs))*((0.23*pow(r,2)*g*M))/2*radii[0]*(pow((radii[0]/abs(m)),4));

        {



        // fdens = (3*visc/R)*d/dr(temp) + lambda + lambda2
        fdens[0] = (3.0*visc/radii[0])*(temp[1] - temp[0])/(radii[1] - radii[0]) + lambda[0] + lambda2[0];
        for (int i = 1; i < N - 1; i++)
            fdens[i] = (3.0*visc/radii[i])*(temp[i + 1] - temp[i - 1]) + lambda[i] + lambda2[i]/
                       (radii[i + 1] - radii[i - 1]);
        fdens[N - 1] = (3.0*visc/radii[N - 1])*(temp[N - 1] - temp[N - 2]) +lambda [N - 1] + lambda2[N-1]/
                       (radii[N - 1] - radii[N - 2]);
    }



    }
//#############################################################################
// Do one time step using Euler method
//#############################################################################

    void DoTimeStep(std::vector<double> radii,
                    std::vector<double> sdens,
                    std::vector<double> fdens,
                    std::vector<double> lambda,
                    std::vector<double> lambda2,
                    double dt,
                    double visc)
    {
        // Calculate dSigma/dt
       CalcFdens(radii, sdens, fdens, visc, lambda, lambda2);

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
        double visc = 1.0e-5; // Viscosity
        double maxT = 100.0;  // Maximum simulation time
        double dt = 0.01;     // Time step

        std::vector<double> radii(N);   // Vector of radii
        std::vector<double> sdens(N);   // Vector of surface densities
        std::vector<double> fdens(N);   // Vector of dSigma/dt
        std::vector<double> lambda(N); //vector of lambda
        std::vector<double> lambda2(N);

        // Fill vector of radii
        for (int i = 0; i < N; i++)
            radii[i] = rMin + (double) i*(rMax - rMin)/(double) (N - 1);

        // Fill vector of surface densities
        for (int i = 0; i < N; i++)
            sdens[i] = 0.001/sqrt(radii[i]);

        //fill lambda vector
        for (int i = 0; i <N; i++)
            lambda [i] = (radii[i] - 1)/0.05*((0.23*(pow((9.5376e-4),2)))/2*radii[i])*(radii[i]/abs(radii[i] -1));

        //fill lambda2
        for(int i = 0; i <N; i++)
            lambda2 [i] = (radii[i] - 1.3)/0.05*1.3*((0.23*(pow((9.5376e-4/3),2)))/2*radii[i])*(radii[i]/abs(radii[i] - 1.3));



        //fill vector of fdens
        for (int i=0; i< N; i++)
            fdens[i] = (1/radii[i])*radii[1]*(3*sqrt(radii[i])*radii[1]*(visc*sdens[i]*sqrt(radii[i])) - (2*lambda[i]*sdens[i]*radii[i]) - (2*lambda2[i]*sdens[i]*radii[i]));


        // Take time steps until t = maxT
        double t = 0.0;
        while (t < maxT) {
            DoTimeStep(radii, sdens, fdens, lambda, lambda2, dt, visc);
            t += dt;
        }


        // Output result to text file
        std::ofstream output_file("./output.txt");
        for (int i = 0; i < N; i++)
            output_file << std::fixed << std::setprecision(20) << fdens[i] << " " << radii[i] << std::endl; //<< sdens[i] << " " << lambda[i] << " " << lambda2[i] << std::endl;
        output_file.close(); //outputs radius of disc and surface area to text file


        return 0;
    } //need to make a plot of surface density vs radius, output to file
//python to read file into graph?
//mathplot lib
//play around with radial location - put close enough together that gaps start to overlap

//use integral (add fdens)