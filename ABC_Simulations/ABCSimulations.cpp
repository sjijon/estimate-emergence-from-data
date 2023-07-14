#include <iostream>
#include <filesystem>
#include <cmath>
#include <numeric>
#include <algorithm>
#include <fstream>
#include <time.h>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>
#include <chrono>
using namespace std::chrono;

// command for compiling
// g++ -std=c++17 ABCSimulations.cpp `pkg-config --libs gsl` 

using namespace std;

////////////////////////////////////////
// parameters
////////////////////////////////////////

// secondary infections distribution - negative binomial
#define kappa 0.1                  // dispersion parameter
#define R0 2.5                       // mean number of secondary infections
#define p_neg kappa/(kappa+R0)      // success probability for the negative-binomial      

// Infectiousness distribution - gamma distribution (mean = 6 days, SD = 2.5 days)
#define shape_inf 6.6
#define scale_inf 0.833

// hospitalization time - gamma distribution
#define shape_detect 6.25
#define scale_detect 1.04

// hospitalization probability
#define p_detect 0.15

// number of observed cases in data
#define N 3072

// Repetitions
#define repeats 10000                             

// Random number generation with Mersenne Twister
gsl_rng * r = gsl_rng_alloc (gsl_rng_mt19937);

////////////////////////////////////////
// Main code
////////////////////////////////////////

int RUN(double, int);

int main(int argc,char *argv[])
{
    // start timer
    auto start = high_resolution_clock::now();

    // parameters to be read in
    double noDays = atof(argv[1]);                    // number of days the simulation needs to be run!
    //int r_ind = 0;                                    // repetition index
    int success = 0;                                  // epidemic counter

    // create folder where simulations are stored
    namespace fs = std::filesystem;
    fs::create_directories("./Days_" + std::to_string((int)noDays));

    // set random seed
    srand(time(0));

    // stochastic simulation loop
    while (success < repeats)
    {
        gsl_rng_set(r,rand());                           // setting the seed
        int curSuc = success;                           // indexing variable for the files
        success += RUN(noDays, curSuc);           // Run stochastic simulation
        //r_ind++;
    }
    
    // end timer
    auto stop = high_resolution_clock::now();
    auto duration = duration_cast<seconds>(stop - start);
    cout << noDays;
    cout << "\n";
    cout << duration.count() << endl;
    return(0);
}

////////////////////////////////////////
// Stochastic simulation
////////////////////////////////////////
int RUN(double daysToEnd, int curSuc)
{
    double t = 0.;                            // start day (t = 0 = day one)
    double dt = 0.1;                          // time step !!! WARNING: Only values of format 10^(-k) allowed
    int length = (int)daysToEnd/dt;           // array length
    int return_value = 0;
    int index = 0;
    
    int I_work[length+1] = {0};                 // number of infecteds working vector (continuous time approximation)
    int I[(int)daysToEnd] = {0};                // number of infecteds at a certain time (days!)
    int D[(int)daysToEnd] = {0};                // number of detected at a certain time (days!)
    
    // Initialization = 1 infected individual at time 0, none detected
    I_work[0] = 1;                  
    I[0] = 1;                       

    // Random distributions necessary throughout
    unsigned int gsl_ran_negative_binomial(const gsl_rng * rrest, double p, double n);  // p = success probability, n = sample size
    double gsl_ran_gamma(const gsl_rng * rrest, double a, double b);        // a = shape, b = scale
    unsigned int gsl_ran_binomial(const gsl_rng * rbin, double p, unsigned int n); // p = prob, n = size
    
    // Loop through time
    int work_ind = 0;
    while (t <= daysToEnd - 1.)
    {   
        // check for detection of individuals added at time t
        int n_samp = gsl_ran_binomial(r, p_detect, I_work[work_ind]);
        
        // loop over detected individuals and add their detection days
        for (int i = 1; i <= n_samp; i++)
        {
            double samp_time = t + gsl_ran_gamma(r,shape_detect,scale_detect);
            int D_ind = ceilf(samp_time);
            if (D_ind < (int)daysToEnd)             // note that index starts at 0!
            {
                D[D_ind]++;
            }  
        }
        
        // no. of offspring of individuals added at time t
        int offspring = 0;
        for (int i = 1; i <= I_work[work_ind]; i++)
        {
            offspring += gsl_ran_negative_binomial(r,p_neg,kappa);
        }
        
        // loop over offsprings from individuals that were added at time t and add their infection times
        for (int i = 1; i <= offspring; i++)
        {   
            double inf_time = t+gsl_ran_gamma(r,shape_inf,scale_inf);
            // IMPORTANT THAT dt = 10^(-k)! otherwise adapt this rounding formula!
            int inf_time_ind = roundf(inf_time/dt);
            int I_ind = ceilf(inf_time);
            
            if (inf_time_ind < length)
            {
                I_work[inf_time_ind]++;
            }
            
            if (I_ind < (int)daysToEnd)
            {
                I[I_ind]++;
            }            
        }
        
        t += dt;
        work_ind++;

        // stopping condition if epidemics grow too fast (number of detection > 10*observed detections)
        if (accumulate(D,D+(int)daysToEnd,0) > 10*N)
            break;
    }
       
    // If epidemic happened (at least 1000 infections) increase success counter and save number of detected per day
    if (accumulate(I,I+(int)daysToEnd,0) > 1000)
    {
        // change success counter
        return_value = 1;

        // initialize file    
        ofstream file ("./Days_" + std::to_string((int)daysToEnd) + "/ABCSim_no_" + std::to_string(curSuc) + ".txt", ios::app);

        // write into file
        for (int i=0; i<(int)daysToEnd;i++)
        {
            file << i+1;
            file << ",";
            file << I[i];
            file << ",";
            file << D[i];
            file << "\n";
        }

        // close file
        file.close();
    }

    return(return_value);
}





