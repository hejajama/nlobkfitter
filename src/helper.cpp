
/*
 * Cache BK solution s.t. the code does not need to solve it again
 * if only sigma_r changes
 *
 * The idea is that the previous bk solution is always saved on the disk
 *
 * NOTE: THIS IS VERY DANGEROUS PIECE OF CODE
 * If MNINUIT is not complied in a single-thread mode, crazy things happen!!!
 */

#include "helper.hpp"
#include <vector>
#include <iostream>
#include <Minuit2/MnUserParameterState.h>
#include <cstdlib>
#include <cmath>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_min.h>


using namespace std;

using namespace ROOT::Minuit2;



// Optimize sigma02 (the coefficient which multiplies SigmaComputer.SigmarLOmass
// Returns optimal sigma02

struct sigma02helper
{
    std::vector<double> expdata;
    std::vector<double> experr;
    std::vector<double> thdata;
    
};

double minimiser_helper_chisqr(double sigma02, void* p)
{
    sigma02helper *par = (sigma02helper*)p;
    
    double sum=0;
    for (unsigned int i=0; i< par->expdata.size(); i++)
    {
        sum += std::pow((sigma02*par->thdata[i]-par->expdata[i])/par->experr[i], 2.0)/par->expdata.size();
    }
    //cout << sigma02 << " " << sum << endl;
    return sum;
}

std::vector<double> FindOptimalSigma02(std::vector<double> expdata, std::vector<double> experr,
                          std::vector<double> thdata)
{
    sigma02helper par;
    par.expdata = expdata;
    par.experr = experr;
    par.thdata = thdata;
    
    if (expdata.size() != experr.size() or expdata.size() != thdata.size())
    {
        cerr << "Optimizer did not get vectors with same size???" << endl;
        exit(1);
    }
    
    int status;
    int iter = 0, max_iter = 100;
    const gsl_min_fminimizer_type *T;
    gsl_min_fminimizer *s;
    double min = 30.0;
    double lower = 1; double upper = 150;;
    gsl_function F;
    
    F.function = &minimiser_helper_chisqr;
    F.params = &par;
    
    T = gsl_min_fminimizer_brent;
    s = gsl_min_fminimizer_alloc (T);
    gsl_min_fminimizer_set (s, &F, min, lower, upper);
    
    
    do
    {
        iter++;
        status = gsl_min_fminimizer_iterate (s);
        
        min = gsl_min_fminimizer_x_minimum (s);
        lower = gsl_min_fminimizer_x_lower (s);
        upper = gsl_min_fminimizer_x_upper (s);
        
        status = gsl_min_test_interval (lower, upper, 0, 0.0001);
    }
    while (status == GSL_CONTINUE && iter < max_iter);
    
    gsl_min_fminimizer_free (s);
    
    //cout << "Minimizer ended, iterations " << iter << " minimum " << min << endl;
    
    std::vector<double> result;
    result.push_back(min);
    result.push_back(minimiser_helper_chisqr(min, &par));
    return result;
}

