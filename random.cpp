#include "random.h"
#include <math.h>
#include <random>

double random::uniform()
{   
    double u = ((double)(rand()) + 1.) / ((double)(RAND_MAX) + 1.);
    return(u);
}

double random::normal()
{
    double u1 = uniform();
    double u2 = uniform();
    double pi = 3.14159265;

    double norm = cos(2*pi*u1)*(sqrt(-2*log(u2)));
    return(norm);
};