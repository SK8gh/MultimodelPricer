#pragma once
// Avoids double inclusions

class random
{
public:
    static double uniform();
    static double normal();
    // The static keyword let us call the methods without instanciating the class
};