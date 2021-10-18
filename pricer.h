#pragma once

#include <iostream>
#include <vector>
#include <string>

using namespace std;

class model{
public:

    // Constructor
    model(const double& rate_arg) : constant_rate(rate_arg){};

    // Virtual methods
    virtual double drift_term(const double& time_step, const double& asset_price) const = 0; 
    virtual double vol_term(const double& time_step, const double& asset_price) = 0;
    
    const double& constant_rate;    
};

class BlackScholes : public model{
public:
    // Constructeur avec paramètres
    BlackScholes(const double& mu_arg, const double& sigma_arg);

    double drift_term(const double& time_step, const double& asset_price) const; 
    double vol_term(const double& time_step, const double& asset_price);

private:
    // Making the default constructor private such that an empty BlackScholes object cannot be instanciated
    BlackScholes();

    // Model parameters
    const double mu_m;
    const double sigma_m;
};

class SABR : public model{
public:
    // Constructeur avec paramètres
    SABR(const double& beta_arg, const double& alpha_arg, const double& rho_arg, const double& vol_spot,
    const double& constant_rate);

    void update_vol(const double& new_vol);

    double drift_term(const double& time_step, const double& asset_price) const; 
    double vol_term(const double& time_step, const double& asset_price);

private:
    // Making the default constructor private such that an empty BlackScholes object cannot be instanciated
    SABR();

    // Model parameters
    const double beta_m;
    const double alpha_m;
    const double rho_m;

    // Non-const volatility value
    double vol_m;
};

class Heston : public model{
public:
    // Constructeur avec paramètres
    Heston(const double& mu_arg, const double& kappa_arg, const double& theta_arg, 
    const double& sigma_arg, const double& rho_arg);

    void update_vol(const double& new_vol);

    double drift_term(const double& time_step, const double& asset_price) const; 
    double vol_term(const double& time_step, const double& asset_price);

private:
    // Making the default constructor private such that an empty BlackScholes object cannot be instanciated
    Heston();

    // Model parameters
    const double mu_m;
    const double kappa_m;
    const double theta_m;
    const double sigma_m;
    const double rho_m;

    // Non-const volatility value
    double vol;
};


class PathSimulator{
public:
    // Euler method for a specified scheme
    static std::vector<double> Euler(model* model_ptr, const double& maturity, const double& spot, const double& vol_spot, 
    const double& model_steps);

    // Option pricing, provide a value estimate, as well as a confidence interval
    static std::vector<double> call(model* model_ptr, const double& maturity, const double& strike, const double& spot, const double& vol_spot, 
    const double& model_steps, const double& number_simulations);

    static std::vector<double> call_spread(model* model_ptr, const double& maturity, const double& strike1, 
    const double& strike2, const double& spot, const double& vol_spot, const double& model_steps, 
    const double& number_simulations);

    static std::vector<double> put(model* model_ptr, const double& maturity, const double& strike, const double& spot, const double& vol_spot, 
    const double& model_steps, const double& number_simulations);

    // On approxime le payoff d'une digitale de strike K par un panier de call options : avec k calls de 
    // strike K - epsilon et k calls de strike K
    static std::vector<double> digitale(model* model_ptr, const double& maturity, const double& strike, 
    const double& epsilon, const double& number_calls, const double& spot, const double& vol_spot, const double& model_steps, 
    const double& number_simulations);

private:
};




