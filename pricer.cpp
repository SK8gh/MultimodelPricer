#include <vector>
#include <math.h>
#include <iostream>
#include <list>

#include "pricer.h"
#include "random.h"
#include "random.cpp"

//Export library
#include <fstream>

using namespace std;

// Helpful functions

double mean(const std::vector<double>& vect){
    double mean_val = 0;
    for (int j=0; j<vect.size(); j++){
        mean_val += vect[j] / vect.size();
    }
    return(mean_val);
}

double standard_dev(const std::vector<double>& vect){
    double mean_val = mean(vect);

    double var = 0;
    for (int j=0; j<vect.size(); j++){
        var += pow((vect[j] - mean_val),2) / vect.size();
    }
    return(sqrt(var));
}


// Option payoffs

double call_value(const double& stock_T, const double& strike){
    if (stock_T > strike){
        return(stock_T - strike);}
    else{return 0;}
}

double put_value(const double& stock_T, const double& strike){
    if (stock_T < strike){
        return(strike - stock_T);}
    else{return 0;}
}

double call_spread_value(const double& stock_T, const double& strike1, const double& strike2){
    return(call_value(stock_T, strike1) - call_value(stock_T, strike2));
}

double digit_value(const double& stock_T, const double& strike, const double& epsilon, const double& number_calls){
    return(number_calls * call_value(stock_T, strike - epsilon) - number_calls * call_value(stock_T, strike));
}

// Implementing the BlackScholes class

BlackScholes::BlackScholes(const double& mu_arg, const double& sigma_arg) : model(mu_arg), mu_m(mu_arg), sigma_m(sigma_arg){};

// Magicien::Magicien() : Personnage(), m_mana(100){};

double BlackScholes::drift_term(const double& time_step, const double& asset_price) const{
    double term = mu_m * asset_price * time_step;
    return(term);
}

double BlackScholes::vol_term(const double& time_step, const double& asset_price){
    double Z = random::normal();
    double term = asset_price * sigma_m * sqrt(time_step) * Z;
    return(term);
}


// Implementing the SABR class

SABR::SABR(const double& beta_arg, const double& alpha_arg, const double& rho_arg, const double& vol_spot_arg, 
const double& rate_arg) : model(rate_arg), beta_m(beta_arg), alpha_m(alpha_arg), rho_m(rho_arg), 
vol_m(vol_spot_arg){};

void SABR::update_vol(const double& new_vol){
    SABR::vol_m = new_vol;
}

double SABR::drift_term(const double& time_step, const double& asset_price) const{
    return(0);
}

double SABR::vol_term(const double& time_step, const double& asset_price){
    // Compute the stock volatility term, using the value of the vol at time t, that we will update afterwards
    double W = random::normal();
    double Z = rho_m * W + (sqrt(1 - pow(rho_m, 2)) * random::normal());

    double term = vol_m * pow(asset_price, beta_m) * sqrt(time_step) * W;

    // Updating the volatility process value
    double new_vol = vol_m + (alpha_m * vol_m * sqrt(time_step) * Z);
    this->update_vol(new_vol);
    
    return(term);
}


// Implementing the Heston class

Heston::Heston(const double& mu_arg, const double& kappa_arg, const double& theta_arg, 
    const double& sigma_arg, const double& rho_arg) : model(mu_arg), mu_m(mu_arg), kappa_m(kappa_arg),
    theta_m(theta_arg), sigma_m(sigma_arg), rho_m(rho_arg){};

void Heston::update_vol(const double& new_vol){
    Heston::vol = new_vol;
}

double Heston::drift_term(const double& time_step, const double& asset_price) const{
    double term = mu_m * asset_price * time_step;
    return(term);
}

double Heston::vol_term(const double& time_step, const double& asset_price){
    double W = random::normal();
    double Z = rho_m * W + (sqrt(1 - pow(rho_m, 2)) * random::normal());

    double term = asset_price * sqrt(time_step) * W * sqrt(vol);

    // Volatility update
    double new_vol = vol + kappa_m * (theta_m - vol) * time_step + sigma_m * sqrt(vol) * sqrt(time_step) * Z;
    this->update_vol(new_vol);

    // Return volatility term
    return(term);
}



// PathSimulator

std::vector<double> PathSimulator::Euler(model* model_ptr, const double& maturity, const double& spot, const double& vol_spot, 
    const double& model_steps){
        // Initializing the stock price vector
        std::vector<double> stock_price;
        stock_price.push_back(spot);
        double time_step = maturity / model_steps;

        for (int j=0; j<model_steps; j++){
            double old_price = stock_price.back();
            double next_price = old_price + model_ptr->drift_term(time_step, old_price) + model_ptr->vol_term(time_step, old_price);
            stock_price.push_back(next_price);
        }

        // Returning the stock price vector
        return(stock_price);
    }


// Option pricing

std::vector<double> PathSimulator::call(model* model_ptr, const double& maturity, const double& strike, 
const double& spot, const double& vol_spot, const double& model_steps, const double& number_simulations){
        std::vector<double> call_values;

        for (int j=0; j<number_simulations; j++){
            std::vector<double> path = PathSimulator::Euler(model_ptr, maturity, spot, vol_spot, model_steps);
            double call_val = call_value(path.back(), strike);
            call_values.push_back(call_val);
        }
        // rate = zero risk rate
        double call_mean = exp(- model_ptr->constant_rate * maturity) * mean(call_values);
        double call_std = standard_dev(call_values);

        std::vector<double> return_vector;
        return_vector.push_back(call_mean);
        return_vector.push_back(2 * call_std / sqrt(number_simulations));

        return(return_vector);
    }

std::vector<double> PathSimulator::put(model* model_ptr, const double& maturity, const double& strike, 
const double& spot, const double& vol_spot, const double& model_steps, const double& number_simulations){
        std::vector<double> put_values;

        for (int j=0; j<number_simulations; j++){
            std::vector<double> path = PathSimulator::Euler(model_ptr, maturity, spot, vol_spot, model_steps);
            double put_val = put_value(path.back(), strike);
            put_values.push_back(put_val);
        }
        // rate = zero risk rate
        double put_mean = exp(- model_ptr->constant_rate * maturity) * mean(put_values);
        double put_std = standard_dev(put_values);

        std::vector<double> return_vector;
        return_vector.push_back(put_mean);
        return_vector.push_back(2 * put_std / sqrt(number_simulations));

        return(return_vector);
    }

std::vector<double> PathSimulator::call_spread(model* model_ptr, const double& maturity, const double& strike1, 
const double& strike2, const double& spot, const double& vol_spot, const double& model_steps, 
const double& number_simulations){
        std::vector<double> call_spread_values;

        for (int j=0; j<number_simulations; j++){
            std::vector<double> path = PathSimulator::Euler(model_ptr, maturity, spot, vol_spot, model_steps);
            double call_spread_val = call_spread_value(path.back(), strike1, strike2);
            call_spread_values.push_back(call_spread_val);
        }
        // rate := zero risk rate
        double call_spread_mean = exp(- model_ptr->constant_rate * maturity) * mean(call_spread_values);
        double call_spread_std = standard_dev(call_spread_values);

        std::vector<double> return_vector;
        return_vector.push_back(call_spread_mean);
        return_vector.push_back(2 * call_spread_std / sqrt(number_simulations));

        return(return_vector);
    }

std::vector<double> PathSimulator::digitale(model* model_ptr, const double& maturity, const double& strike, 
const double& epsilon, const double& number_calls, const double& spot, const double& vol_spot, 
const double& model_steps, const double& number_simulations){
        std::vector<double> digit_values;

        for (int j=0; j<number_simulations; j++){
            std::vector<double> path = PathSimulator::Euler(model_ptr, maturity, spot, vol_spot, model_steps);
            double digit_val = digit_value(path.back(), strike, epsilon, number_calls);
            digit_values.push_back(digit_val);
        }
        // rate = zero risk rate
        double digit_mean = exp(- model_ptr->constant_rate * maturity) * mean(digit_values);
        double digit_std = standard_dev(digit_values);

        std::vector<double> return_vector;
        return_vector.push_back(digit_mean);
        return_vector.push_back(2 * digit_std / sqrt(number_simulations));

        return(return_vector);
    }


