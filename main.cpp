#include <iostream>

#include "pricer.h"
#include "pricer.cpp"

#include <vector>
#include <list>
#include <math.h>

//Export library
#include <fstream>

using namespace std;

void save_to_csv(std::string path, std::vector<double> vect){
    std::ofstream file(path);

    for (int j=0; j<vect.size(); j++){
        file << vect[j] << std::endl;
    }

    file.close();
}

void print_payoff(std::string option, std::vector<double> option_summary, double strike){
    std::cout << option << " strike = " << strike << ". value = " << option_summary[0] << " +/- " << option_summary[1] << std::endl;
}

int main(){

    double maturity(1), spot(100), vol_spot(0.3), number_steps(100), number_MC(100000);

    // SABR
    double alpha(.3), beta(1);

    // SABR, BS, HESTON
    double mu(0.05), kappa(1.22), theta(0.3), sigma(0.3), rho(-0.66);

    Heston HestonModel(mu, kappa, theta, sigma, rho);
    SABR SABRmodel(beta, alpha, rho, vol_spot, mu);
    BlackScholes BSmodel(mu, sigma);

    // Digitales
    double epsilon = 1;
    double number_calls = 10;

    double strike = 100;
    double strike2 = 110;


    // Options 
    // Calls
    std::vector<double> call_BS = PathSimulator::call(&BSmodel, maturity, strike, spot, vol_spot, number_steps, number_MC);
    std::vector<double> call_SABR = PathSimulator::call(&SABRmodel, maturity, strike, spot, vol_spot, number_steps, number_MC);
    std::vector<double> call_Heston = PathSimulator::call(&HestonModel, maturity, strike, spot, vol_spot, number_steps, number_MC);

    print_payoff("Black Scholes Call", call_BS, strike);
    print_payoff("SABR Call", call_SABR, strike);
    print_payoff("Heston Call", call_Heston, strike);
    std::cout << "\n";


    // Puts
    std::vector<double> put_BS = PathSimulator::put(&BSmodel, maturity, strike, spot, vol_spot, number_steps, number_MC);
    std::vector<double> put_Heston = PathSimulator::put(&HestonModel, maturity, strike, spot, vol_spot, number_steps, number_MC);

    print_payoff("Black Scholes Put", put_BS, strike);
    print_payoff("Heston Put", put_Heston, strike);
    std::cout << "\n";


    // Call spread
    std::vector<double> call_spread_BS = PathSimulator::call_spread(&BSmodel, maturity, strike, strike2, spot, vol_spot, number_steps, number_MC);
    std::vector<double> call_spread_Heston = PathSimulator::call_spread(&HestonModel, maturity, strike, strike2, spot, vol_spot, number_steps, number_MC);

    print_payoff("Black Scholes Call Spread", call_spread_BS, strike);
    print_payoff("Heston Call Spread", call_spread_Heston, strike);
    std::cout << "\n";


    // Digitale
    std::vector<double> digitale_BS = PathSimulator::digitale(&BSmodel, maturity, strike, epsilon, number_calls, spot, vol_spot, number_steps, number_MC);
    std::vector<double> digitale_Heston = PathSimulator::digitale(&HestonModel, maturity, strike, epsilon, number_calls, spot, vol_spot, number_steps, number_MC);

    print_payoff("Black Scholes Digitale", digitale_BS, strike);
    print_payoff("Heston Digitale", digitale_Heston, strike);
    std::cout << "\n";
}
