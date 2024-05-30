#include "../include/matplotlibcpp.h"
#include<iostream>
#include<fstream>
#include<string>
#include<sstream>

namespace plt = matplotlibcpp;

int main()
{
    double E = 7.31 * pow(10, 10);
    double nu = 0.325;
    double rho = 2770;
    std::fstream file;
    file.open("dispersion_data", std::ios::in);
    std::string line;
    std::vector<double>kfem, wfem;
    while(file.is_open() && !file.eof())
    {
        std::getline(file, line);
        std::stringstream s(line);
        double temp1, temp2;
        s >> temp1 >> temp2;
        kfem.push_back(temp1);
        wfem.push_back(temp2);
    }

    
    plt::figure(); 

    plt::ylim(0, 6);


    std::vector<double> ka, wa;
    double lambda = nu * E / ((1 + nu) * ( 1 - 2*  nu));
    double mu = E / (2 * (1 + nu));

    double c = sqrt((mu)/rho);
    
    double n = 0;
    double a = 0.04;
    double m;
    for(int i = 0; i < 3; i++)
    {
        m = i;
        for(int iter = 0; iter <= 20; iter++)
        {
            double ky = iter * 0.05 * M_PI / a;
            double temp = c * sqrt(pow((n * M_PI * 2 / a),2) + pow((ky + (m * M_PI * 2) / a),2));
            wa.push_back(temp * a / sqrt(E/(3*(1 - 2*nu)*rho)));
            ka.push_back(ky * a / M_PI);
            std::cout << ka[iter] << "  "  << wa[iter] << std::endl;
        }
    }

    plt::plot(ka, wa, {{"color", "black"}});
    
    
    ka.clear();
    wa.clear(); 
    
    n = 1;
    for(int i = 0; i < 3; i++)
    {
        m = i;
        for(int iter = 0; iter <= 20; iter++)
        {
            double ky = iter * 0.05 * M_PI / a;
            double temp = c * sqrt(pow((n * M_PI * 2 / a),2) + pow((ky + (m * M_PI * 2) / a),2));
            wa.push_back(temp * a / sqrt(E/(3*(1 - 2*nu)*rho)));
            ka.push_back(ky * a / M_PI);
            std::cout << ka[iter] << "  "  << wa[iter] << std::endl;
        }
    }
    plt::plot(ka, wa, {{"color", "black"}});

    ka.clear();
    wa.clear(); 
    
    c = sqrt((lambda + 2 * mu)/rho);
    n = 0;
    for(int i = 0; i < 3; i++)
    {
        m = i;
        for(int iter = 0; iter <= 20; iter++)
        {
            double ky = iter * 0.05 * M_PI / a;
            double temp = c * sqrt(pow((n * M_PI * 2 / a),2) + pow((ky + (m * M_PI * 2) / a),2));
            wa.push_back(temp * a / sqrt(E/(3*(1 - 2*nu)*rho)));
            ka.push_back(ky * a / M_PI);
            std::cout << ka[iter] << "  "  << wa[iter] << std::endl;
        }
    }

    plt::plot(ka, wa, {{"color", "black"}, {"label", "analytical"}});
    
    plt::scatter(kfem, wfem, {{"color", "red"}, {"label", "FEM"}});
    
    plt::legend();
    plt::title("Dispersion diagram");

    plt::xlabel("wave number");
    plt::ylabel("angular frequency");

    plt::savefig("Dispersion_Curve.png");

    return 0;
}