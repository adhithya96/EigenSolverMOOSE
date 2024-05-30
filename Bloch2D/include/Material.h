#pragma once

#include "eigen-master/Eigen/Eigen"

class Material
{
private:
    double E, nu, rho;

public:
    Material(){};
    double get_E();
    double get_nu();
    double get_rho();
    double get_wavenum(int choice);
    void set_wavenum(double kx, double ky);
    double get_latticevec();
    Material(double E, double nu, double rho, double kx, double ky, double a);
    Eigen::MatrixXd get_Dmatrix();
    double kx, ky, a;
};
