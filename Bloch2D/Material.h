#pragma once

class Material
{
private:
    double E, nu, rho;

public:
    Material(){};
    double get_E(){return E;};
    double get_nu(){return nu;};
    double get_rho(){return rho;};
    double get_wavenum(int choice);
    void set_wavenum(double kx, double ky);
    double get_latticevec(){return a;};
    Material(double E, double nu, double rho, double kx, double ky, double a);
    Eigen::MatrixXd get_Dmatrix();
    double kx, ky, a;
};

void  Material::set_wavenum(double kx, double ky)
{
    this->kx = kx;
    this->ky = ky;
}

double Material::get_wavenum(int choice)
{
    if(choice == 1)
        return kx;
    else 
        return ky;
}

Material::Material(double E, double nu, double rho, double kx, double ky, double a)
{
    this->E = E;
    this->nu =  nu;
    this->rho = rho;
    this->kx = kx;
    this->ky = ky;
    this->a = a;
}

Eigen::MatrixXd Material::get_Dmatrix()
{
    //Lame parameters
    double lambda = nu * E / ((1 + nu) * ( 1 - 2*  nu));
    double mu = E / (2 * (1 + nu));

    Eigen::MatrixXd D(3, 3);
    D << lambda + 2*mu, lambda, 0, 
        lambda, lambda + 2 * mu, 0, 
        0, 0, mu;

    return D;
}