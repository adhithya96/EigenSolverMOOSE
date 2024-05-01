#pragma once

#include "Material.h"
#include "Quad.h"

class QuadBilinear : Material, Quad
{

public:
    QuadBilinear(const int nelex, const int neley, const double lx, const double ly);

    int get_bnodes(int nodenum, std::string type);
    int get_bnodes_size(std::string type);
    Eigen::VectorXd get_coordinates(int  nodenum);
    Eigen::VectorXd get_connectivity(int elenum);
    int get_nnode(){return nnode;};
    int get_nelem(){return nelem;};
    Eigen::VectorXd get_boundarynodes();
    Eigen::MatrixXd get_localstiffness(Eigen::MatrixXd D, Eigen::MatrixXd x);
    double get_N(int i, double exi, double eta);
    double get_dNdx(int i, int j, double exi, double eta);
    Eigen::MatrixXd get_localmass(double rho, Eigen::MatrixXd node);
    Eigen::MatrixXd get_constraintmatrix(double kx, double ky, double a, int nnodex, int nnodey);
    Eigen::MatrixXd get_jacobianmat(Eigen::MatrixXd x, double exi, double eta);
};
