#pragma once

#include "Quad.h"

class Quad2 : public Material
{
protected:
    int nnode, nelem, nen, ndof, ngp;
    Eigen::MatrixXd node;
    Eigen::MatrixXd elem;
    Eigen::MatrixXd Ke, Me;
    
public:
    Quad2(const int nelex, const int neley, const double lx, const double ly);
    Quad2(std::string fname);
    int get_nnode();
    int get_nelem();
    int get_ndof();
    int get_nen();
    int get_ngp();
    Eigen::VectorXd get_coordinates(int  nodenum);
    Eigen::VectorXd get_connectivity(int elenum);
    Eigen::VectorXd RefElementCoord(int nodenum);
    Eigen::VectorXd GaussIntegrationWeights(int nodenum);
    Eigen::VectorXd get_boundarynodes();
    Eigen::MatrixXd get_localstiffness(Eigen::MatrixXd D, Eigen::MatrixXd x);
    double get_N(int i, double exi, double eta);
    double get_dNdx(int i, int j, double exi, double eta);
    Eigen::MatrixXd get_localmass(double rho, Eigen::MatrixXd node);
    Eigen::MatrixXd get_constraintmatrix(double kx, double ky, double a, int nnodex, int nnodey);
    Eigen::MatrixXd get_jacobianmat(Eigen::MatrixXd x, double exi, double eta);
};
