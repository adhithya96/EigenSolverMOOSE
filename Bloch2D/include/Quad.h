#pragma once

#include<iostream>
#include<fstream>
#include<string>
#include<sstream>
#include<cassert>
#include "eigen-master/Eigen/Eigen"
#include "Material.h"

class Quad : public Material
{
private:
    int nnode, nelem, nen, ndof, order;
    Eigen::MatrixXd node;
    Eigen::MatrixXd elem;
    Eigen::MatrixXd Ke, Me;
    std::vector<int> cut_elem;
    
public:
    Quad(){};
    Quad(const int nelex, const int neley, const double lx, const double ly);
    Quad(std::string fname);
    int get_nen();
    Eigen::VectorXd get_coordinates(int  nodenum);
    Eigen::VectorXd get_connectivity(int elenum);
    int get_order();
    int get_nnode();
    int get_nelem();
    int get_ndof();
    Eigen::VectorXd get_boundarynodes();
    Eigen::MatrixXd get_localstiffness(Eigen::MatrixXd D, Eigen::MatrixXd x);
    double get_N(int i, double exi, double eta);
    double get_dNdx(int i, int j, double exi, double eta);
    Eigen::MatrixXd get_localmass(double rho, Eigen::MatrixXd node);
    Eigen::MatrixXd get_constraintmatrix(double kx, double ky, double a, int nnodex, int nnodey);
    Eigen::MatrixXd get_jacobianmat(Eigen::MatrixXd x, double exi, double eta);
    void add_cut_elements(int elem);
    int  get_cut_elements(int pos);
};

