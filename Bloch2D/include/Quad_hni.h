#pragma once

#include<iostream>
#include<fstream>
#include<string>
#include<sstream>
#include<cassert>
#include "eigen-master/Eigen/Eigen"
#include "Material.h"
#include "Quad.h"
#include "gauss_points.h"
#include "Levelset.h"

class Quad_hni : public Quad
{
private:
    int nnode, nelem, nen, ndof, order;
    Eigen::MatrixXd node;
    Eigen::MatrixXd elem;
    Eigen::MatrixXd Ke, Me;
    std::vector<int> cut_elem;
    
public:
    Quad_hni(const int nelex, const int neley, const double lx, const double ly){};
    void add_cut_elements(int elem);
    int  get_cut_elements(int pos);
    Eigen::VectorXd get_N(int i, double exi, double eta);
    Eigen::VectorXd get_dNdx(int i, int j, double exi, double eta);
    //local stiffness - affine curves
    Eigen::MatrixXd get_localstiffness(Eigen::MatrixXd D, Eigen::MatrixXd x);
    //local stiffness - parametric curves
    Eigen::MatrixXd get_localstiffness(Eigen::MatrixXd D, Eigen::MatrixXd node, Levelset *ls, 
        std::vector<double> t0, std::vector<double> t1, std::vector<double> p0, std::vector<double>p1, 
        double m0, double m1);
    //local mass - parametric curves    
    Eigen::MatrixXd get_localmass(double rho, Eigen::MatrixXd node, Levelset *ls, 
        std::vector<double> t0, std::vector<double> t1, std::vector<double> p0, std::vector<double>p1, 
        double m0, double m1);
    //local stiffness - affine curves
    Eigen::MatrixXd get_localmass(double rho, Eigen::MatrixXd x);
};
