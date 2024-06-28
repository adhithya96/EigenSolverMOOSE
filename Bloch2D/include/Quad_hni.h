#pragma once

#include<iostream>
#include<fstream>
#include<string>
#include<sstream>
#include<cassert>
#include "eigen-master/Eigen/Eigen"
#include "Material.h"
#include "Quad.h"
#include "Quad2.h"
#include "Quad3.h"
#include "Quad4.h"
#include "Quad5.h"
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
    int get_nen();
    int get_nnode();
    int get_nelem();
    int get_ndof();
    int get_order();
    Eigen::VectorXd get_coordinates(int nodenum);
    Eigen::VectorXd get_connectivity(int elenum);
    Quad_hni(const int nelex, const int neley, const double lx, const double ly, int order);
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
    template <typename T>
    void copymesh(T* mesh);

};


template <typename T>
void Quad_hni::copymesh(T* mesh)
{
    nen = mesh->get_nen();
    nnode = mesh->get_nnode();
    nelem = mesh->get_nelem();
    ndof = mesh->get_ndof();
    order = mesh->get_order();

    node = Eigen::MatrixXd::Zero(nnode, 2);
    for(int i = 0; i < nnode; i++)
        for(int j = 0; j < 2; j++)
            node(i, j) = mesh->get_coordinates(i)(j);

    elem = Eigen::MatrixXd::Zero(nelem, nen);
    for(int i = 0; i < nelem; i++)
        for(int j = 0; j < nen; j++)
            elem(i, j) = mesh->get_connectivity(i)(j);
}