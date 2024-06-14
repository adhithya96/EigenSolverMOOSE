#pragma once
#include <math.h>
#include "eigen-master/Eigen/Eigen"
#include<iostream>

class Circle
{
private:
    double a, b, r;
public: 
    Circle(double  a, double b, double r);
    double phi(double x, double y);
    double dphidx(double x, double y);
    double dphidy(double x, double y);
    Eigen::VectorXd dfdm(std::vector<double> t0, std::vector<double> t1, Eigen::VectorXd m);
    Eigen::MatrixXd ddfdm(std::vector<double> t0, std::vector<double> t1);
    double norm(std::vector<double> vec);
};