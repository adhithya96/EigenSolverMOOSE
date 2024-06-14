#pragma once
#include "eigen-master/Eigen/Eigen"

class gauss_points
{
private:
    Eigen::VectorXd w, x;
public:
    Eigen::VectorXd get_weights();
    Eigen::VectorXd get_points();
    gauss_points(int order);
};