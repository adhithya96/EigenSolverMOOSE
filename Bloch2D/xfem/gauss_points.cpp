#include "../include/gauss_points.h"

gauss_points::gauss_points(int order)
{
    w = Eigen::VectorXd::Zero(order);
    x = Eigen::VectorXd::Zero(order);
    if(order == 1)
    {
        x(0) = 0;
        w(0) = 2;
    }
    if(order == 2)
    {
        x(0) = 0.57735;
        x(1) = -0.57735;

        w(0) = 1;
        w(1) = 1;
    }
    if(order == 3)
    {
        x(0) = -0.77459;
        x(1) = 0;
        x(2) = 0.77459;

        w(0) =  0.55556;
        w(1) =  0.88889;
        w(2) =  0.55556;
    }
    if(order ==  4)
    {
        x(0) = -0.86114;
        x(1) = -0.33998;
        x(2) = 0.33998;
        x(3) = 0.86114;

        w(0) = 0.65214;
        w(1) = 0.34785;
        w(2) = 0.34785;
        w(3) = 0.65214;
    }
    if(order == 5)
    {
        x(0) = -0.90618;
        x(1) = -0.53847;
        x(2) = 0;
        x(3) = 0.53847;
        x(4) = 0.90618;  

        w(0) = 0.236927;
        w(1) = 0.478629;
        w(2) = 0.568889;
        w(3) = 0.478629;
        w(4) = 0.236927;
    }
}

Eigen::VectorXd gauss_points::get_points()
{
    return x;
}

Eigen::VectorXd gauss_points::get_weights()
{
    return w;
}