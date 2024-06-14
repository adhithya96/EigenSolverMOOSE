#include "../include/Circle.h"


Circle::Circle(double a, double b, double r)
{
    this->a = a;
    this->b = b;
    this->r = r;
}

double Circle::phi(double x, double y)
{
    return sqrt(pow((x - a), 2) + pow((y - b),2)) - r;
}

double Circle::dphidx(double x, double y)
{
    return (x - a) / sqrt(pow(x - a, 2) + pow(y - b, 2));
}

double Circle::dphidy(double x, double y)
{
    return (y - b) / sqrt(pow(x - a, 2) + pow(y - b, 2));
}

double Circle::norm(std::vector<double> vec)
{
    return sqrt(pow(vec[0], 2) + pow(vec[1], 2));
}

Eigen::VectorXd Circle::dfdm(std::vector<double> t0, std::vector<double> t1, Eigen::VectorXd m)
{
    Eigen::VectorXd ans(2);
    ans(0) = (12.0/35.0) * pow(norm(t0), 2) * 2 * m(0);
    ans(1) = (1.0/105.0) * pow(norm(t1), 2) * 2 * m(1);

    std::cout << ans << std::endl;
    return ans;
}

Eigen::MatrixXd Circle::ddfdm(std::vector<double> t0, std::vector<double> t1)
{
    Eigen::MatrixXd ans =  Eigen::MatrixXd::Zero(2, 2);
    ans(0, 0) = (12.0/35.0) * pow(norm(t0), 2) * 2;
    ans(1, 1) = (1.0/105.0) * pow(norm(t1), 2) * 2;

    std::cout << ans << std::endl;
    return ans;
}