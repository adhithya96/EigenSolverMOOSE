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

Eigen::VectorXd Circle::dfdm(std::vector<double> t0, std::vector<double> t1, 
                            Eigen::VectorXd m, std::vector<double> p0, std::vector<double> p1)
{
    Eigen::VectorXd ans(2);
    
    Rosetta::LegendrePolynomial gp(4);
    // std::cout << "roots " << std::endl;
    // for(int i = 1; i <= 4; i++)
    //     std::cout << gp.root(i) << std::endl;
    // std::cout << "weights" << std::endl;
    // for(int i = 1; i <= 4; i++)
    //     std::cout << gp.weight(i) << std::endl;
    for(int i = 1; i <= 4; i++)
    {
        double t = (1 - gp.root(i)) * 0 / 2 + (1 + gp.root(i)) * 1 / 2;
        //dfdm0
        Eigen::VectorXd dcxdm = Levelset::dcxdm(t, t0, t1);
        Eigen::VectorXd dcydm = Levelset::dcydm(t, t0, t1);
        Eigen::VectorXd c = Levelset::c(t, t0, t1, p0, p1, m(0), m(1));
    
        double phi = Circle::phi(c(0), c(1));

        ans(0) += 2 * gp.weight(i) * phi * (1 / c.norm()) * ((c(0) - a) * dcxdm(0) + (c(1) - b) * dcydm(0));
        ans(1) += 2 * gp.weight(i) * phi * (1 / c.norm()) * ((c(0) - a) * dcxdm(1) + (c(1) - b) * dcydm(1)); 
    }
    std::cout << ans << std::endl;
    
    return ans * 0.5;
}

Eigen::MatrixXd Circle::ddfdm(std::vector<double> t0, std::vector<double> t1, 
                            Eigen::VectorXd m, std::vector<double> p0, std::vector<double> p1)
{
    Eigen::MatrixXd ans =  Eigen::MatrixXd::Zero(2, 2);

    Rosetta::LegendrePolynomial gp(2);

    for(int i = 1; i <= 2; i++)
    {
        //d2F(m)/dm0
        double t = (1 - gp.root(i)) * 0 / 2 + (1 + gp.root(i)) * 1 / 2;
        
        //dA/dm0
        Eigen::VectorXd dcxdm = Levelset::dcxdm(t, t0, t1);
        std::cout << dcxdm << std::endl;
        Eigen::VectorXd dcydm = Levelset::dcydm(t, t0, t1);
        std::cout << dcydm << std::endl;
        Eigen::VectorXd c = Levelset::c(t, t0, t1, p0, p1, m(0), m(1));
        std::cout << c << std::endl;

        double phi = Circle::phi(c(0), c(1));

        //A
        ans(0, 0) += gp.weight(i) * (2 * pow((dcxdm(0) * (c(0) - a) + dcydm(0) * (c(1) - b)), 2) / (pow(c(0) - a, 2) + pow((c(1) - b), 2)) - 
                    2 * phi * pow((c(0) - a) * dcxdm(0) +  (c(1) - b) * dcydm(0), 2) / pow(pow(c(0) - a, 2) + pow((c(1) - b), 2), 1.5) + 
                    2 * phi * (pow(dcxdm(0), 2) + pow(dcydm(0), 2)) / sqrt((c(0) - a, 2) + pow(c(1) - b, 2))); 
        std::cout << ans(0, 0) << std::endl;
        
        ans(1, 1) += gp.weight(i) * (2 * pow((dcxdm(1) * (c(0) - a) + dcydm(1) * (c(1) - b)), 2) / (pow(c(0) - a, 2) + pow((c(1) - b), 2)) - 
                2 * phi * pow((c(0) - a) * dcxdm(1) +  (c(1) - b) * dcydm(1), 2) / pow(pow(c(0) - a, 2) + pow((c(1) - b), 2), 1.5) + 
                2 * phi * (pow(dcxdm(1), 2) + pow(dcydm(1), 2)) / sqrt((c(0) - a, 2) + pow(c(1) - b, 2))); 
        std::cout << ans(1, 1) << std::endl;

        ans(0, 1) += gp.weight(i) * (2 * ((dcxdm(0) + dcxdm(1)) * (c(0) - a) + (dcydm(0) + dcydm(1)) * (c(1) - b)) / (pow(c(0) - a, 2) + pow(c(1) - b, 2)) - 
                2 * phi * ((c(0) - a) * (dcxdm(0) + dcxdm(1)) +  (c(1) - b) * (dcydm(0) + dcydm(1))) / pow(pow(c(0) - a, 2) + pow((c(1) - b), 2), 1.5) + 
                2 * phi * (dcxdm(0) * dcxdm(1) + dcydm(0) * dcydm(1)) / sqrt((c(0) - a, 2) + pow(c(1) - b, 2))); 
        std::cout << ans(0, 1) << std::endl;
    }

    ans(1, 0) = ans(0, 1);
    
    std::cout << ans << std::endl;
    return ans * 0.5;
}