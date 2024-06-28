#include "../include/Levelset.h"
#include "../include/gauss_points.h"

Levelset::Levelset()
{
    
}

std::vector<double> Levelset::get_p0(int pos)
{
    return p0[pos];
}

std::vector<double> Levelset::get_p1(int pos)
{
    return p1[pos];
}

std::vector<double> Levelset::get_t0(int pos)
{
    return t0[pos];
}

std::vector<double> Levelset::get_t1(int pos)
{
    return t1[pos];
}

double Levelset::get_m0(int pos)
{
    return m0[pos];
}

double Levelset::get_m1(int pos)
{
    return m1[pos];
}

double Levelset::norm(Eigen::VectorXd vec)
{
    return sqrt(pow(vec(0), 2) + pow(vec(1), 2));
}

double Levelset::evaluate_area_affine(Eigen::MatrixXd node)
{
    Rosetta::LegendrePolynomial gp(1);

    double area = 0;
    for(int i = 0; i < 4; i++)
    {
        double b;
        double jacob;
        if(i == 0)
        {
            b = -node(0, 1);
            jacob =  (node(1, 0) - node(0, 0)) / 2.0;
        }
        if(i == 1)
        {
            b = node(1, 0);
            jacob =  (node(2, 1) - node(1, 1)) / 2.0;
        }
        if(i == 2)
        {
            b = node(2, 1);
            jacob =  (node(2, 0) - node(3, 0)) / 2.0;
        }
        if(i == 3)
        {
            b = -node(3, 0);  
            jacob =  (node(3, 1) - node(0, 1)) / 2.0;
        }

        area = (1.0 / 2.0) * b * gp.weight(1) * jacob;
        
        std::cout << area << std::endl;
    }

    return area;
}

Eigen::VectorXd Levelset::c(double t, std::vector<double> t0, std::vector<double> t1, 
                            std::vector<double> p0, std::vector<double> p1, double m0, double m1)
{
    Eigen::VectorXd ans(2);
    for(int i = 0; i < 2; i++)
        ans(i) = (2*pow(t, 3) - 3*pow(t, 2) + 1) * p0[i] + (pow(t, 3) - 
        2 * pow(t, 2) + t) * m0 * t0[i] + (-2 * pow(t, 3) + 3 * pow(t, 2)) * p1[i] +
        (pow(t, 3) - pow(t, 2)) * m1 * t1[i];
    
    return ans;
}

Eigen::VectorXd Levelset::dcdt(double t, std::vector<double> t0, std::vector<double> t1, 
                            std::vector<double> p0, std::vector<double> p1, double m0, double m1)
{
    Eigen::VectorXd ans(2);
    for(int i = 0; i < 2; i++)
        ans(i) = (6 * pow(t, 2) - 6 * t) * p0[i] + (3 * pow(t, 2) - 4 * t + 1) * m0 * t0[i]
        + (-6 * pow(t, 2) + 6 * t) * p1[i] + (3 * pow(t, 2) - 2 * t) * m1 * t1[i];
    
    return ans;
}

Eigen::VectorXd Levelset::dcxdm(double t, std::vector<double> t0, std::vector<double> t1)
{
    Eigen::VectorXd ans(2);
    //dcxdm0
    ans(0) = (pow(t, 3) - 2 * pow(t, 2) + t) * t0[0];
    //dcydm1
    ans(1) = (pow(t, 3) - pow(t, 2)) * t1[0];

    return ans;
}



Eigen::VectorXd Levelset::dcydm(double t, std::vector<double> t0, std::vector<double> t1)
{
    Eigen::VectorXd ans(2);
    //dcydm0
    ans(0) = (pow(t, 3) - 2 * pow(t, 2) + t) * t0[1];
    //dcydm1
    ans(1) = (pow(t, 3) - pow(t, 2)) * t1[1];

    return ans;
}

double Levelset::evaluate_area_parametric(std::vector<double> t0, std::vector<double> t1, 
                            std::vector<double> p0, std::vector<double> p1, double m0, double m1)
{
    Rosetta::LegendrePolynomial gp(2);
    double area = 0; 
    for(int i = 1; i <= 2; i++)
    {
        double t = (1 - gp.root(i)) * 0 / 2 + (1 + gp.root(i)) * 1 / 2; 
        Eigen::VectorXd c = Levelset::c(t, t0, t1, p0, p1, m0, m1);
        Eigen::VectorXd dcdt = Levelset::dcdt(t, t1, t1, p0, p1, m0, m1);

        area += (c(0) * dcdt(1) - c(1) * dcdt(0)) * gp.weight(i);
    }

    return area;

}