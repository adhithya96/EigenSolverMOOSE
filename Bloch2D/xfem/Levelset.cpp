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

double Levelset::norm(Eigen::VectorXd vec)
{
    return sqrt(pow(vec(0), 2) + pow(vec(1), 2));
}

double Levelset::evaluate_area_affine(Eigen::MatrixXd node)
{
    gauss_points gp = gauss_points(1);
    Eigen::VectorXd w = gp.get_weights();
    std::cout << w << std::endl;
    Eigen::VectorXd x = gp.get_points();
    std::cout << x << std::endl;
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
        for(int j = 0; j < w.size(); j++)
            area += (1.0 / 2.0) * b * w(j) * jacob;
        
        std::cout << area << std::endl;
    }

    return area;
}

Eigen::VectorXd Levelset::c(double t, int pos)
{
    Eigen::VectorXd ans(2);
    for(int i = 0; i < 2; i++)
        ans(i) = (2*pow(t, 3) - 3*pow(t, 2) + 1) * p0[pos][i] + (pow(t, 3) - 
        2 * pow(t, 2) + t) * m0[pos] * t0[pos][i] + (-2 * pow(t, 3) + 3 * pow(t, 2)) * p1[pos][i] +
        (pow(t, 3) - pow(t, 2)) * m1[pos] * t1[pos][i];
    
    return ans;
}

Eigen::VectorXd Levelset::dcdt(double t, int pos)
{
    Eigen::VectorXd ans(2);
    for(int i = 0; i < 2; i++)
        ans(i) = (6 * pow(t, 2) - 6 * t) * p0[pos][i] + (3 * pow(t, 2) - 4 * t + 1) * m0[pos] * t0[pos][i]
        + (-6 * pow(t, 2) + 6 * t) * p1[pos][i] + (3 * pow(t, 2) - 2 * t) * m1[pos] * t1[pos][i];
    
    return ans;
}

double Levelset::evaluate_area_parametric(int pos)
{
    gauss_points gp = gauss_points(2);
    Eigen::VectorXd x = gp.get_points();
    Eigen::VectorXd w = gp.get_weights();
    double area = 0; 
    for(int i = 0; i < 2; i++)
    {
        Eigen::VectorXd c = Levelset::c(x(i), pos);
        Eigen::VectorXd dcdt = Levelset::dcdt(x(i), pos);

        area += (c(0) * dcdt(1) - c(1) * dcdt(0)) * w(i);
    }

    return area;

}