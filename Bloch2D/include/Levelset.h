#pragma once

#include "eigen-master/Eigen/Eigen"
#include "gauss_points.h"
#include <iostream>

class Levelset
{
private:
    std::vector<double> m0, m1;
    std::vector<std::vector<double>> p0, p1, t0, t1;

public:
    std::vector<double> get_p0(int pos);
    std::vector<double> get_p1(int pos);
    std::vector<double> get_t0(int pos);
    std::vector<double> get_t1(int pos);
    double get_m1(int pos);
    double get_m0(int pos);
    Levelset();
    Eigen::VectorXd c(double t, std::vector<double> t0, std::vector<double> t1, 
                            std::vector<double> p0, std::vector<double> p1, double m0, double m1);
    Eigen::VectorXd dcdt(double t, std::vector<double> t0, std::vector<double> t1, 
                            std::vector<double> p0, std::vector<double> p1, double m0, double m1);
    Eigen::VectorXd dcxdm(double t, std::vector<double> t0, std::vector<double> t1);
    Eigen::VectorXd dcydm(double t, std::vector<double> t0, std::vector<double> t1);
    template <typename T>
    int cut_element(Eigen::MatrixXd node, T *geom);
    // Find p0, p1, t0, t1
    template <typename T>
    void find_intersection(Eigen::MatrixXd node, T *geom);
    // Find m0, m1
    template <typename T>
    void find_tangent_mag(std::vector<double> p0, std::vector<double> p1, std::vector<double> t0, std::vector<double> t1, T *geom);
    template <typename T>
    double bisection(int dir, T *geom, std::vector<double> node1, std::vector<double> node2);
    double norm(Eigen::VectorXd vec);
    double evaluate_area_affine(Eigen::MatrixXd node);
    double evaluate_area_parametric(std::vector<double> t0, std::vector<double> t1, 
                            std::vector<double> p0, std::vector<double> p1, double m0, double m1);
};

template <typename T>
void Levelset::find_tangent_mag(std::vector<double> p0, std::vector<double> p1, std::vector<double> t0, std::vector<double> t1, T *geom)
{
    Eigen::VectorXd temp(2);
    temp << geom->dphidx(p0[0], p0[1]), geom->dphidy(p0[0], p0[1]);
    double dp0 = norm(temp);
    temp.setZero();
    temp << geom->dphidx(p1[0], p1[1]), geom->dphidy(p1[0], p1[1]);
    double dp1 = norm(temp);
    // initial guess
    Eigen::VectorXd m(2);
    m(0) = dp0;
    m(1) = dp1;
    std::cout << m << std::endl;
    int iter = 0;
    do
    {
        Eigen::VectorXd dfdm(2);
        Eigen::MatrixXd ddfdm(2, 2);
        
        dfdm = geom->dfdm(t0, t1, m, p0, p1);
        ddfdm = geom->ddfdm(t0, t1, m, p0, p1);

        //assert(ddfdm.determinant() > 0);
        m -= ddfdm.inverse() * dfdm;

        std::cout << m << std::endl;
        std::cout << norm(dfdm) << std::endl; 
        iter++;
    } while ((norm(geom->dfdm(t0, t1, m, p0, p1)) > pow(10, -3)) && (iter < 100));

    assert(norm(geom->dfdm(t0, t1, m, p0, p1)) < pow(10, -3));

    this->m0.push_back(m(0));
    this->m1.push_back(m(1));
}

// template <typename T>
// Eigen::VectorXd Levelset::dfdm(int pos, T* geom)
// {
//     std::vector<double> t0 = t0[pos];
//     std::vector<double> t1 = t1[pos];
//     std::vector<double> p0 = p0[pos];
//     std::vector<double> p1 = p1[pos];

//     Eigen::VectorXd c = c()
// }

template <typename T>
double Levelset::bisection(int dir, T *geom, std::vector<double> node1, std::vector<double> node2)
{
    double x;
    if (dir == 0)
    {
        x = (node1[0] + node2[0]) / 2.0;
        while (abs(geom->phi(x, node1[1])) > pow(10, -6))
        {
            if(geom->phi(x, node1[1]) > 0)
                x = (node1[0] + x) / 2.0;
            else
                x = (x + node2[0]) / 2.0;

            std::cout << x << std::endl;
            std::cout << geom->phi(x, node1[1]) << std::endl;
        }
    }
    else
    {
        x = (node1[1] + node2[1]) / 2.0;
        while (abs(geom->phi(node1[0], x)) > pow(10, -6))
        {
            if(geom->phi(node1[0], x) > 0)
                x = (node1[1] + x) / 2.0;
            else
                x = (x + node2[1]) / 2.0;

            std::cout << x << std::endl;
            std::cout << geom->phi(x, node1[1]) << std::endl;
        }
    }

    return x;
}

template <typename T>
void Levelset::find_intersection(Eigen::MatrixXd node, T *geom)
{
    Eigen::VectorXd phii(4);

    for (int i = 0; i < 4; i++)
        phii(i) = geom->phi(node(i, 0), node(i, 1));
    int iter = 0;
    if (phii(0) * phii(1) < 0)
    {
        // evaluate the intersection point
        double temp[2];
        temp[1] = node(0, 1);
        temp[0] = bisection(0, geom, {node(0, 0), node(0, 1)}, {node(1, 0), node(1, 1)});
        this->p0.push_back({temp[0], temp[1]});
        iter = 1;
        // evaluate the intersection tangent
        double dphi[2];
        dphi[0] = geom->dphidx(temp[0], temp[1]);
        dphi[1] = geom->dphidy(temp[0], temp[1]);
        this->t0.push_back({-dphi[1], dphi[0]});
    }
    if (phii(2) * phii(3) < 0)
    {
        double temp[2];
        temp[1] = node(2, 1);
        temp[0] = bisection(0, geom, {node(2, 0), node(2, 1)}, {node(3, 0), node(3, 1)});
        // evaluate the intersection tangent
        double dphi[2];
        dphi[0] = geom->dphidx(temp[0], temp[1]);
        dphi[1] = geom->dphidy(temp[0], temp[1]);
        if (iter == 0)
        {
            this->p0.push_back({temp[0], temp[1]});
            this->t0.push_back({-dphi[1], dphi[0]});
        }
        else
        {
            this->p1.push_back({temp[0], temp[1]});
            this->t1.push_back({-dphi[1], dphi[0]});
        }
    }
    if (phii(1) * phii(2) < 0)
    {
        double temp[2];
        temp[0] = node(2, 0);
        temp[1] = bisection(1, geom, {node(1, 0), node(1, 1)}, {node(2, 0), node(2, 1)});
        // evaluate the intersection tangent
        double dphi[2];
        dphi[0] = geom->dphidx(temp[0], temp[1]);
        dphi[1] = geom->dphidy(temp[0], temp[1]);
        if (iter == 0)
        {
            this->p0.push_back({temp[0], temp[1]});
            this->t0.push_back({-dphi[1], dphi[0]});
        }
        else
        {
            this->p1.push_back({temp[0], temp[1]});
            this->t1.push_back({-dphi[1], dphi[0]});
        }
    }
    if (phii(0) * phii(3) < 0)
    {
        double temp[2];
        temp[0] = node(0, 0);
        temp[1] = bisection(1, geom, {node(0, 0), node(0, 1)}, {node(3 , 0), node(3, 1)});
        // evaluate the intersection tangent
        double dphi[2];
        dphi[0] = geom->dphidx(temp[0], temp[1]);
        dphi[1] = geom->dphidy(temp[0], temp[1]);
        if (iter == 0)
        {
            this->p0.push_back({temp[0], temp[1]});
            this->t0.push_back({-dphi[1], dphi[0]});
        }
        else
        {
            this->p1.push_back({temp[0], temp[1]});
            this->t1.push_back({-dphi[1], dphi[0]});
        }
    }
}

template <typename T>
int Levelset::cut_element(Eigen::MatrixXd node, T *geom)
{
    for (int i = 0; i < node.size(); i++)
    {
        double xi = node(i, 0);
        double yi = node(i, 1);
        double phii = geom->phi(xi, yi);
        for (int j = i + 1; j < node.size(); j++)
        {
            double xj = node(j, 0);
            double yj = node(j, 1);
            double phij = geom->phi(xj, yj);
            if (phii * phij < 0)
                return 1;
        }
    }

    return 0;
}
