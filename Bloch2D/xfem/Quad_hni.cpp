#include "../include/Quad_hni.h"


int Quad_hni::get_nen() 
{
    return nen;

}

int Quad_hni::get_nnode()
{
    return nnode;
}
    
int Quad_hni::get_nelem()
{
    return nelem;
}

int Quad_hni::get_ndof()
{
    return ndof;
}

Eigen::VectorXd Quad_hni::get_coordinates(int nodenum)
{
    Eigen::VectorXd coord(2);
    coord(0) = node(nodenum, 0);
    coord(1) = node(nodenum, 1); 

    return coord;
}

int Quad_hni::get_order()
{
    return order;
}

Eigen::VectorXd Quad_hni::get_connectivity(int elenum)
{
    Eigen::VectorXd ele(get_nen());
    for(int i = 0; i < get_nen(); i++)
        ele(i) = elem(elenum, i);

    return ele;
}

Quad_hni::Quad_hni(const int nelex, const int neley, const double lx, const double ly, int order)
{
    if(order == 1)
    {
        Quad mesh = Quad(nelex, neley, lx, ly);
        copymesh(&mesh);
    }
    if(order == 2)
    {
        Quad2 mesh = Quad2(nelex, neley, lx, ly);
        copymesh(&mesh);
    }
    if(order == 3)
    {
        Quad3 mesh = Quad3(nelex, neley, lx, ly);
        copymesh(&mesh);
    }
    if(order == 4)
    {
        Quad4 mesh = Quad4(nelex, neley, lx, ly);
        copymesh(&mesh);
    }
    if(order == 5)
    {
        Quad5 mesh = Quad5(nelex, neley, lx, ly);
        copymesh(&mesh);
    }
}

void Quad_hni::add_cut_elements(int elem)
{
    cut_elem.push_back(elem);
}

int  Quad_hni::get_cut_elements(int pos)
{
    return cut_elem[pos];
}


Eigen::VectorXd Quad_hni::get_N(int i, double exi, double eta)
{
    Eigen::VectorXd N = Eigen::VectorXd::Zero(2 * order + 1);
    if(i == 1)
    {
        N(0) = 0.25;
        N(1) = -0.25 * (exi + eta);
        N(2) = 0.25 * exi * eta;
    }
    if(i == 2)
    {
        N(0) = 0.25;
        N(1) = 0.25 * (exi - eta);
        N(2) = -0.25 * exi * eta;
    }
    if(i == 3)
    {
        N(0) = 0.25;
        N(1) = 0.25 * (exi + eta);
        N(2) = 0.25 * exi * eta;
    }
    if(i == 4)
    {
        N(0) = 0.25;
        N(1) = 0.25 * (-exi + eta);
        N(2) = -0.25 * exi * eta;
    }

    return N;
}

Eigen::VectorXd Quad_hni::get_dNdx(int i, int j, double exi, double eta)
{
    Eigen::VectorXd dNdx = Eigen::VectorXd::Zero(2 * this->order);
    if(i == 1 && j== 1)
    {
        dNdx(0) = -0.25 ;
        dNdx(1) = 0.25 * eta;
    }
    else if(i == 1 && j == 2)
    {
        dNdx(0) = -0.25;
        dNdx(1) = 0.25 * exi;
    }
    else if(i == 2 && j== 1)
    {
        dNdx(0) = 0.25; 
        dNdx(1) = -0.25 * eta;
    }
    else  if(i == 2 && j == 2)
    {
        dNdx(0) = -0.25;
        dNdx(1) = -0.25 * exi;
    }
    else if(i == 3 && j== 1)
    {
        dNdx(0) = 0.25;
        dNdx(1) = 0.25 * eta;
    }
    else if(i == 3 && j == 2)
    {
        dNdx(0) = 0.25;
        dNdx(1) = 0.25 * exi;
    }
    else if(i == 4 && j == 1)
    {
        dNdx(0) = -0.25;
        dNdx(1) = -0.25 * eta; 
    }
    else if(i == 4 && j == 2)
    {
        dNdx(0) = 0.25;
        dNdx(1) = -0.25 * exi;
    }
    else
    {
        assert("Wrong choice");
    }

    return dNdx;
}

//local stiffness matrix 
//integrate over affine curve
Eigen::MatrixXd Quad_hni::get_localstiffness(Eigen::MatrixXd D, Eigen::MatrixXd x)
{
    //std::cout << x << std::endl;
    Eigen::MatrixXd ke = Eigen::MatrixXd::Zero(8, 8);

    for(int i =  0; i <= 2 * order - 1; i++)
    {
        for(int j =  0; j <= 2 * order - 1; j++)
        {
            //find no. of gauss points
            int temp = pow((i + j), 2);
            int gporder = (temp - 1) / 2; 
            if(gporder <= 0)
                gporder = 1;
            Rosetta::LegendrePolynomial gp(gporder);

            Eigen::MatrixXd B = Eigen::MatrixXd::Zero(3, 8);

            Eigen::MatrixXd ktemp = Eigen::MatrixXd::Zero(8, 8);
            //loop through each affine curve
            for(int k = 0; k < 4; k++)
            {
                double b =  0;
                double jacob;
                //loop through each gauss point
                for(int l = 0; l < gporder; l++)
                {
                    double t = (1 + gp.root(l)) / 2;
                    Eigen::VectorXd hni1(2 * order - 1), hni2(2 * order - 1), temp(2);
                    for(int m = 0; m < get_nen(); m++)
                    {
                        if(k == 0)
                        {
                            b = -node(0, 1);
                            jacob =  (x(1, 0) - x(0, 0)) / 2.0;
                            hni1 = get_dNdx(m + 1, 1, t, x(l, 1)); 
                            hni2 = get_dNdx(m + 1, 2, t, x(l, 1));
                        }
                        if(k == 1)
                        {
                            b = node(1, 0);
                            jacob = (x(2, 1) - x(1, 1)) / 2.0;
                            hni1 = get_dNdx(m + 1, 1, x(l, 0), t); 
                            hni2 = get_dNdx(m + 1, 2, x(l, 0), t);
                        }
                        if(k == 2)
                        {
                            b += node(2, 1);
                            jacob =  (x(2, 0) - x(3, 0)) / 2.0;
                            hni1 = get_dNdx(m + 1, 1, t, x(l, 1)); 
                            hni2 = get_dNdx(m + 1, 2, t, x(l, 1));
                        }
                        if(k = 3)
                        {
                            b = -node(3, 0);  
                            jacob = (node(3, 1) - node(0, 1)) / 2.0;
                            hni1 = get_dNdx(m + 1, 1, x(l, 0), t); 
                            hni2 = get_dNdx(m + 1, 2, x(l, 0), t);
                        }
                        //Ke_(i + j) = B^i D B^j 
                        //See equation 52b in the paper
                        //Modeling curved interfaces without element-partitioning in the extended finite element method
                        temp(0) = hni1(i);
                        temp(1) = hni2(j);

                        //std::cout << jmat.transpose().inverse() << std::endl;
                        Eigen::VectorXd temp2 = temp / jacob; 

                        B(0, 2 * m) += temp2(0);    
                        B(1, 2 * m + 1) += temp2(1);    
                        B(2, 2 * m) += temp2(1); 
                        B(2, 2 * m + 1) += temp2(0);
                    }
                    
                    ktemp += B.transpose() * D * B * gp.weight(l) * jacob;
                }
                ktemp *= b;
            }

            ke += ktemp / (i + j + 2);
        }
    }

    return ke;
}

//local stiffness matrix
//integrate over parametric curve
Eigen::MatrixXd Quad_hni::get_localstiffness(Eigen::MatrixXd D, Eigen::MatrixXd node, Levelset *ls, 
        std::vector<double> t0, std::vector<double> t1, std::vector<double> p0, std::vector<double>p1, 
        double m0, double m1)
{
    //std::cout << x << std::endl;
    Eigen::MatrixXd ke = Eigen::MatrixXd::Zero(get_nen() * 2, get_nen() * 2);
    
    for(int i =  0; i <= 2 * order - 1; i++)
    {
        for(int j =  0; j <= 2 * order - 1; j++)
        {
            //find no. of gauss points
            int temp = pow((i + j + 5), 2);
            int gporder = (temp - 1) / 2; 
            if(gporder <= 0)
                gporder = 1;
            Rosetta::LegendrePolynomial gp(gporder);

            Eigen::MatrixXd B = Eigen::MatrixXd::Zero(3, get_nen() * 2);

            Eigen::MatrixXd ktemp = Eigen::MatrixXd::Zero(get_nen() * 2, get_nen() * 2);
            //loop through each gauss point
            for(int l = 0; l < gporder; l++)
            {
                double t = (1 + gp.root(l)) / 2;

                Eigen::VectorXd c = ls->c(t, t0, t1, p0, p1, m0, m1);
                Eigen::VectorXd dcdt = ls->dcdt(t, t0, t1, p0, p1, m0, m1);

                for(int m = 0; m < get_nen(); m++)
                {
                    Eigen::VectorXd hni1(2 * order - 1), hni2(2 * order - 1), temp(2);

                    hni1 = get_dNdx(m + 1, 1, c(0), c(1)); 
                    hni2 = get_dNdx(m + 1, 2, c(0), c(1));

                    temp(0) = hni1(i);
                    temp(1) = hni2(j);

                    //std::cout << jmat.transpose().inverse() << std::endl;
                    Eigen::VectorXd temp2 = temp; 

                    B(0, 2 * m) += temp2(0);    
                    B(1, 2 * m + 1) += temp2(1);    
                    B(2, 2 * m) += temp2(1); 
                    B(2, 2 * m + 1) += temp2(0);
                }

                ktemp += (c(0) * dcdt(1) - c(1) * dcdt(0)) * B.transpose() * D * B * 0.5 * gp.weight(l);
            }
            //std::cout << B << std::endl;
            //std::cout << B.transpose() << std::endl;
            //std::cout << D * B << std::endl;
            //std::cout << B.transpose() * D * B << std::endl;
            //std::cout << B.transpose() * D << std::endl;
            ke += ktemp / (i + j + 2);
            //std::cout << k << std::endl;
        }
    }

    return ke;
}

//local mass matrix
//integrate over parametric curve
Eigen::MatrixXd Quad_hni::get_localmass(double rho, Eigen::MatrixXd node, Levelset *ls, 
        std::vector<double> t0, std::vector<double> t1, std::vector<double> p0, std::vector<double>p1, 
        double m0, double m1)
{
    //std::cout << x << std::endl;
    Eigen::MatrixXd me = Eigen::MatrixXd::Zero(get_nen() * 2, get_nen() * 2);

    for(int i =  0; i <= 2 * order; i++)
    {
        for(int j =  0; j <= 2 * order; j++)
        {
            //find no. of gauss points
            int temp = pow((i + j + 5), 2);
            int gporder = (temp - 1) / 2; 
            if(gporder <= 0)
                gporder = 1;
            Rosetta::LegendrePolynomial gp(gporder);

            Eigen::MatrixXd B = Eigen::MatrixXd::Zero(3, get_nen() * 2);

            Eigen::MatrixXd mtemp = Eigen::MatrixXd::Zero(get_nen() * 2, get_nen() * 2);
            //loop through each gauss point
            for(int l = 0; l < gporder; l++)
            {
                double t = (1 + gp.root(l)) / 2;

                Eigen::VectorXd c = ls->c(t, t0, t1, p0, p1, m0, m1);
                Eigen::VectorXd dcdt = ls->dcdt(t, t0, t1, p0, p1, m0, m1);

                Eigen::VectorXd hni1(2 * order), hni2(2 * order);    
                for(int m  = 0; m < get_nen(); m++)
                {
                    hni1 = get_N(m + 1, c(0), c(1)); 
                    
                    for(int n = 0; n < get_nen(); n++)
                    {
                        hni2 = get_N(n + 1, c(0), c(1));
                        mtemp(2 * m, 2 * n) = hni1(i) * hni2(j);
                        mtemp(2 * m + 1, 2 * n + 1) = mtemp(2 * m, 2 * n); 
                    }
                }

                mtemp = (c(0) * dcdt(1) - c(1) * dcdt(0)) * mtemp * 0.5 * gp.weight(l);
            }
            //std::cout << B << std::endl;
            //std::cout << B.transpose() << std::endl;
            //std::cout << D * B << std::endl;
            //std::cout << B.transpose() * D * B << std::endl;
            //std::cout << B.transpose() * D << std::endl;
            me += mtemp / (i + j + 2);
            //std::cout << k << std::endl;
        }
    }

    return rho * me;
}


//local mass matrix
//integrate over affine curve
Eigen::MatrixXd Quad_hni::get_localmass(double rho, Eigen::MatrixXd x)
{
    //std::cout << x << std::endl;
    Eigen::MatrixXd me = Eigen::MatrixXd::Zero(get_nen() * 2, get_nen() * 2);

    for(int i =  0; i <= 2 * order; i++)
    {
        for(int j =  0; j <= 2 * order; j++)
        {
            //find no. of gauss points
            int temp = pow((i + j), 2);
            int gporder = (temp - 1) / 2; 
            if(gporder <= 0)
                gporder = 1;
            Rosetta::LegendrePolynomial gp(gporder);

            Eigen::MatrixXd mtemp = Eigen::MatrixXd::Zero(get_nen() *  2, get_nen() * 2);
            //loop through each affine curve
            for(int k = 0; k < 4; k++)
            {
                double b =  0;
                double jacob;
                //loop through each gauss point
                for(int l = 0; l < gporder; l++)
                {
                    double t = (1 + gp.root(l)) / 2;
                    Eigen::VectorXd hni1(2 * order), hni2(2 * order), temp(2);
                    for(int m = 0; m < get_nen(); m++)
                    {
                        for(int n = 0; n < get_nen(); n++)
                        {
                            if(k == 0)
                            {
                                b = -node(0, 1);
                                jacob =  (x(1, 0) - x(0, 0)) / 2.0;
                                hni1 = get_N(m + 1, t, x(l, 1)); 
                                hni2 = get_N(n + 1, t, x(l, 1));
                            }
                            if(k == 1)
                            {
                                b = node(1, 0);
                                jacob = (x(2, 1) - x(1, 1)) / 2.0;
                                hni1 = get_N(m + 1, x(l, 0), t); 
                                hni2 = get_N(n + 1, x(l, 0), t);
                            }
                            if(k == 2)
                            {
                                b += node(2, 1);
                                jacob =  (x(2, 0) - x(3, 0)) / 2.0;
                                hni1 = get_N(m + 1, t, x(l, 1)); 
                                hni2 = get_N(n + 1, t, x(l, 1));
                            }
                            if(k = 3)
                            {
                                b = -node(3, 0);  
                                jacob = (node(3, 1) - node(0, 1)) / 2.0;
                                hni1 = get_N(m + 1, x(l, 0), t); 
                                hni2 = get_N(n + 1, x(l, 0), t);
                            }
                            //Ke_(i + j) = B^i D B^j 
                            //See equation 52b in the paper
                            //Modeling curved interfaces without element-partitioning in the extended finite element method
                            temp(0) = hni1(i);
                            temp(1) = hni2(j);

                            mtemp(2 * m, 2 * n) = hni1(i) * hni2(j);
                            mtemp(2 * m + 1, 2 * n + 1) = mtemp(2 * m, 2 * n); 
                        }
                    }
                    
                    mtemp = mtemp * 0.5 * gp.weight(l);
                }

                mtemp *= b;
            }

            me += mtemp / (i + j + 2);
        }
    }

    return me;
}
