
#include "../include/Quad2.h"
//For higher order elements Gauss Lobatto points are used 
//This includes all the nodal points in the isoparametric element.


Eigen::VectorXd Quad2::get_coordinates(int nodenum)
{
    Eigen::VectorXd coord(2);
    coord(0) = node(nodenum, 0);
    coord(1) = node(nodenum, 1); 

    return coord;
}

Eigen::VectorXd Quad2::get_connectivity(int elenum)
{
    Eigen::VectorXd ele(get_nen());
    for(int i = 0; i < get_nen(); i++)
        ele(i) = elem(elenum, i);

    return ele;
}

Quad2::Quad2(const int nelex, const int neley, const double lx, const double ly)
{
    const double dx = lx / (2 * nelex);
    const double dy = ly / (2 * neley);
    int nnodex = 2 * nelex + 1;
    int nnodey = 2 * neley + 1;
    nnode = (2 * nelex + 1) * (2 * neley + 1);
    nelem = nelex * neley;
    nen = 9;
    ndof = 4;
    ngp = 9;

    node = Eigen::MatrixXd::Zero(nnode, 2);
    int count = 0;
    for(int j = 0;  j < nnodey; j++)
    {
        for(int i = 0; i < nnodex; i++)
        {
            node(j * nnodex + i, 0) =  dx * i;
            node(j * nnodex + i, 1) =  dy * j;
        }
    }

    elem = Eigen::MatrixXd::Zero(nelem, nen);
    for(int j = 0; j < neley; j++)
    {
        for(int i = 0; i < nelex; i++)
        {
            elem(i + j * nelex, 0) = 2 * j * nnodex + 2 * i;
            elem(i + j * nelex, 1) = 2 * j * nnodex + 2 * i + 2;
            elem(i + j * nelex, 2) = 2 * (j + 1) * nnodex + 2 * i + 2;   
            elem(i + j * nelex, 3) = 2 * (j + 1) * nnodex + 2 * i;
            elem(i + j * nelex, 4) = 2 * j * nnodex + 2 * i + 1;
            elem(i + j * nelex, 5) = (2 * j + 1) * nnodex + 2 * i + 2;
            elem(i + j * nelex, 6) = 2 * (j + 1) * nnodex + 2 * i + 1;   
            elem(i + j * nelex, 7) = (2 * j + 1) * nnodex + 2 * i;
            elem(i + j * nelex, 8) = (2 * j + 1) * nnodex + 2 * i + 1;
        }
    }
}

Eigen::VectorXd Quad2::RefElementCoord(int nodenum)
{
    Eigen::VectorXd x(2);
    switch(nodenum)
    {
        case 1: 
            x(0) = -1;
            x(1) = -1;
            break;
        case 2:
            x(0) = 1;
            x(1) = -1;
            break;
        case 3:
            x(0) = 1;
            x(1) = 1;
            break;
        case 4:
            x(0) = -1;
            x(1) = 1;
            break;
        case 5:
            x(0) = 0;
            x(1) = -1;
            break;
        case 6:
            x(0) = 1;
            x(1) = 0;
            break;
        case 7:
            x(0) = 0;
            x(1) = 1;
            break;
        case 8:
            x(0) = -1;
            x(1) = 0;
            break;
        case 9:
            x(0) = 0;
            x(1) = 0;
            break;
        default:
            x(0) = -99;
            x(1) = -99;
            break;
    }

    return x;
}

int Quad2::get_nen() 
{
    return nen;

}

int Quad2::get_nnode()
{
    return nnode;
}
    
int Quad2::get_nelem()
{
    return nelem;
}

int Quad2::get_ndof()
{
    return ndof;
}

int Quad2::get_ngp()
{
    return ngp;
}

Eigen::VectorXd Quad2::GaussIntegrationWeights(int nodenum)
{
    Eigen::VectorXd w(2);
    switch(nodenum)
    {
        case 1: 
            w(0) = 0.333;
            w(1) = 0.333;
            break;
        case 2:
            w(0) = 0.333;
            w(1) = 0.333;
            break;
        case 3:
            w(0) = 0.333;
            w(1) = 0.333;
            break;
        case 4:
            w(0) = 0.333;
            w(1) = 0.333;
            break;
        case 5:
            w(0) = 1.333;
            w(1) = 0.333;
            break;
        case 6:
            w(0) = 0.333;
            w(1) = 1.333;
            break;
        case 7:
            w(0) = 1.333;
            w(1) = 0.333;
            break;
        case 8:
            w(0) = 0.333;
            w(1) = 1.333;
            break;
        case 9:
            w(0) = 1.333;
            w(1) = 1.333;
            break;
        default:
            w(0) = -99;
            w(1) = -99;
            break;
    }
    return w;
}

//Evaluating value of shape function at node i
double Quad2::get_N(int i, double exi, double eta)
{
    //Get the reference  coordinates
    Eigen::VectorXd x = RefElementCoord(i);


    Eigen::VectorXd  refx(3);
    refx << -1,0,1;
    Eigen::VectorXd refy(3); 
    refy << -1, 0, 1;
    
    double val = 1;
    for(int k = 0; k < 3; k++)
    {
        if(refx(k) != x(0))
            val *= (exi - refx(k)) / (x(0) - refx(k));
        if(refy(k) != x(1))
            val *= (eta - refy(k)) / (x(1) - refy(k));
    }

    return val;
}

//dNi/dexi if j = 1 at exi and eta
//dNi/deta if j = 2 at exi and eta
double Quad2::get_dNdx(int i, int j, double exi, double eta)
{
    Eigen::VectorXd x = RefElementCoord(i);
    double exii = x(0), etaj = x(1);
//    std::cout << "Element coordinates" << std::endl;
//    std::cout << x << std::endl;
    //std::cout << x << std::endl;
    Eigen::VectorXd exik(3);
    exik << -1,0,1;
    Eigen::VectorXd etal(3);
    etal << -1, 0, 1;
    
    double val;
    double denexi = 1, deneta = 1, numexi = 1, numeta = 1;
    for(int iter = 0; iter < 3; iter++)
    {
        if(exik(iter) != exii)
        {
            denexi *= (exii - exik(iter));
            numexi *= (exi - exik(iter));
        }
        if(etal(iter) != etaj)
        {
            deneta *= (etaj - etal(iter));
            numeta *= (eta - etal(iter));
        }
    }

//    std::cout << "exi denominator" << std::endl;
//    std::cout << denexi << std::endl;
//    std::cout << "eta denominator" << std::endl;
//    std::cout << deneta << std::endl;

    if(j == 1)
    {
        numexi = 0;
        for(int iter =  0; iter < 3; iter++)
        {
            if(exik(iter) != exii)
            {
                double temp = 1;
                for(int iter2 = 0; iter2  < 3; iter2++)
                {
                    if((iter2 != iter)&&(exik(iter2) != exii))
                        temp *= (exi - exik(iter2));
                }
                numexi += temp;
            }
        }
    }
    
    if(j == 2)
    {
        numeta = 0;
        for(int iter = 0; iter < 3; iter++)
        {
            if(etal(iter) != etaj)
            {
                double temp = 1;
                for(int iter2 = 0; iter2  < 3; iter2++)
                {
                    if((iter2 != iter)&&(etal(iter2) != etaj))
                        temp *= (eta - etal(iter2));
                }
                numeta += temp;
            }
        }

    }  

//    std::cout << "exi numerator" << std::endl;
//    std::cout << numexi << std::endl;
//    std::cout << "eta numerator" << std::endl;
//    std::cout << numeta << std::endl;

    val = numexi * numeta / (denexi * deneta);
    
    return val;
}

Eigen::MatrixXd Quad2::get_jacobianmat(Eigen::MatrixXd x, double exi, double eta)
{
    Eigen::VectorXd xe(get_nen()), ye(get_nen());
    for(int i = 0; i < get_nen(); i++)
    {
        xe(i) = x(i, 0);
        ye(i) = x(i, 1);
    }
//    std::cout << "x coordinate"  << std::endl;
//    std::cout << xe << std::endl;
//    std::cout << "y coordinate" << std::endl;
//    std::cout << ye << std::endl;
    double dxdexi = 0., dxdeta = 0., dydexi = 0., dydeta = 0.; 
    Eigen::MatrixXd jmat = Eigen::MatrixXd::Zero(2, 2);
    for(int j = 0; j < get_nen(); j++)
    {
        dxdexi += get_dNdx(j + 1, 1, exi, eta) * xe(j);
        dxdeta += get_dNdx(j + 1, 2, exi, eta) * xe(j);
        dydexi += get_dNdx(j + 1, 1, exi, eta) * ye(j);
        dydeta += get_dNdx(j + 1, 2, exi, eta) * ye(j);
    }

    jmat(0, 0) = dxdexi;
    jmat(0, 1) = dxdeta;
    jmat(1, 0) = dydexi;
    jmat(1, 1) = dydeta;

//    std::cout << "jacobian matrix at quadrature point " << exi << "  " << eta << std::endl;
//    std::cout << jmat << std::endl;

    assert(jmat.determinant() > 0);

    return jmat;
}

Eigen::MatrixXd Quad2::get_localstiffness(Eigen::MatrixXd D, Eigen::MatrixXd node)
{
    //gaussian points
    Eigen::VectorXd x(2);

    //std::cout << x << std::endl;
    Eigen::MatrixXd k = Eigen::MatrixXd::Zero(get_nen() * 2, get_nen() * 2);
   
    for(int gp =  0; gp < get_ngp(); gp++)
    {
        Eigen::MatrixXd B = Eigen::MatrixXd::Zero(3, get_nen() * 2);

        //get gauss lobatto point locations
        x =  RefElementCoord(gp + 1);
        std::cout << x << std::endl;

        Eigen::MatrixXd jmat = get_jacobianmat(node, x(0), x(1));
        
        Eigen::VectorXd temp(2);
        
        //get gauss lobatto point weights
        Eigen::VectorXd w =  GaussIntegrationWeights(gp + 1);
        std::cout << w << std::endl;

        for(int j = 0; j < get_nen(); j++)
        {
            temp(0) = get_dNdx(j + 1, 1, x(0), x(1)); 
            temp(1) = get_dNdx(j + 1, 2, x(0), x(1));

            std::cout << jmat.transpose().inverse() << std::endl;
            Eigen::VectorXd temp2 = jmat.transpose().inverse() * temp;

            B(0, 2 * j) += temp2(0); 
            B(1, 2 * j + 1) += temp2(1); 
            B(2, 2 * j) += temp2(1); 
            B(2, 2 * j + 1) += temp2(0); 
        }
        std::cout << B << std::endl;
        //std::cout << B.transpose() << std::endl;
        //std::cout << D * B << std::endl;
         //std::cout << B.transpose() * D * B << std::endl;
        //std::cout << B.transpose() * D << std::endl;
        k += w(0) * w(1) * B.transpose() * D * B * jmat.determinant();
        //std::cout << k << std::endl;
    }

    return k;
}

Eigen::MatrixXd Quad2::get_localmass(double rho, Eigen::MatrixXd node)
{
    //std::cout << x << std::endl;
    Eigen::MatrixXd me = Eigen::MatrixXd::Zero(get_nen() * 2, get_nen() * 2);

    for(int gp =  0; gp < get_ngp();  gp++)
    {
        Eigen::VectorXd x(2);
        //gauss lobatto points
        x = RefElementCoord(gp + 1);
        
        //Get Gauss lobatto weights
        Eigen::VectorXd  w(2);
        w = GaussIntegrationWeights(gp + 1);

        Eigen::MatrixXd m = Eigen::MatrixXd::Zero(get_nen() * 2, get_nen() * 2);

        for(int j  = 0; j < get_nen(); j++)
        {
            for(int k = 0; k < get_nen(); k++)
            {
                m(2 * j, 2 * k) = get_N(j + 1, x(0), x(1)) * get_N(k + 1, x(0), x(1));
                m(2 * j + 1, 2 * k + 1) = m(2 * j, 2 * k); 
            }
        }
        
        Eigen::MatrixXd jmat = get_jacobianmat(node, x(0), x(1));

        me += w(0) * w(1) * jmat.determinant() * m;
     
    }

    std::cout  << me << std::endl;

    return rho * me;
}


Eigen::MatrixXd Quad2::get_constraintmatrix(double kx, double ky, double a, int nnodex, int nnodey)
{
    Eigen::MatrixXd Q = Eigen::MatrixXd::Zero(get_ndof() * nnode, get_ndof() * (nnodex - 1) * (nnodey  - 1));
    int iter = 0, ndim = 2;

    std::cout << "wave  number " << " " << "kx" << "  " << "ky" << std::endl;
    std::cout << kx << "   " << ky << std::endl;
    std::cout << "lattice length" << std::endl;
    std::cout << a << std::endl;
    int temp = 0;
    for(int i =  0; i < nnodex; i++)
    {
        for(int j = 0; j < nnodey; j++)
        {
            if(j == nnodex - 1 && i != nnodey - 1) 
            {
                int xrur = ndim * (j + nnodex * i);
                int yrur = ndim * ((nnodex - 1) * i);
                int yiur = yrur + ndim * (nnodex - 1) * (nnodey - 1);
                int xrui = xrur +  ndim * nnodex * nnodey;

                std::cout << xrur << std::endl;
                std::cout << yrur << std::endl;
                std::cout << yiur << std::endl;
                std::cout << xrui << std::endl;

                //ur
                Q(xrur, yrur) = cos(kx * a);
                Q(xrur, yiur) = -sin(kx * a);
                //ui
                Q(xrui, yrur) = sin(kx * a);            
                Q(xrui, yiur) = cos(kx * a);
                //vr
                Q(xrur + 1,yrur + 1) = cos(kx * a);
                Q(xrur + 1,yiur + 1) = -sin(kx * a);
                //vi
                Q(xrui + 1, yrur + 1) = sin(kx * a);            
                Q(xrui + 1, yiur + 1) = cos(kx * a);

                std::cout << Q << std::endl;
            }
            else if (i == nnodey - 1  && j != nnodex - 1)
            {
                int xrur = ndim * (i * nnodex + j);
                int yrur = ndim * j;
                int yiur = yrur + ndim * (nnodex - 1) * (nnodey - 1);
                int xrui = xrur +  ndim * nnodex * nnodey;

                std::cout << xrur << std::endl;
                std::cout << yrur << std::endl;
                std::cout << yiur << std::endl;
                std::cout << xrui << std::endl;

                //ur
                Q(xrur, yrur) = cos(ky * a);
                Q(xrur, yiur) = -sin(ky * a);
                //ui
                Q(xrui, yrur) = sin(ky * a);            
                Q(xrui, yiur) = cos(ky * a);
                //vr
                Q(xrur + 1,yrur + 1) = cos(ky * a);
                Q(xrur + 1,yiur + 1) = -sin(ky * a);
                //vi
                Q(xrui + 1, yrur + 1) = sin(ky * a);            
                Q(xrui + 1, yiur + 1) = cos(ky * a);
            }
            else if(i == nnodey - 1 && j == nnodex - 1)
            {
                int xrur = ndim * (i * nnodex + j);
                int yrur = 0;
                int yiur = yrur + ndim * (nnodex - 1) * (nnodey - 1);
                int xrui = xrur +  ndim * nnodex * nnodey;

                std::cout << xrur << std::endl;
                std::cout << yrur << std::endl;
                std::cout << yiur << std::endl;
                std::cout << xrui << std::endl;

                //ur
                Q(xrur, yrur) = cos((kx + ky) * a);
                Q(xrur, yiur) = -sin((kx + ky) * a);
                //ui
                Q(xrui, yrur) = sin((kx + ky) * a);            
                Q(xrui, yiur) = cos((kx + ky) * a);
                //vr
                Q(xrur + 1,yrur + 1) = cos((kx + ky) * a);
                Q(xrur + 1,yiur + 1) = -sin((kx + ky) * a);
                //vi
                Q(xrui + 1, yrur + 1) = sin((kx + ky) * a);            
                Q(xrui + 1, yiur + 1) = cos((kx + ky) * a);

                std::cout << Q << std::endl;            
            }
            else
            {
                int xrur = ndim * (j + nnodey * i);
                int yrur = ndim * temp;
                //ur
                Q(xrur, yrur) = 1;
                //ui
                int xiur = xrur + ndim * nnodex * nnodey;
                int yiur = yrur + ndim * (nnodex - 1) * (nnodey - 1);
                Q(xiur, yiur) = 1;            
                //vr
                Q(xrur + 1, yrur + 1) = 1;
                //vi
                Q(xiur + 1, yiur + 1) = 1;
                temp++;

                std::cout << xrur << std::endl;
                std::cout << yrur << std::endl;
                std::cout << yiur << std::endl;
                std::cout << xiur << std::endl;

                std::cout << Q << std::endl;

            }


        }
    }

    return Q;
}
