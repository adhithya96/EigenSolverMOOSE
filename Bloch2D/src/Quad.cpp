#include "../include/Quad.h"

Quad::Quad(const int nelex, const int neley, const double lx, const double ly)
{
    nen = 4;
    ndof = 4;
    order = 1;
    const double dx = lx / nelex;
    const double dy = ly / neley;
    nnode = (nelex + 1) * (neley + 1);
    nelem = nelex * neley;

    node = Eigen::MatrixXd::Zero(nnode, 2);
    int count = 0;
    for(int j = 0;  j < neley + 1; j++)
    {
        for(int i = 0; i < nelex + 1; i++)
        {
            node(j * (nelex + 1) + i, 0) =  dx * i;
            node(j * (nelex + 1) + i, 1) =  dy * j;
        }
    }

    elem = Eigen::MatrixXd::Zero(nelem, nen);
    for(int j = 0; j < neley; j++)
    {
        for(int i = 0; i < nelex; i++)
        {
            elem(i + j * nelex, 0) = j * (nelex + 1) + i;
            elem(i + j * nelex, 1) = j * (nelex + 1) + i + 1;
            elem(i + j * nelex, 2) = (j + 1) * (nelex + 1)  + i + 1;   
            elem(i + j * nelex, 3) = (j + 1) * (nelex + 1) + i;
        }
    }
}

Quad::Quad(std::string file_name)
{
    std::fstream file;
    file.open(file_name);
    std::string line;
    //read nodal and elemental data from mesh file in vtk format 
    while(file.is_open() && !file.eof())
    {
        std::getline(file, line);
        std::stringstream s(line);

        std::string temp;
        s >> temp;
        //std::cout << line << std::endl;
        //std::cout << temp << std::endl;
        if(temp == "POINTS")
        {
            s >> nnode;
            //std::cout << nnode << std::endl;
            node = Eigen::MatrixXd(nnode, 2);
            for(int i = 0; i < nnode; i++)
            {
                getline(file, line);
                std::stringstream s(line);
                s >> node(i, 0) >> node(i, 1);
                //std::cout << node(i, 0) << node(i, 1) << std::endl;
            }
        }
        else if(temp == "CELLS")
        {
            s >> nelem;
            elem = Eigen::MatrixXd(nelem, 4);
            for(int i = 0; i < nelem; i++)
            {
                getline(file, line);
                std::stringstream s(line);
                int temp;
                s >> temp >> elem(i, 0) >> elem(i, 1) >> elem(i, 2) >> elem(i, 3);
                //std::cout << temp << elem(i, 0) << elem(i, 1) << elem(i, 2) << elem(i, 3);
            }
        }
        else
            continue;
    }
    if(!file.is_open())
        std::cout << "Unable to open file " << std::endl;
}

int Quad::get_order()
{
    return order;
}

Eigen::VectorXd Quad::get_coordinates(int nodenum)
{
    Eigen::VectorXd coord(2);
    coord(0) = node(nodenum, 0);
    coord(1) = node(nodenum, 1); 

    return coord;
}

int Quad::get_nen()
{
    return nen;

}

int Quad::get_nnode()
{
    return nnode;
}
    
int Quad::get_nelem()
{
    return nelem;
}

int Quad::get_ndof()
{
    return ndof;
}

Eigen::VectorXd Quad::get_connectivity(int elenum)
{
    Eigen::VectorXd ele(4);
    ele(0) = elem(elenum, 0);
    ele(1) = elem(elenum, 1);
    ele(2) = elem(elenum, 2);
    ele(3) = elem(elenum, 3);

    return ele;
}

double Quad::get_N(int i, double exi, double eta)
{
    if(i == 1)
        return 0.25 * (1 - exi) * (1 - eta);
    else if (i == 2)
        return 0.25 * (1 + exi) * (1 - eta);
    else if (i == 3)
        return 0.25 * (1 + exi) * (1 + eta);
    else if (i == 4)
        return 0.25 * (1 - exi) * (1 + eta);
    else
    {
        assert("wrong choice");
        return -99;
    }
}

double Quad::get_dNdx(int i, int j, double exi, double eta)
{
    if(i == 1 && j== 1)
        return -0.25 * (1 - eta);
    else if(i == 1 && j == 2)
        return -0.25 * (1 - exi);
    else if(i == 2 && j== 1)
        return 0.25 * (1 - eta);
    else  if(i == 2 && j == 2)
        return -0.25 * (1 + exi);
    else if(i == 3 && j== 1)
        return 0.25 * (1 + eta);
    else if(i == 3 && j == 2)
        return 0.25 * (1 + exi);
    else if(i == 4 && j == 1)
        return -0.25 * (1 + eta);
    else if(i == 4 && j == 2)
        return 0.25 * (1 - exi);
    else
    {
        assert("Wrong choice");
        return -99;
    }
}

void Quad::add_cut_elements(int elem)
{
    cut_elem.push_back(elem);
}

int  Quad::get_cut_elements(int pos)
{
    return cut_elem[pos];
}


Eigen::MatrixXd Quad::get_jacobianmat(Eigen::MatrixXd x, double exi, double eta)
{

    Eigen::VectorXd xe(4), ye(4);
    for(int i = 0; i < 4; i++)
    {
        xe(i) = x(i, 0);
        ye(i) = x(i, 1);
    }
    //std::cout << "x coordinate"  << std::endl;
    //std::cout << xe << std::endl;
    //std::cout << "y coordinate" << std::endl;
    //std::cout << ye << std::endl;
    double dxdexi = 0., dxdeta = 0., dydexi = 0., dydeta = 0.; 
    Eigen::MatrixXd jmat = Eigen::MatrixXd::Zero(2, 2);
    for(int j = 0; j < 4; j++)
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

    //std::cout << "jacobian matrix at quadrature point " << exi << "  " << eta << std::endl;
    //std::cout << jmat << std::endl;

    return jmat;
}

Eigen::MatrixXd Quad::get_localstiffness(Eigen::MatrixXd D, Eigen::MatrixXd node)
{
    //weights
    double w1 = 1., w2 =  1., w3 = 1., w4 = 1.;
    //gaussian points
    Eigen::MatrixXd x(4, 2);
    x << -0.5774, -0.5774,
        0.5744, -0.5744,
        0.5744, 0.5744,
        -0.5744, 0.5744;
    //std::cout << x << std::endl;
    Eigen::MatrixXd k = Eigen::MatrixXd::Zero(8, 8);
   
    for(int i =  0; i < 4; i++)
    {
        Eigen::MatrixXd B = Eigen::MatrixXd::Zero(3, 8);

        Eigen::MatrixXd jmat = get_jacobianmat(node, x(i, 0), x(i, 1));

        for(int j = 0; j < 4; j++)
        {
            Eigen::VectorXd temp(2);
            temp(0) = get_dNdx(j + 1, 1, x(i, 0), x(i, 1)); 
            temp(1) = get_dNdx(j + 1, 2, x(i, 0), x(i, 1));

            //std::cout << jmat.transpose().inverse() << std::endl;
            Eigen::VectorXd temp2 = jmat.transpose().inverse() * temp;

            B(0, 2 * j) += temp2(0); 
            B(1, 2 * j + 1) += temp2(1); 
            B(2, 2 * j) += temp2(1); 
            B(2, 2 * j + 1) += temp2(0); 
        }
        //std::cout << B << std::endl;
        //std::cout << B.transpose() << std::endl;
        //std::cout << D * B << std::endl;
         //std::cout << B.transpose() * D * B << std::endl;
        //std::cout << B.transpose() * D << std::endl;
        k += B.transpose() * D * B * jmat.determinant();
        //std::cout << k << std::endl;
    }

    return k;
}

Eigen::MatrixXd Quad::get_localmass(double rho, Eigen::MatrixXd node)
{
        //weights
    double w1 = 1., w2 =  1., w3 = 1., w4 = 1.;
    //gaussian points
    Eigen::MatrixXd x(4, 2);
    x << -0.5774, -0.5774,
        0.5744, -0.5744,
        0.5744, 0.5744,
        -0.5744, 0.5744;
    Eigen::MatrixXd me =  Eigen::MatrixXd::Zero(8, 8);

    for(int i  = 0; i < 4; i++)
    {
        Eigen::MatrixXd m = Eigen::MatrixXd::Zero(8, 8);

        m(0, 0) = pow(get_N(1, x(i, 0), x(i, 1)), 2);
        m(0, 2) = get_N(1, x(i, 0), x(i, 1)) * get_N(2, x(i, 0), x(i, 1));
        m(0, 4) = get_N(1, x(i, 0), x(i, 1)) * get_N(3, x(i, 0), x(i, 1));
        m(0, 6) = get_N(1, x(i, 0), x(i, 1)) * get_N(4, x(i, 0), x(i, 1));

        m(1, 1) = pow(get_N(1, x(i, 0), x(i, 1)), 2);
        m(1, 3) = get_N(1, x(i, 0), x(i, 1)) * get_N(2, x(i, 0), x(i, 1));
        m(1, 5) = get_N(1, x(i, 0), x(i, 1)) * get_N(3, x(i, 0), x(i, 1));
        m(1, 7) = get_N(1, x(i, 0), x(i, 1)) * get_N(4, x(i, 0), x(i, 1));

        m(2, 0) = get_N(1, x(i, 0), x(i, 1)) * get_N(2, x(i, 0), x(i, 1));
        m(2, 2) = pow(get_N(2, x(i, 0), x(i, 1)), 2);
        m(2, 4) = get_N(2, x(i, 0), x(i, 1)) * get_N(3, x(i, 0), x(i, 1));
        m(2, 6) = get_N(2, x(i, 0), x(i, 1)) * get_N(4, x(i, 0), x(i, 1));

        m(3, 1) = get_N(1, x(i, 0), x(i, 1)) * get_N(2, x(i, 0), x(i, 1));
        m(3, 3) = pow(get_N(2, x(i, 0), x(i, 1)), 2);
        m(3, 5) = get_N(2, x(i, 0), x(i, 1)) * get_N(3, x(i, 0), x(i, 1));
        m(3, 7) = get_N(2, x(i, 0), x(i, 1)) * get_N(4, x(i, 0), x(i, 1));

        m(4, 0) = get_N(1, x(i, 0), x(i, 1)) * get_N(3, x(i, 0), x(i, 1));
        m(4, 2) = get_N(2, x(i, 0), x(i, 1)) * get_N(3, x(i, 0), x(i, 1));
        m(4, 4) = pow(get_N(3, x(i, 0), x(i, 1)), 2);
        m(4, 6) = get_N(3, x(i, 0), x(i, 1)) * get_N(4, x(i, 0), x(i, 1));

        m(5, 1) = get_N(1, x(i, 0), x(i, 1)) * get_N(3, x(i, 0), x(i, 1));
        m(5, 3) = get_N(2, x(i, 0), x(i, 1)) * get_N(3, x(i, 0), x(i, 1));
        m(5, 5) = pow(get_N(3, x(i, 0), x(i, 1)), 2);
        m(5, 7) = get_N(3, x(i, 0), x(i, 1)) * get_N(4, x(i, 0), x(i, 1));

        m(6, 0) = get_N(1, x(i, 0), x(i, 1)) * get_N(4, x(i, 0), x(i, 1));
        m(6, 2) = get_N(2, x(i, 0), x(i, 1)) * get_N(4, x(i, 0), x(i, 1));
        m(6, 4) = get_N(3, x(i, 0), x(i, 1)) * get_N(4, x(i, 0), x(i, 1));
        m(6, 6) = pow(get_N(4, x(i, 0), x(i, 1)), 2);

        m(7, 1) = get_N(1, x(i, 0), x(i, 1)) * get_N(4, x(i, 0), x(i, 1));
        m(7, 3) = get_N(2, x(i, 0), x(i, 1)) * get_N(4, x(i, 0), x(i, 1));
        m(7, 5) = get_N(3, x(i, 0), x(i, 1)) * get_N(4, x(i, 0), x(i, 1));
        m(7, 7) = pow(get_N(4, x(i, 0), x(i, 1)), 2);

        Eigen::MatrixXd jmat = get_jacobianmat(node, x(i, 0), x(i, 1));

        me += jmat.determinant() * m;
    }
    
    //std::cout  << me << std::endl;

    return rho * me;
}

Eigen::MatrixXd Quad::get_constraintmatrix(double kx, double ky, double a, int nnodex, int nnodey)
{
    Eigen::MatrixXd Q = Eigen::MatrixXd::Zero(4 * nnode, 4 * (nnodex - 1) * (nnodey  - 1));
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

                //std::cout << xrur << std::endl;
                //std::cout << yrur << std::endl;
                //std::cout << yiur << std::endl;
                //std::cout << xrui << std::endl;

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

                //std::cout << Q << std::endl;
            }
            else if (i == nnodey - 1  && j != nnodex - 1)
            {
                int xrur = ndim * (i * nnodex + j);
                int yrur = ndim * j;
                int yiur = yrur + ndim * (nnodex - 1) * (nnodey - 1);
                int xrui = xrur +  ndim * nnodex * nnodey;

                //std::cout << xrur << std::endl;
                //std::cout << yrur << std::endl;
                //std::cout << yiur << std::endl;
                //std::cout << xrui << std::endl;

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

                //std::cout << xrur << std::endl;
                //std::cout << yrur << std::endl;
                //std::cout << yiur << std::endl;
                //std::cout << xrui << std::endl;

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

                //std::cout << Q << std::endl;
            
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

                //std::cout << xrur << std::endl;
                //std::cout << yrur << std::endl;
                //std::cout << yiur << std::endl;
                //std::cout << xiur << std::endl;

                //std::cout << Q << std::endl;

            }


        }
    }

    return Q;
}
