//Calculate area of the circle inside a square using a level set method
#include "../include/Quad_hni.h"
#include "../include/Circle.h"
#include "../include/Levelset.h"
#include "../include/matplotlibcpp.h"
#include "../include/Quad2.h"
#include "../include/Quad3.h"
#include "../include/Quad4.h"
#include "../include/Quad5.h"
#include<vector>

int main()
{
    int nelex =  1;
    int neley = 1;
    Quad_hni mesh =  Quad_hni(nelex, neley, 1, 1);
    std::cout << "Nodal data" << std::endl;
    std::cout << mesh.get_nnode() << std::endl;
    for(int i = 0; i < mesh.get_nnode(); i++)
        std::cout << mesh.get_coordinates(i)(0) << " " << mesh.get_coordinates(i)(1) << std::endl;
    std::cout << "Number of elements per node" << std::endl;
    std::cout << mesh.get_nen() << std::endl;
    std::cout << "Number of dofs per node" << std::endl;
    std::cout << mesh.get_ndof() << std::endl;
    std::cout << "Connectivity data" << std::endl;
    for(int i = 0; i < mesh.get_nelem(); i++)
    {
        for(int j = 0; j < mesh.get_nen(); j++)
            std::cout << mesh.get_connectivity(i)(j)  << "  ";
        std::cout << std::endl;
    }

    //Get material properties 
    Material _mat = Material(7.31 * pow(10,10), 0.325, 2770, 0, 0, 0.04);

    Eigen::SparseMatrix<double, Eigen::ColMajor> K(mesh.get_nnode() * mesh.get_ndof(), mesh.get_nnode() * mesh.get_ndof()), 
    M(mesh.get_nnode() * mesh.get_ndof(), mesh.get_nnode() * mesh.get_ndof());
    
    Levelset ls = Levelset();

    Circle circle = Circle(0, 0, 0.5);

    int cut_element_size = 0;

    for(int i = 0; i < mesh.get_nelem(); i++)
    {
        Eigen::VectorXd ele = mesh.get_connectivity(i);
        std::cout << ele << std::endl;
        Eigen::MatrixXd x(ele.size(), 2);
        Eigen::VectorXd temp(2);
        for(int j = 0; j < ele.size(); j++)
        {
            temp = mesh.get_coordinates(ele(j));
            for(int k = 0; k < 2; k++)
                x(j, k) = temp(k);
        }

        //Get material data
        Eigen::MatrixXd mat = _mat.get_Dmatrix();
        std::cout << "Material matrix" << std::endl;
        std::cout << mat << std::endl;

        //find out the cut elements and parameterize the curves using splines
        int ans = ls.cut_element(x, &circle);

        if(ans == 1)
        {
            mesh.add_cut_elements(i);
            //calculate the cubic spline parameterization of cut elements.
            //p0, p1, t0, t1
            ls.find_intersection(x, &circle);
            //m0, m1
            std::vector<double> p0 = ls.get_p0(cut_element_size);
            std::cout << "p0" << std::endl;
            for(int j = 0; j < 2; j++)
                std::cout << p0[j] << " ";
            std::cout << std::endl;
            std::vector<double> p1 = ls.get_p1(cut_element_size);
            std::cout << "p1" << std::endl;
            for(int j = 0; j < 2; j++)
                std::cout << p1[j] << " ";
            std::cout << std::endl;
            std::vector<double> t0 = ls.get_t0(cut_element_size);
            std::cout << "t0" << std::endl;
            for(int j = 0; j < 2; j++)
                std::cout << t0[j] << " ";
            std::cout << std::endl;
            std::vector<double> t1 = ls.get_t1(cut_element_size);
            std::cout << "t1" << std::endl;
            for(int j = 0; j < 2; j++)
                std::cout << t1[j] << " ";
            std::cout << std::endl;
            ls.find_tangent_mag(p0, p1, t0, t1, &circle);\
            double m0 = ls.get_m0(cut_element_size);
            double m1 = ls.get_m1(cut_element_size);
            std::cout << "value of m0" << std::endl;
            std::cout << m0 << std::endl;
            std::cout << "value of m1" << std::endl;
            std::cout << m1 << std::endl;
            //print the curve points
            std::cout << "Curve points" << std::endl;
            for(int j = 0; j < 10; j++)
            {
                std::cout << "Curve point " << j + 1 << std::endl;
                std::cout << ls.c( 1.0 * j / 10.0, t0, t1, p0, p1, m0, m1) << std::endl;
            }

            //Evaluate area of mesh using HNI
            double area_affine = ls.evaluate_area_affine(x);
            std::cout << area_affine << std::endl;
            //evaluate area of curve using HNI
            double area_param = ls.evaluate_area_parametric(t0, t1, p0, p1, m0, m1);
            std::cout << area_param << std::endl;
            //Integrate using HNI
            //stiffness matrix
            Eigen::MatrixXd kaffine(8, 8), kparam(8, 8);
            kaffine = mesh.get_localstiffness(mat, x);
            kparam = mesh.get_localstiffness(mat, x, &ls, t0, t1, p0, p1, m0, m1);
            //mass matrix
            Eigen::MatrixXd maffine(8, 8), mparam(8, 8);
            maffine = mesh.get_localmass(_mat.get_rho(), x);
            mparam = mesh.get_localmass(_mat.get_rho(), x, &ls, t0, t1, p0, p1, m0, m1);
            //Assembly
            
        }
        else
        {

        }
        
        //Assembly

    }
    
    std::cout << "Assembly complete" << std::endl;
    double E = _mat.get_E();
    double nu = _mat.get_nu();
    double rho = _mat.get_rho();
    std::vector<double> wavnum, omega;
    std::fstream file;
    file.open("dispersion_data", std::ios::out | std::ios::app);
    if(file.is_open())
    {
        for(int iter = 0; iter <= 20; iter++)
        {
            //Get constraint matix and obtain reduced order mass and stiffness matrix
            //_mat.set_wavenum(0., M_PI / (4 * _mat.get_latticevec()));
            _mat.set_wavenum(0, iter * 0.05 * M_PI / _mat.get_latticevec());
            double kx = _mat.get_wavenum(1);
            double ky = _mat.get_wavenum(2);
            double a = _mat.get_latticevec();
            Eigen::MatrixXd Q = mesh.get_constraintmatrix(kx, ky, a, mesh.get_order() * nelex + 1, mesh.get_order() * neley + 1);
            
            //std::cout << "Constraint matrix" << std::endl;
            //std::cout << Q << std::endl;

            Eigen::MatrixXd Kd = Eigen::MatrixXd(K);
            Eigen::MatrixXd Md  = Eigen::MatrixXd(M);

            //std::cout << "Global stiffness matrix" << std::endl;
            //std::cout << Kd << std::endl;
            //std::cout << "Global mass matrix" << std::endl;
            //std::cout << Md << std::endl;

            Eigen::MatrixXd Kr = Q.transpose() * Kd * Q;
            Eigen::MatrixXd Mr = Q.transpose() * Md * Q;

            //std::cout << "Reduced stiffness matrix" << std::endl;
            //std::cout << Kr << std::endl;
            //std::cout << "Reduced mass matrix" << std::endl;
            //std::cout << Mr << std::endl;

            //Solve for eigenvalues
            Eigen::EigenSolver<Eigen::MatrixXd> ges(Kr * Mr.inverse());
        
            // Get the eigenvalues
            Eigen::VectorXcd eigenvalues = ges.eigenvalues();

            std::cout << eigenvalues << std::endl;

            for(int j = 0; j < eigenvalues.size(); j++)
            {
                wavnum.push_back(ky * a /  M_PI);
                double temp =  sqrt(eigenvalues(j, 0).real()) * a / sqrt(E / (3 * (1 - 2 * nu) * rho));
                omega.push_back(temp);
                file << ky * a /  M_PI << "  " << sqrt(eigenvalues(j, 0).real()) * a / sqrt(E / (3 * (1 - 2 * nu) * rho)) << std::endl;
            }


            // Output the eigenvalues
            //std::cout << "Eigenvalues: " << std::endl << eigenvalues << std::endl;
        }
    }
    else
        std::cout << "File is not open" << std::endl;

    return 0;
}
