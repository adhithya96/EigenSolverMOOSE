#include "../include/Quad.h"
#include "../include/matplotlibcpp.h"
#include "../include/Quad2.h"
#include "../include/Quad3.h"
#include "../include/Quad4.h"
#include "../include/Quad5.h"

namespace plt = matplotlibcpp;

int main()
{
    int choice;
    
    //Get mesh data
    //Quad mesh = Quad("square_mesh_coarse.vtk");
    int nelex =  1;
    int neley = 1;
    Quad5 mesh =  Quad5(nelex, neley, 0.04, 0.04);
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

    /*for(int j = 0; j < mesh.get_nen(); j++)
    {
        std::cout << "shape function of node " << j + 1 <<  " at node " << j + 1 << std::endl;
        std::cout << mesh.get_N(j + 1, mesh.RefElementCoord(j + 1)(0), mesh.RefElementCoord(j + 1)(1)) << std::endl;   

        std::cout << "Derivative of shape function of node " << j + 1 << " at node " << j + 1 << " in direction " << 1 << std::endl;
        std::cout << mesh.get_dNdx(j + 1, 1, mesh.RefElementCoord(j + 1)(0), mesh.RefElementCoord(j + 1)(1)) << std::endl;

        std::cout << "Derivative of shape function of node " << j + 1 << " at node " << j + 1 << " in direction " << 2 << std::endl;
        std::cout << mesh.get_dNdx(j + 1, 2, mesh.RefElementCoord(j + 1)(0), mesh.RefElementCoord(j + 1)(1)) << std::endl;

    }*/

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

        std::cout << x << std::endl;
    
        //Get material data
        Eigen::MatrixXd mat = _mat.get_Dmatrix();
        std::cout << "Material matrix" << std::endl;
        std::cout << mat << std::endl;

        //Get local stiffness matrix
        Eigen::MatrixXd ke = mesh.get_localstiffness(mat, x);
        std::cout << "Local stiffness matrix" << std::endl;
        std::cout << ke << std::endl;

        //Get local mass matrix
        Eigen::MatrixXd me = mesh.get_localmass(_mat.get_rho(), x); 
        std::cout << "local mass matrix" << std::endl;
        std::cout << me << std::endl;

        //Assembly
        int nnode = mesh.get_nnode();
        for(int j = 0; j < ele.size(); j++)
            for(int k = 0; k < ele.size(); k++)
                for(int l = 0; l < 2; l++)
                    for(int m = 0; m < 2; m++)
                    {
                        //Kreal
                        K.coeffRef(2 * ele(j) + l, 2 * ele(k) + m) += ke(l + j * 2, m + k * 2);
                        //Mreal
                        M.coeffRef(2 * ele(j) + l, 2 * ele(k) + m) += me(l + j * 2, m + k  * 2);
                        //Kimag
                        K.coeffRef(2 * ele(j) + l + nnode * 2, 2 * ele(k) + m + nnode * 2) += ke(l + j * 2, m + k * 2);
                        //Mimag
                        M.coeffRef(2 * ele(j) + l + nnode * 2, 2 * ele(k) + m + nnode * 2) += me(l + j * 2, m + k * 2);
                    }           
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

    std::cout << "Finished evaluating the dispersion curve" << std::endl;
    std::vector<double> ka, wa;
    double c = sqrt(E/rho);
    double n = 1.;

    plt::figure(); 

    plt::ylim(0, 6);

    /*double m = -1;
    double a = _mat.get_latticevec(); 
    for(int iter = 0; iter <= 20; iter++)
    {
        _mat.set_wavenum(0, iter * 0.05 * M_PI / a);
        double ky = _mat.get_wavenum(2);
        double temp = c * sqrt(pow((n * M_PI * 2 / a),2) + pow((ky + (m * M_PI * 2) / a),2));
        wa.push_back(temp * a / sqrt(E/(3*(1 - 2*nu)*rho)));
        ka.push_back(ky * a / M_PI);
    }

    plt::plot(ka, wa);*/


    plt::scatter(wavnum, omega);                        

    plt::title("Dispersion diagram");

    plt::xlabel("wave number");
    plt::ylabel("angular frequency");

    plt::savefig("Dispersion_Curve.pdf");

    return 0;

//    Eigen::SparseLU<Eigen::SparseMatrix<double, Eigen::ColMajor>> solver;
//    solver.analyzePattern(M);
//    solver.factorize(M); //LU decomposition

//    assert(solver.info() == Eigen::Success);

//    Eigen::VectorXd lambda = solver.solve(K);

//    std::cout << lambda << std::endl;
}