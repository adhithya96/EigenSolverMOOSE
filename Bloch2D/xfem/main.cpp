//Calculate area of the circle inside a square using a level set method
#include "../include/Quad.h"
#include "../include/Circle.h"
#include "../include/Levelset.h"
#include<vector>

int main()
{
    int nelex =  1;
    int neley = 1;
    Quad mesh =  Quad(nelex, neley, 1, 1);
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

        Levelset ls = Levelset();
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
            ls.find_tangent_mag(p0, p1, t0, t1, &circle);
            //print the curve points
            std::cout << "Curve points" << std::endl;
            for(int j = 0; j < 10; j++)
            {
                std::cout << "Curve point " << j + 1 << std::endl;
                std::cout << ls.c( 1.0 * j / 10.0, cut_element_size) << std::endl;
            }

            //Evaluate area of mesh using HNI
            double area_affine = ls.evaluate_area_affine(x);
            std::cout << area_affine << std::endl;
            //evaluate area of curve usnig HNI
            double area_param = ls.evaluate_area_parametric(cut_element_size);
            std::cout << area_param << std::endl;

        }
    }


    return 0;
}
