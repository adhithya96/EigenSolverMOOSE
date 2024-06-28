#pragma once

#include <iostream>
#include <iomanip>
#include <cmath>
#include "eigen-master/Eigen/Eigen"

namespace Rosetta {

    /*! Implementation of Gauss-Legendre quadrature
    *  http://en.wikipedia.org/wiki/Gaussian_quadrature
    *  http://rosettacode.org/wiki/Numerical_integration/Gauss-Legendre_Quadrature
    * 
    */
    class LegendrePolynomial {
    public:
        LegendrePolynomial (int eDEGREE) {
            _r.reserve(eDEGREE + 1);
            _w.reserve(eDEGREE + 1);
            // Solve roots and weights
            for (int i = 0; i <= eDEGREE; ++i) {
                double dr = 1;

                // Find zero
                Evaluation eval(cos(M_PI * (i - 0.25) / (eDEGREE + 0.5)), eDEGREE);
                do {
                    dr = eval.v() / eval.d();
                    eval.evaluate(eval.x() - dr, eDEGREE);
                } while (fabs (dr) > 2e-16);

                this->_r[i] = eval.x();
                this->_w[i] = 2 / ((1 - eval.x() * eval.x()) * eval.d() * eval.d());
            }
        }

        double root(int i) const { return this->_r[i]; }
        double weight(int i) const { return this->_w[i]; }
        std::vector<double> roots() const { return this->_r; }
        std::vector<double> weights() const { return this->_w; }
    private:
        std::vector<double> _r;
        std::vector<double> _w;

        /*! Evaluate the value *and* derivative of the
        *   Legendre polynomial
        */
        class Evaluation {
        public:
            explicit Evaluation (double x, int eDEGREE) : _x(x), _v(1), _d(0) {
                this->evaluate(x, eDEGREE);
            }

            void evaluate(double x, int eDEGREE) {
                this->_x = x;

                double vsub1 = x;
                double vsub2 = 1;
                double f     = 1 / (x * x - 1);
            
                for (int i = 2; i <= eDEGREE; ++i) {
                    this->_v = ((2 * i - 1) * x * vsub1 - (i - 1) * vsub2) / i;
                    this->_d = i * f * (x * this->_v - vsub1);

                    vsub2 = vsub1;
                    vsub1 = this->_v;
                }
            }

            double v() const { return this->_v; }
            double d() const { return this->_d; }
            double x() const { return this->_x; }

        private:
            double _x;
            double _v;
            double _d;
        };
    };
};
=======
class gauss_points
{
private:
    Eigen::VectorXd w, x;
public:
    Eigen::VectorXd get_weights();
    Eigen::VectorXd get_points();
    gauss_points(int order);
};
