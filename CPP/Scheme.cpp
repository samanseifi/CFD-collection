#include <eigen3/Eigen/Dense>

#include "Scheme.h"


Scheme::Scheme() {}

Scheme::Scheme(double a){
        CFL = a;
}

Scheme::~Scheme() {}

Eigen::VectorXd Scheme::Up_Wind(Eigen::VectorXd u_old) {

        int m = u_old.size();
        Eigen::VectorXd u_n(m);

        for (int i = 1; i < m; i++) {
                u_n(i) = u_old(i) - CFL * (u_old(i) - u_old(i - 1));
        }

        // Periodic boundary condition
        u_n(0) = u_n(m - 1);

        return u_n;


}
