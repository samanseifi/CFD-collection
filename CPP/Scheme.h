#ifndef SCHEME_H_
#define SCHEME_H_

#include <eigen3/Eigen/Dense>

class Scheme {

public:
        Scheme();
        Scheme(double);
        virtual ~Scheme();

        Eigen::VectorXd Up_Wind(Eigen::VectorXd);

        Eigen::VectorXd u_new;

private:
        double CFL;

};

#endif
