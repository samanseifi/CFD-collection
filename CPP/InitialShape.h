#ifndef INITIALSHAPE_H_
#define INITIALSHAPE_H_

#include <eigen3/Eigen/Dense>

#include "Grid.h"

class InitialShape {

public:
        InitialShape();
        InitialShape(Grid, int);
        virtual ~InitialShape();

        Eigen::VectorXd WaveFunction(Grid);
        Eigen::VectorXd StepFunction(Grid);
        void Print();

        Eigen::VectorXd u_init;

private:

};

#endif
