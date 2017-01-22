#include <eigen3/Eigen/Dense>
#include <iostream>
#include <string>

#include "InitialShape.h"

InitialShape::InitialShape() {

}

InitialShape::InitialShape(Grid grid, int code): u_init(grid.num_nodes()) {
        // code = 1: Wave function
        // code = 2: Step function
        if (code == 1) {
                u_init = WaveFunction(grid);
        } else if (code == 2) {
                u_init = StepFunction(grid);
        }
}

InitialShape::~InitialShape() {


}

Eigen::VectorXd InitialShape::WaveFunction(Grid grid) {
        Eigen::VectorXd u_0(grid.num_nodes());
        for (int i = 0; i < grid.num_nodes(); i++) {
                if (grid.x_init(i) > 0.25 && grid.x_init(i) < 0.75) {
                        u_0(i) = 1.0;
                } else {
                        u_0(i) = 0.0;
                }
        }
        return u_0;
}

Eigen::VectorXd InitialShape::StepFunction(Grid grid) {
        Eigen::VectorXd u_0(grid.num_nodes());
        for (int i = 0; i < grid.num_nodes(); i++) {
                if (grid.x_init(i) < 0.75) {
                        u_0(i) = 1.0;
                } else {
                        u_0(i) = 0.0;
                }
        }
        return u_0;
}

void InitialShape::Print() {
        //std::cout << u_init << std::endl;
}
