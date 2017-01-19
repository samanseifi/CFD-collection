#include <iostream>
#include <eigen3/Eigen/Dense>

#include "InitialShape.h"
#include "Grid.h"

int main() {

      Grid grid(20);
      InitialShape IC(grid, 1);
      std::cout << IC.u_init << std::endl;
      std::cout << grid.x_init << std::endl;

      return 0;

}
