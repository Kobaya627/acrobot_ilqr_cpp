#include "acrobot_ilqr.h"

int main(int argc, char **argv)
{
    acrobot_ilqr controller;
    Eigen::Vector4d x,xdot,xnext;
    // x << PI , 0.0, 0.0, 0.0 ;
    // double u = 0.0;
    // xdot = controller.dynamics(x,u);
    // xnext = controller.dynamics_RK4(x,u);
    // std::cout << "x = " << x.transpose() <<std::endl;
    // std::cout << "xdot = "<<xdot.transpose() <<std::endl;
    // std::cout << "xnext = "<<xnext.transpose() <<std::endl<<std::endl;

    // // std::vector<double> xx(5);       // domain space vector
    // // std::vector<double> yy(4);       // range space vector
    // Eigen::VectorXd xx(5),yy(4);

    // xx[0] = x(0);
    // xx[1] = x(1);
    // xx[2] = x(2);
    // xx[3] = x(3);
    // xx[4] = u;

    // yy = controller.dynamics(xx);
    // std::cout << "yy = " << yy[0] << " "<< yy[1] << " "<< yy[2] << " "<< yy[3] << std::endl;

    controller.run_ilqr();
    // controller.plot();   
    controller.record();

    return 0;
}
