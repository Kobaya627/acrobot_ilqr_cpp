#ifndef SRC_ACROBOT_ILQR_H
#define SRC_ACROBOT_ILQR_H

#define THETA1  0
#define THETA2  1
#define THETA1DOT   2
#define THETA2DOT   3
#define CONTROL 4
#define PI  3.141592653

#include "Eigen/Core"
#include "vector"
#include "cppad/cppad.hpp"
// #include "matplotlib/matplotlibcpp.h"
#include "fstream"
#include <filesystem>
namespace fs = std::filesystem;

class acrobot_ilqr
{
private:
    /* data */
    int Nx_,Nu_,Nt_;
    double Tfinal_,h_,R_, *utraj_,J_;
    Eigen::Matrix4d Q_,Qn_;
    Eigen::Vector4d xgoal_, *xtraj_;
    Eigen::Matrix4d *A_;
    Eigen::Vector4d *B_;
    Eigen::Vector4d *K_;

    double stage_cost(Eigen::Vector4d x, double u);
    double terminal_cost(Eigen::Vector4d x);
    double cost(Eigen::Vector4d *xtraj, double *utraj);

    double wrapAngle(double angle);

public:

    template <class Type>
    Type dynamics( Type ax);

    Eigen::Vector4d dynamics(Eigen::Vector4d x, const double &u);

    template <class Type>
    Type dynamics_RK4( Type ax);

    Eigen::Vector4d dynamics_RK4( Eigen::Vector4d x, const double &u);


    Eigen::MatrixXd dfdx(const Eigen::Vector4d &x, const double &u);

    acrobot_ilqr(/* args */);
    ~acrobot_ilqr();
    void run_ilqr();
    void plot();
    void record();
    
};

#endif