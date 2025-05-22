#include "acrobot_ilqr.h"

acrobot_ilqr::acrobot_ilqr(){
    Nx_ = 4;
    Nu_ = 1;
    Tfinal_ = 10.0;
    h_ = 0.05;
    Nt_ = (int)Tfinal_/h_+1;
    Q_ <<   1.0,  0,  0,  0,
            0,  1.0,    0,  0,
            0,  0,  0.1,    0,
            0,  0,  0,  0.1;

    Qn_ <<   100.0,  0,  0,  0,
            0,  100.0,    0,  0,
            0,  0,  100.0,    0,
            0,  0,  0,  100.0;

    R_ = 0.01;
    J_ = 0.0;
    xgoal_ << PI, 0.0, 0.0, 0.0;
    xtraj_ = new Eigen::Vector4d[Nt_];
    utraj_ = new double[Nt_-1];

    A_ = new Eigen::Matrix4d[Nt_-1];
    B_ = new Eigen::Vector4d[Nt_-1];
    K_ = new Eigen::Vector4d[Nt_-1];

}

acrobot_ilqr::~acrobot_ilqr(){
    delete[] xtraj_;
    delete[] utraj_;
    delete[] A_;
    delete[] B_;
    delete[] K_;
}

template <class Type>
Type acrobot_ilqr::dynamics( Type ax)     //由于要使用CppAD导致dynamics参数有点奇怪，ax前四个元素为x，第五个元素为u
{   
    // parameters
    Type x_dot(4);
    const static double g = 9.8;
    const static double l1 = 1.0;
    const static double lc1 = 0.5;
    const static double m1 = 1.0;
    const static double I1 = m1 * l1 * l1 / 12.0;
    const static double l2 = 2.0;
    const static double lc2 = 1.0;
    const static double m2 = 1.0;
    const static double I2 = m2 * l2 * l2 / 12.0;

    // intermediate state values
    CppAD::AD<double> d1 = m1*lc1*lc1 + m2*(l1*l1 + lc2*lc2 + 2*l1*lc2*cos(ax[THETA2])) + I1 + I2;
    CppAD::AD<double> d2 = m2*(lc2*lc2 + l1*lc2*cos(ax[THETA2])) + I2;
    CppAD::AD<double> phi2 = m2*lc2*g*cos(ax[THETA1] + ax[THETA2] - M_PI/2);
    CppAD::AD<double> phi1 = -m2*l1*lc2*ax[THETA2DOT]*ax[THETA2DOT]*sin(ax[THETA2]) -
                  2*m2*l1*lc2*ax[THETA2DOT]*ax[THETA1DOT]*sin(ax[THETA2]) +
                  (m1*lc1 + m2*l1)*g*cos(ax[THETA1]-M_PI/2) + phi2;

    // dynamics
    x_dot[THETA1]    = ax[THETA1DOT];
    x_dot[THETA2]    = ax[THETA2DOT];
    x_dot[THETA2DOT] = (ax[CONTROL] + phi1*d2/d1 - m2*l1*lc2*ax[THETA1DOT]*ax[THETA1DOT]*sin(ax[THETA2]) - phi2)/(m2*lc2*lc2 + I2 - d2*d2/d1);
    x_dot[THETA1DOT] = -(d2*x_dot[THETA2DOT] + phi1)/d1;
    // std::cout << "xdot = " << x_dot[0] << " "<< x_dot[1] << " "<< x_dot[2] << " "<< x_dot[3] << std::endl;
    return x_dot;
}

Eigen::Vector4d acrobot_ilqr::dynamics(Eigen::Vector4d x, const double &u)
{
    // parameters
    const static double g = 9.8;
    const static double l1 = 1.0;
    const static double lc1 = 0.5;
    const static double m1 = 1.0;
    const static double I1 = m1 * l1 * l1 / 12.0;
    const static double l2 = 2.0;
    const static double lc2 = 1.0;
    const static double m2 = 1.0;
    const static double I2 = m2 * l2 * l2 / 12.0;

    // intermediate state values
    double d1 = m1*lc1*lc1 + m2*(l1*l1 + lc2*lc2 + 2*l1*lc2*cos(x(THETA2))) + I1 + I2;
    double d2 = m2*(lc2*lc2 + l1*lc2*cos(x(THETA2))) + I2;
    double phi2 = m2*lc2*g*cos(x(THETA1) + x(THETA2) - M_PI/2);
    double phi1 = -m2*l1*lc2*x(THETA2DOT)*x(THETA2DOT)*sin(x(THETA2)) -
                    2*m2*l1*lc2*x(THETA2DOT)*x(THETA1DOT)*sin(x(THETA2)) +
                    (m1*lc1 + m2*l1)*g*cos(x(THETA1)-M_PI/2) + phi2;

    // dynamics
    Eigen::Vector4d x_dot;
    x_dot(THETA1)    = x(THETA1DOT);
    x_dot(THETA2)    = x(THETA2DOT);
    x_dot(THETA2DOT) = (u + phi1*d2/d1 - m2*l1*lc2*x(THETA1DOT)*x(THETA1DOT)*sin(x(THETA2)) - phi2)/(m2*lc2*lc2 + I2 - d2*d2/d1);
    x_dot(THETA1DOT) = -(d2*x_dot(THETA2DOT) + phi1)/d1;
    return x_dot;
}

Eigen::Vector4d acrobot_ilqr::dynamics_RK4(Eigen::Vector4d x, const double &u){
    Eigen::Vector4d f1,f2,f3,f4,xnext;

    f1 = dynamics(x,u);
    f2 = dynamics(x+0.5*h_*f1,u);
    f3 = dynamics(x+0.5*h_*f2,u);
    f4 = dynamics(x+h_*f3,u);

    xnext = x + (h_/6.0)*(f1+2*f2+2*f3+f4);

    return xnext;
}


template <class Type>
Type acrobot_ilqr::dynamics_RK4( Type ax){

    // wrap angle 
    // while (ax[0] > PI || ax[0] < -PI ){
    //     ax[0] = (ax[0] > PI) ? (ax[0]-2*PI) : ax[0];
    //     ax[0] = (ax[0] < -PI) ? (ax[0]+2*PI) : ax[0];   
    // }
    // while (ax[1] > PI || ax[1] < -PI ){
    //     ax[1] = (ax[1] > PI) ? (ax[1]-2*PI) : ax[1];
    //     ax[1] = (ax[1] < -PI) ? (ax[1]+2*PI) : ax[1];
    // }
    
    Type f1(4),f2(4),f3(4),f4(4),ax_temp(5);
    ax_temp = ax;
    f1 = dynamics(ax);
    for(int i =0;i<4;i++)
        ax_temp[i] = ax[i]+0.5*h_*f1[i];
    f2 = dynamics(ax_temp);

    for(int i =0;i<4;i++)
        ax_temp[i] = ax[i]+0.5*h_*f2[i];
    f3 = dynamics(ax_temp);

    for(int i =0;i<4;i++)
        ax_temp[i] = ax[i]+h_*f3[i];
    f4 = dynamics(ax_temp);

    Type xnext(4);
    for (int i=0;i<4;i++)
        xnext[i] = ax[i] + (h_/6.0)*(f1[i]+2*f2[i]+2*f3[i]+f4[i]);

    // wrap angle 
    // while (xnext[0] > PI || xnext[0] < -PI ){
    //     xnext[0] = (xnext[0] > PI) ? (xnext[0]-2*PI) : xnext[0];
    //     xnext[0] = (xnext[0] < -PI) ? (xnext[0]+2*PI) : xnext[0];   
    // }
    // while (xnext[1] > PI || xnext[1] < -PI ){
    //     xnext[1] = (xnext[1] > PI) ? (xnext[1]-2*PI) : xnext[1];
    //     xnext[1] = (xnext[1] < -PI) ? (xnext[1]+2*PI) : xnext[1];
    // }


    // std::cout << "xnext = " << xnext[0] << " "<< xnext[1] << " "<< xnext[2] << " "<< xnext[3] << std::endl;
    
    return xnext;
}

Eigen::MatrixXd acrobot_ilqr::dfdx(const Eigen::Vector4d &x, const double &u){
    using CppAD::AD;
    using std::vector;

    // domain space vector
     size_t n = 5;               // number of domain space variables
     vector< AD<double> > ax(5); // vector of domain space variables

    // declare independent variables and start recording operation sequence
    CppAD::Independent(ax);

    // range space vector
     size_t m = 4;               // number of ranges space variables
     vector< AD<double> > ay(4); // vector of ranges space variables

    ay = dynamics_RK4(ax);     
    // store operation sequence in f: X -> Y and stop recording
     CppAD::ADFun<double> f(ax, ay);

    // compute derivative using operation sequence stored in f
     vector<double> jac(m * n),hes(n*n); // Jacobian of f (m by n matrix)
     vector<double> xx(n);       // domain space vector

    xx[0] = x(0);
    xx[1] = x(1);
    xx[2] = x(2);
    xx[3] = x(3);
    xx[4] = u;
    jac  = f.Jacobian(xx);      // Jacobian for operation sequence

    Eigen::MatrixXd Jac_xu(4,5);
    for (int i = 0;i<4;i++)
        for (int j = 0;j<5;j++)
            Jac_xu(i,j) = jac[5*i+j];
    return Jac_xu;
}

double acrobot_ilqr::stage_cost(Eigen::Vector4d x, double u){
    return 0.5*(x.transpose()-xgoal_.transpose())*Q_*(x-xgoal_)+0.5*R_*u*u;
}

double acrobot_ilqr::terminal_cost(Eigen::Vector4d x){
    return 0.5*(x.transpose()-xgoal_.transpose())*Qn_*(x-xgoal_);
}

double acrobot_ilqr::cost( Eigen::Vector4d *xtraj, double *utraj){
    double J = 0.0;
    for (int i = 0;i<(Nt_-2);i++)
        J = J + stage_cost(xtraj[i],utraj[i]);
    J = J + terminal_cost(xtraj[Nt_-1]);
    return J;
}

void acrobot_ilqr::run_ilqr(){
    Eigen::Vector4d x0;
    x0 << 0.0, 0.0, 0.0, 0.0;
    // double u0 =1.0;
    xtraj_[0] = x0;
    for (int i = 0;i<(Nt_-1);i++)
        utraj_[i] = 0.1;

    // initial rollout 
    for (int i = 0;i<Nt_-1;i++)
        xtraj_[i+1] = dynamics_RK4(xtraj_[i],utraj_[i]);
    J_ = cost(xtraj_,utraj_);

    // DDP algorithm
    Eigen::Vector4d *p = new Eigen::Vector4d[Nt_];
    Eigen::Matrix4d *P = new Eigen::Matrix4d[Nt_];
    double *d = new double[Nt_-1];
    double *d_abs = new double[Nt_-1];
    d_abs[0] = 1.0;

    double delta_J = 0.0;

    Eigen::Vector4d *xn = new Eigen::Vector4d[Nt_];
    double *un = new double[Nt_-1];

    Eigen::Vector4d gx,Gxu,Gux,q;
    Eigen::Matrix4d Gxx;
    double gu=0.0,Guu=0.0,r;

    Eigen::MatrixXd derivative(4,5);

    double alpha ,Jn;

    int iter = 0, iter_max = 150;

    std::ofstream write_data("../data/acrobot_ilqr/iteration",std::ios::out);

    while ( *std::max_element(d_abs,d_abs+Nt_-1) > 1e-2 ){
        iter += 1;
        delta_J = 0.0;

        p[Nt_-1] = Qn_*(xtraj_[Nt_-1]-xgoal_);
        P[Nt_-1] = Qn_;

        // backward pass
        for (int i = (Nt_-2);i>-1;i--){
            q = Q_*(xtraj_[i]-xgoal_);
            r = R_*utraj_[i];
            

            derivative = dfdx(xtraj_[i],utraj_[i]);
            A_[i] = derivative.leftCols(4);
            B_[i] = derivative.rightCols(1);

            gx = q + A_[i].transpose() * p[i+1];
            gu = r + B_[i].transpose() * p[i+1];

            Gxx = Q_ + A_[i].transpose()*P[i+1]*A_[i];
            Guu = R_ + B_[i].transpose()*P[i+1]*B_[i];
            Gxu = A_[i].transpose()*P[i+1]*B_[i];
            Gux = B_[i].transpose()*P[i+1]*A_[i];

            d[i] = gu/Guu;
            d_abs[i] = abs(d[i]);
            K_[i] = Gux/Guu;

            p[i] = gx - K_[i]*gu+K_[i]*Guu*d[i] - Gxu*d[i];
            P[i] = Gxx + K_[i]*Guu*K_[i].transpose() - Gxu*K_[i].transpose() - K_[i]*Gux.transpose();

            delta_J = delta_J + gu * d[i];
        }
        // forward rollout with line search
        xn[0] = xtraj_[0];
        alpha = 1.0;

        for (int i = 0;i<Nt_-1;i++){
            xn[i+1] = dynamics_RK4(xn[i],un[i]);
        }
            Jn = cost(xn,un);
        
        while (std::isnan(Jn) || Jn > (J_ - 1e-2*alpha*delta_J)){
            alpha = 0.5*alpha;

            for (int i = 0;i<Nt_-1;i++){
                un[i] = utraj_[i] - alpha*d[i] - K_[i].transpose()*(xn[i] - xtraj_[i]);
                xn[i+1] = dynamics_RK4(xn[i],un[i]);
            }
            Jn = cost(xn,un);
        }

        J_ = Jn;
        for (int i =0;i<Nt_-1;i++){
            xtraj_[i] = xn[i];
            utraj_[i] = un[i];
        }
        xtraj_[Nt_-1] = xn[Nt_-1];

        std::cout << "iter = " << iter <<"\t delta_J = "<< delta_J << "\t J = " << J_ << "\t max element = " <<*std::max_element(d_abs,d_abs+Nt_-1) << std::endl;

        write_data << iter << " " << J_ << " " <<*std::max_element(d_abs,d_abs+Nt_-1) << std::endl;
        
    }

    delete[] p,P,d,d_abs,xn,un;
}

void acrobot_ilqr::plot(){
    // namespace plt = matplotlibcpp;    
    // std::vector<double> theta1(Nt_),theta2(Nt_),theta1_dot(Nt_),theta2_dot(Nt_),u(Nt_-1);

    // for (int i = 0;i<Nt_-1;i++){
    //     theta1[i] = xtraj_[i](0);
    //     theta2[i] = xtraj_[i](1);
    //     theta1_dot[i] = xtraj_[i](2);
    //     theta2_dot[i] = xtraj_[i](3);
    //     u[i] = utraj_[i];
    // }

    // theta1[Nt_-1] = xtraj_[Nt_-1](0);
    // theta2[Nt_-1] = xtraj_[Nt_-1](1);
    // theta1_dot[Nt_-1] = xtraj_[Nt_-1](2);
    // theta2_dot[Nt_-1] = xtraj_[Nt_-1](3);

    // // 绘图
    // plt::subplot(2, 1, 1);  // 2行1列，第1个子图
    // plt::plot(theta1, {{"label", "theta1"}});
    // plt::plot(theta2, {{"label", "theta2"}});
    // plt::plot(theta1_dot, {{"label", "theta1_dot"}});
    // plt::plot(theta2_dot, {{"label", "theta2_dot"}});
    // plt::title("State");  
    // plt::legend();  


    // plt::subplot(2, 1, 2);  
    // plt::plot(u, {{"label", "u"}});

    // // 添加标题和图例
    // plt::title("Control");  // 子图标题
    // plt::legend();  // 显示图例

    // // 显示图形
    // plt::show();
}

void acrobot_ilqr::record(){
    std::filesystem::create_directories("../data/acrobot_ilqr/");
    std::ofstream write_traj("../data/acrobot_ilqr/traj",std::ios::out);
    std::ofstream write_A("../data/acrobot_ilqr/A",std::ios::out);
    std::ofstream write_B("../data/acrobot_ilqr/B",std::ios::out);
    std::ofstream write_K("../data/acrobot_ilqr/K",std::ios::out);
    for (int i = 0;i<Nt_-1;i++){
        write_traj << xtraj_[i](0) << " "<< xtraj_[i](1) << " "<< xtraj_[i](2) << " "<< xtraj_[i](3) << " " << utraj_[i]<<std::endl;

        write_A << A_[i](0,0) << " " << A_[i](0,1) << " " << A_[i](0,2) << " " << A_[i](0,3) << " " << A_[i](1,0) << " " << A_[i](1,1) << " " << A_[i](1,2) << " " 
        << A_[i](1,3) << " " << A_[i](2,0) << " " << A_[i](2,1) << " " << A_[i](2,2) << " " << A_[i](2,3) << " " << A_[i](3,0) << " " << A_[i](3,1) << " " << 
        A_[i](3,2) << " " << A_[i](3,3) << std::endl;

        write_B << B_[i](0) << " " << B_[i](1) << " " << B_[i](2) << " " << B_[i](3) << std::endl;
        
        write_K << K_[i](0) << " " << K_[i](1) << " " << K_[i](2) << " " << K_[i](3) << std::endl;
         

    }
    write_traj << xtraj_[Nt_-1](0) << " "<< xtraj_[Nt_-1](1) << " "<< xtraj_[Nt_-1](2) << " "<< xtraj_[Nt_-1](3) << " " << 0.0 << std::endl;

    write_traj.close();
    write_A.close();
    write_B.close();
    write_K.close();
}

double acrobot_ilqr::wrapAngle(double angle){
    while(angle > PI || angle < -PI)
    {
        if ( angle > PI )
            angle = angle - 2*PI;
        else
            angle = angle + 2*PI;
    }
    return angle;
}
