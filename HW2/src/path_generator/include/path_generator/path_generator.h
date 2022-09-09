#ifndef _PATH_GENERATOR_H_
#define _PATH_GENERATOR_H_

#include "lbfgs.hpp"
#include <iostream>
#include <algorithm>
#include <ros/ros.h>
#include <geometry_msgs/Point.h>
#include <geometry_msgs/PoseStamped.h>
#include <nav_msgs/Path.h>
#include <Eigen/Eigen>

#include <visualization_msgs/MarkerArray.h>
#include <visualization_msgs/Marker.h>

class PathGenerator{
public:
    PathGenerator(){}
    ~PathGenerator(){}
    void init(ros::NodeHandle &nh);

private:
    bool need_path_;
    int pt_cnt_ = 0;
    ros::Subscriber start_goal_sub_;
    ros::Publisher path_pub_, map_pub_;
    ros::Timer exec_timer_;

    Eigen::Matrix2d start_goal_;
    Eigen::MatrixXd obs_info_, AD_, Ac_, Ad_, AE_;

    visualization_msgs::Marker map_marker_;
    visualization_msgs::MarkerArray map_markers_;

    void execCallback(const ros::TimerEvent &e);
    void StartGoalCallback(const geometry_msgs::PoseStamped::ConstPtr &msg);
    void VisObs(Eigen::MatrixXd obs_info);
    void VisPath(Eigen::VectorXd x);
    int run(const int N);
    void getTrajMat(const int N, Eigen::MatrixXd &AD, Eigen::MatrixXd &Ac, Eigen::MatrixXd &Ad, Eigen::MatrixXd &AE);

    static double costFunction(void *instance, const Eigen::VectorXd &x, Eigen::VectorXd &g){
        PathGenerator &obj = *(PathGenerator *)instance;
        const int n = x.size();
        const int mid_pt_num = n / 2; // mid point num
        const int N = mid_pt_num + 1; // N of x_N
        g = Eigen::VectorXd::Zero(n);
        
        Eigen::VectorXd all_x = Eigen::VectorXd::Zero(N + 1), all_y = Eigen::VectorXd::Zero(N + 1);
        all_x(0) = obj.start_goal_(0, 0);
        all_y(0) = obj.start_goal_(0, 1);
        all_x(N) = obj.start_goal_(1, 0);
        all_y(N) = obj.start_goal_(1, 1);
        all_x.segment(1, mid_pt_num) = x.segment(0, mid_pt_num);
        all_y.segment(1, mid_pt_num) = x.segment(mid_pt_num, mid_pt_num);
        // std::cout << "3333333333" << std::endl;
        double fx = 0.0, potential = 0.0, energy = 0.0;
        // std::cout << all_x.size() << std::endl;
        // std::cout << AE.size() << std::endl;
        energy += all_x.transpose() * obj.AE_ * all_x;
        energy += all_y.transpose() * obj.AE_ * all_y;
        // std::cout << "4444444444" << std::endl;
        Eigen::VectorXd temp11 = (obj.AE_ + obj.AE_.transpose()) * all_x;
        // std::cout << "5555555555" << std::endl;
        g.segment(0, mid_pt_num) += temp11.segment(1, mid_pt_num);
        temp11 = (obj.AE_ + obj.AE_.transpose()) * all_y;
        g.segment(mid_pt_num, mid_pt_num) += temp11.segment(1, mid_pt_num);

        // std::cout << "666666666" << std::endl;
        double temp, obs_num = obj.obs_info_.rows();
        Eigen::Vector2d pt, vec;
        for(int i = 0; i < mid_pt_num; ++i){
            if(i > 0){
                pt(0) = x(i);
                pt(1) = x(i + mid_pt_num);
                for(int j = 0; j < obs_num; ++j){
                    vec(0) = pt(0) - obj.obs_info_(j, 0);
                    vec(1) = pt(1) - obj.obs_info_(j, 1);
                    temp = obj.obs_info_(j, 2) - vec.norm();
                    if(temp > 0){
                        potential += temp;
                        g(i)              -= 1000 * vec(0) / vec.norm();
                        g(i + mid_pt_num) -= 1000 * vec(1) / vec.norm();
                    }
                }
            }
        }
        // std::cout << "777777777" << std::endl;
        fx = 1000 * potential + energy;
        return fx;
    }

    static int monitorProgress(void *instance,
                               const Eigen::VectorXd &x,
                               const Eigen::VectorXd &g,
                               const double fx,
                               const double step,
                               const int k,
                               const int ls)
    {
        // std::cout << std::setprecision(4)
        //           << "================================" << std::endl
        //           << "Iteration: " << k << std::endl
        //           << "Function Value: " << fx << std::endl
        //           << "Gradient Inf Norm: " << g.cwiseAbs().maxCoeff() << std::endl
        //           << "Variables: " << std::endl
        //           << x.transpose() << std::endl;
        return 0;
    }

};

#endif