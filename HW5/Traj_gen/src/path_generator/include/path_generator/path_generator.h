#ifndef _PATH_GENERATOR_H_
#define _PATH_GENERATOR_H_

#include "lbfgs.hpp"
#include <iostream>
#include <algorithm>
#include <vector>
#include <ros/ros.h>
#include <geometry_msgs/Point.h>
#include <geometry_msgs/PoseStamped.h>
#include <nav_msgs/Path.h>
#include <Eigen/Eigen>
#include <pcl/point_cloud.h>
#include <pcl/point_types.h>
#include <pcl/kdtree/kdtree_flann.h>
#include <pcl_conversions/pcl_conversions.h>
#include <sensor_msgs/PointCloud2.h>

#include <visualization_msgs/MarkerArray.h>
#include <visualization_msgs/Marker.h>

#include <path_generator/matplotlibcpp.h>

class PathGenerator{
public:
    PathGenerator(){}
    ~PathGenerator(){}
    void init(ros::NodeHandle &nh);

private:
    bool need_path_;
    int pt_cnt_ = 0, edge_num_;
    double amax_ = 5.0, vmax_ = 3.0, rho_ = 1, beta_ = 1e3, gamma_ = 1;
    ros::Subscriber start_goal_sub_;
    ros::Publisher path_pub_, map_pub_, cloud_pub_;
    ros::Timer exec_timer_;

    Eigen::Matrix2d start_goal_;
    Eigen::MatrixXd obs_type_, obs_info_, obs_half_space_, AD_, Ac_, Ad_, AE_, G_, h_, lambda_;
    Eigen::VectorXd s_, px_, py_, px_dot_, py_dot_, px_dot_dot_, py_dot_dot_, f_;

    std::vector<Eigen::MatrixXd> A1_, b1_, A2_, b2_, A3_, b3_, A41_, b41_, A42_, b42_, A43_, b43_, A44_, b44_, A51_, b51_, A52_, b52_, A53_, b53_, A54_, b54_, mu1_, mu2_, mu3_, mu41_, mu42_, mu43_, mu44_, mu51_, mu52_, mu53_, mu54_;

    visualization_msgs::Marker map_marker_;
    visualization_msgs::MarkerArray map_markers_;
    sensor_msgs::PointCloud2 globalMap_pcd_;

    void execCallback(const ros::TimerEvent &e);
    void StartGoalCallback(const geometry_msgs::PoseStamped::ConstPtr &msg);
    void VisObs(Eigen::MatrixXd obs_info);
    void VisObsHpre(Eigen::MatrixXd obs_info);
    void VisPath(Eigen::VectorXd x);
    void TOPP(Eigen::VectorXd x);
    int run(const int N);
    void getTrajMat(const int N, Eigen::MatrixXd &AD, Eigen::MatrixXd &Ac, Eigen::MatrixXd &Ad, Eigen::MatrixXd &AE);
    double pt2line(Eigen::Vector2d pt, Eigen::Vector2d line1, Eigen::Vector2d line2, Eigen::Vector2d &ob){
        Eigen::Vector2d ap, ab;
        ap(0) = pt(0) - line1(0);
        ap(1) = pt(1) - line1(1);
        ab(0) = line2(0) - line1(0);
        ab(1) = line2(1) - line1(1);
        double r = ap.dot(ab) / ab.norm(), dis;
        if(r <= 0){
            dis = ap.norm();
            ob = line1;
        }else if (r >= 1){
            Eigen::Vector2d bp(pt(0) - line2(0), pt(1) - line2(1));
            dis = bp.norm();
            ob = line2;
        }else{
            Eigen::Vector2d c(line1(0) + ab(0) * r, line1(1) + ab(1) * r);
            Eigen::Vector2d cp(pt(0) - c(0), pt(1) - c(1));
            dis = cp.norm();
            ob = c;
        }
        return dis;
    }

    bool inside_obj(Eigen::Vector2d pt, int obj_idx){
        double temp;
        for(int i = 0; i < edge_num_; ++i){
            temp = obs_half_space_(obj_idx, i * 3) * pt(0) + obs_half_space_(obj_idx, i * 3 + 1) * pt(1) - obs_half_space_(obj_idx, i * 3 + 2);
            if(temp > 0)
                return false; 
        }
        return true;
    }

    Eigen::VectorXd proj2SOC(Eigen::VectorXd v){
        int n = v.size();
        Eigen::VectorXd proj;
        if(n >= 2){
            double v0 = v(0);
            int N = v.size();
            Eigen::VectorXd v1 = v.segment(1, N - 1);
            double v1norm = v1.norm();
            if(v0 <= -v1norm){
                proj = Eigen::VectorXd::Zero(N);
            }else if (v0 >= v1norm){
                proj = v;
            }else{
                double coeff = (v0 + v1norm) / 2 / v1norm;
                proj = v * coeff;
                proj(0) = v1norm * coeff;
            }
        }else{
            if(v(0) >= 0)
                proj = v;
            else
                proj = Eigen::VectorXd::Zero(1);
        }
        return proj;
    }

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
        double obs_num = obj.obs_info_.rows();
        Eigen::Vector2d pt, vec;
        // std::cout << obj.obs_half_space_ << std::endl;
        for(int i = 0; i < mid_pt_num; ++i){
            if(i > 0){
                pt(0) = x(i);
                pt(1) = x(i + mid_pt_num);
                for(int j = 0; j < obs_num; ++j){
                    if(obj.inside_obj(pt, j)){
                        
                        // double min_dis = 1e8, dis;
                        // Eigen::Vector2d min_ob, ob;
                        // for(int k = 0; k < obj.edge_num_; ++k){
                        //     int k_1 = (k + 1) % obj.edge_num_;
                        //     Eigen::Vector2d line1(obj.obs_info_(j, 2 * k_1), obj.obs_info_(j, 2 * k_1 + 1));
                        //     Eigen::Vector2d line2(obj.obs_info_(j, 2 * k), obj.obs_info_(j, 2 * k + 1));
                        //     dis = obj.pt2line(pt, line1, line2, ob);

                        //     if(dis < min_dis){
                        //         min_dis = dis;
                        //         min_ob = ob;
                        //     }
                        // }
                        // potential += min_dis;
                        // Eigen::Vector2d vec = pt - min_ob;
                        // g(i)              += 1000 * vec(0) / vec.norm();
                        // g(i + mid_pt_num) += 1000 * vec(1) / vec.norm();

                        double dis, epsilon = 1e-2, potential_ij = 0;
                        Eigen::Vector2d temp_g = Eigen::Vector2d::Zero(), ob, vec;
                        for(int k = 0; k < obj.edge_num_; ++k){
                            int k_1 = (k + 1) % obj.edge_num_;
                            Eigen::Vector2d line1(obj.obs_info_(j, 2 * k_1), obj.obs_info_(j, 2 * k_1 + 1));
                            Eigen::Vector2d line2(obj.obs_info_(j, 2 * k), obj.obs_info_(j, 2 * k + 1));
                            dis = obj.pt2line(pt, line1, line2, ob);

                            potential_ij += exp(-dis / epsilon);

                            vec = pt - ob;
                            temp_g(0) += exp(-dis / epsilon) * vec(0) / vec.norm();
                            temp_g(1) += exp(-dis / epsilon) * vec(1) / vec.norm();
                        }
                        potential += -epsilon * log(potential_ij);
                        g(i)              += 1000 * temp_g(0) / potential_ij;
                        g(i + mid_pt_num) += 1000 * temp_g(1) / potential_ij;

                        break;
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

    static double costFunction_TOPP(void *instance, const Eigen::VectorXd &x, Eigen::VectorXd &g){
        PathGenerator &obj = *(PathGenerator *)instance;
        // g = 
        double cost = (obj.G_ * x - obj.h_ + obj.lambda_ / obj.rho_).squaredNorm();
        g = obj.f_ + obj.G_.transpose() * (obj.lambda_ + obj.rho_ * (obj.G_ * x - obj.h_));
        int K = obj.A1_.size();
        for(int k = 0; k <= K; ++k){
            // std::cout << "1===================" << std::endl;
            if(k <= K - 1){
                cost += (obj.proj2SOC(obj.mu1_[k] / obj.rho_ - obj.A1_[k] * x - obj.b1_[k])).squaredNorm();
                cost += (obj.proj2SOC(obj.mu51_[k] / obj.rho_ - obj.A51_[k] * x - obj.b51_[k])).squaredNorm();
                cost += (obj.proj2SOC(obj.mu52_[k] / obj.rho_ - obj.A52_[k] * x - obj.b52_[k])).squaredNorm();
                cost += (obj.proj2SOC(obj.mu53_[k] / obj.rho_ - obj.A53_[k] * x - obj.b53_[k])).squaredNorm();
                cost += (obj.proj2SOC(obj.mu54_[k] / obj.rho_ - obj.A54_[k] * x - obj.b54_[k])).squaredNorm();
                g -= obj.A1_[k].transpose() * obj.proj2SOC(obj.mu1_[k] - obj.rho_ * (obj.A1_[k] * x + obj.b1_[k]));
                g -= obj.A51_[k].transpose() * obj.proj2SOC(obj.mu51_[k] - obj.rho_ * (obj.A51_[k] * x + obj.b51_[k]));
                g -= obj.A51_[k].transpose() * obj.proj2SOC(obj.mu52_[k] - obj.rho_ * (obj.A52_[k] * x + obj.b52_[k]));
                g -= obj.A51_[k].transpose() * obj.proj2SOC(obj.mu53_[k] - obj.rho_ * (obj.A53_[k] * x + obj.b53_[k]));
                g -= obj.A51_[k].transpose() * obj.proj2SOC(obj.mu54_[k] - obj.rho_ * (obj.A54_[k] * x + obj.b54_[k]));
            }
            // std::cout << "2===================" << std::endl;
            cost += (obj.proj2SOC(obj.mu2_[k] / obj.rho_ - obj.A2_[k] * x - obj.b2_[k])).squaredNorm();
            cost += (obj.proj2SOC(obj.mu3_[k] / obj.rho_ - obj.A3_[k] * x - obj.b3_[k])).squaredNorm();
            cost += (obj.proj2SOC(obj.mu41_[k] / obj.rho_ - obj.A41_[k] * x - obj.b41_[k])).squaredNorm();
            cost += (obj.proj2SOC(obj.mu42_[k] / obj.rho_ - obj.A42_[k] * x - obj.b42_[k])).squaredNorm();
            cost += (obj.proj2SOC(obj.mu43_[k] / obj.rho_ - obj.A43_[k] * x - obj.b43_[k])).squaredNorm();
            cost += (obj.proj2SOC(obj.mu44_[k] / obj.rho_ - obj.A44_[k] * x - obj.b44_[k])).squaredNorm();
            // std::cout << "3===================" << std::endl;
            g -= obj.A2_[k].transpose() * obj.proj2SOC(obj.mu2_[k] - obj.rho_ * (obj.A2_[k] * x + obj.b2_[k]));
            g -= obj.A3_[k].transpose() * obj.proj2SOC(obj.mu3_[k] - obj.rho_ * (obj.A3_[k] * x + obj.b3_[k]));
            g -= obj.A41_[k].transpose() * obj.proj2SOC(obj.mu41_[k] - obj.rho_ * (obj.A41_[k] * x + obj.b41_[k]));
            g -= obj.A42_[k].transpose() * obj.proj2SOC(obj.mu42_[k] - obj.rho_ * (obj.A42_[k] * x + obj.b42_[k]));
            g -= obj.A43_[k].transpose() * obj.proj2SOC(obj.mu43_[k] - obj.rho_ * (obj.A43_[k] * x + obj.b43_[k]));
            g -= obj.A44_[k].transpose() * obj.proj2SOC(obj.mu44_[k] - obj.rho_ * (obj.A44_[k] * x + obj.b44_[k]));
        }
        return obj.f_.dot(x) + obj.rho_ / 2 * cost;
    }

};

#endif