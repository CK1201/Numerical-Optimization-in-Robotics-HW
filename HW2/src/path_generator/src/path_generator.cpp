#include "path_generator/path_generator.h"

void PathGenerator::init(ros::NodeHandle &nh){
    need_path_ = false;
    
    start_goal_sub_ = nh.subscribe<geometry_msgs::PoseStamped>("/move_base_simple/goal", 1, &PathGenerator::StartGoalCallback, this);
    
    map_pub_ = nh.advertise<visualization_msgs::MarkerArray>("/path_generator/map", 1);
    path_pub_ = nh.advertise<visualization_msgs::Marker>("/path_generator/path", 1);
    // path_pub_ = nh.advertise<nav_msgs::Path>("/path_generator/path", 1);
    
    exec_timer_ = nh.createTimer(ros::Duration(0.1), &PathGenerator::execCallback, this);

    // map
    obs_info_ = Eigen::MatrixXd::Zero(10, 3);
    // FIXME obstacle
    obs_info_ << 0, 0, 5, // x,y,radius
                 9, 0, 3,
                 0, -8, 2,
                 3, 15, 5,
                 15, 10, 6,
                 -15, -2, 3,
                 -10, -12, 4,
                 6, -18, 3,
                 15, -5, 4,
                 -15, 10, 4;
    VisObs(obs_info_);
}

void PathGenerator::execCallback(const ros::TimerEvent &e){
    if(need_path_){
        std::cout << start_goal_ << std::endl;
        double max_vel = 1.0;
        double t1 = start_goal_(0, 0) - start_goal_(1, 0);
        double t2 = start_goal_(0, 1) - start_goal_(1, 1);
        int mid_pt_num = sqrt(t1 * t1 + t2 * t2) / max_vel;
        
        int ret = run(2 * mid_pt_num);
        need_path_ = false;
        std::cout << ret << std::endl;
    }
    map_pub_.publish(map_markers_);

}

void PathGenerator::StartGoalCallback(const geometry_msgs::PoseStamped::ConstPtr &msg){
    if(need_path_)
        return;
    start_goal_(pt_cnt_, 0) = msg->pose.position.x;
    start_goal_(pt_cnt_, 1) = msg->pose.position.y;
    ++pt_cnt_;
    pt_cnt_ %= 2;
    if(pt_cnt_ == 0)
        need_path_ = true;
}

int PathGenerator::run(const int N){
    double finalCost;
    Eigen::VectorXd x(N);
    const int mid_pt_num = N / 2;
    // start_goal_static_ = start_goal_;
    getTrajMat(mid_pt_num + 1, AD_, Ac_, Ad_, AE_);

    /* Set the initial guess */
    for (int i = 0; i < mid_pt_num; ++i){
        x(i)     = (start_goal_(1, 0) - start_goal_(0, 0)) / (mid_pt_num + 2) * (i + 1) + start_goal_(0, 0);
        x(i + mid_pt_num) = (start_goal_(1, 1) - start_goal_(0, 1)) / (mid_pt_num + 2) * (i + 1) + start_goal_(0, 1);
    }
    
    /* Set the minimization parameters */
    lbfgs::lbfgs_parameter_t params;
    params.g_epsilon = 1.0e-8;
    params.past = 3;
    params.delta = 1.0e-8;

    /* Start minimization */
    int ret = lbfgs::lbfgs_optimize(x,
                                    finalCost,
                                    costFunction,
                                    monitorProgress,
                                    this,
                                    params);

    VisPath(x);

    /* Report the result. */
    std::cout << std::setprecision(4)
                << "================================" << std::endl
                << "L-BFGS Optimization Returned: " << ret << std::endl
                << "Minimized Cost: " << finalCost << std::endl
                << "Optimal Variables: " << std::endl
                << x.transpose() << std::endl;
    return ret;
}

void PathGenerator::getTrajMat(const int N, Eigen::MatrixXd &AD, Eigen::MatrixXd &Ac, Eigen::MatrixXd &Ad, Eigen::MatrixXd &AE){
    Eigen::MatrixXd temp_D1, temp_D2, temp_c1, temp_c2, temp_d1, temp_d2, Dd;
        
    temp_D1 = Eigen::MatrixXd::Zero(N + 1, N - 1);
    temp_D2 = Eigen::MatrixXd::Zero(N - 1, N + 1);
    temp_c1 = Eigen::MatrixXd::Zero(N, N + 1);
    temp_c2 = Eigen::MatrixXd::Zero(N, N + 1);
    temp_d1 = Eigen::MatrixXd::Zero(N, N + 1);
    temp_d2 = Eigen::MatrixXd::Zero(N, N + 1);

    Dd = Eigen::MatrixXd::Identity(N - 1, N - 1) * 4;
    for(int i = 0; i < N; ++i){
        if(i < N - 2){
            Dd(i, i + 1) = 1;
            Dd(i + 1, i) = 1;
            temp_D2(i, i) = -1;
            temp_D2(i, i + 2) = 1;
        }
        temp_c1(i, i) = -1;
        temp_c1(i, i + 1) = 1;
        temp_c2(i, i) = -2;
        temp_c2(i, i + 1) = -1;

        temp_d1(i, i) = 1;
        temp_d1(i, i + 1) = -1;
        temp_d2(i, i) = 1;
        temp_d2(i, i + 1) = 1;
    }
    temp_D2(N - 2, N - 2) = -1;
    temp_D2(N - 2, N) = 1;
    
    temp_D1.block(1, 0, N - 1, N - 1) = Dd.inverse();
    AD = 3 * temp_D1 * temp_D2;

    Ac = 3 * temp_c1 + temp_c2 * AD;
    Ad = 2 * temp_d1 + temp_d2 * AD;
    AE = 4 * Ac.transpose() * Ac + 12 * Ac.transpose() * Ad + 12 * Ad.transpose() * Ad;
}

void PathGenerator::VisObs(Eigen::MatrixXd obs_info){
    int obs_size = obs_info.rows();
    map_marker_.header.frame_id = "map";
    map_marker_.ns = "path_generator";
    map_marker_.type = visualization_msgs::Marker::CYLINDER;
    map_marker_.action = visualization_msgs::Marker::ADD;
    map_marker_.pose.position.z = 1;
    map_marker_.pose.orientation.x = 0.0;
    map_marker_.pose.orientation.y = 0.0;
    map_marker_.pose.orientation.z = 0.0;
    map_marker_.pose.orientation.w = 1.0;
    map_marker_.scale.z = 2.00;
    map_marker_.color.a = 1.0;
    map_marker_.color.r = 1.0;
    map_marker_.color.g = 0.0;
    map_marker_.color.b = 0.0;

    for(int i = 0; i < obs_size; ++i){
        map_marker_.id = i;
        map_marker_.pose.position.x = obs_info(i, 0);
        map_marker_.pose.position.y = obs_info(i, 1);
        map_marker_.scale.x = obs_info(i, 2) * 2;
        map_marker_.scale.y = obs_info(i, 2) * 2;
        map_markers_.markers.push_back(map_marker_);
    }
    map_pub_.publish(map_markers_);
}

void PathGenerator::VisPath(Eigen::VectorXd x){
    const int n = x.size();
    const int mid_pt_num = n / 2;
    const int N = mid_pt_num + 1;
    Eigen::VectorXd all_x = Eigen::VectorXd::Zero(N + 1), all_y = Eigen::VectorXd::Zero(N + 1);
    all_x.segment(1, mid_pt_num) = x.segment(0, mid_pt_num);
    all_y.segment(1, mid_pt_num) = x.segment(mid_pt_num, mid_pt_num);

    all_x(0) = start_goal_(0, 0);
    all_y(0) = start_goal_(0, 1);
    all_x(N) = start_goal_(1, 0);
    all_y(N) = start_goal_(1, 1);

    // Eigen::MatrixXd AD, Ac, Ad, AE, D;
    Eigen::VectorXd a_x, b_x, c_x, d_x, a_y, b_y, c_y, d_y;
            
    // getTrajMat(N, AD, Ac, Ad, AE);
    a_x = all_x.segment(0, N);
    b_x = (AD_ * all_x).segment(0, N);
    c_x = Ac_ * all_x;
    d_x = Ad_ * all_x;
    a_y = all_y.segment(0, N);
    b_y = (AD_ * all_y).segment(0, N);
    c_y = Ac_ * all_y;
    d_y = Ad_ * all_y;

    int resolution = 10;
    double s;
    // nav_msgs::Path path;
    visualization_msgs::Marker path;
    path.header.stamp = ros::Time::now();
    path.header.frame_id = "map";
    path.ns = "codo_fsm_node/history_traj";
    path.id = 0;
    path.type = visualization_msgs::Marker::SPHERE_LIST;
    path.action = visualization_msgs::Marker::ADD;
    path.scale.x = 0.5;
    path.scale.y = 0.5;
    path.scale.z = 0.5;
    path.pose.orientation.w = 1.0;

    path.color.a = 1.0;
    path.color.r = 0.0;
    path.color.g = 1.0;
    path.color.b = 1.0;


    geometry_msgs::Point pt;
    pt.z = 1;

    for(int i = 0; i < N; ++i){
        for(int j = 0; j < resolution; ++j){
            s = j / resolution;
            pt.x = a_x(i) + b_x(i) * s + c_x(i) * s * s + d_x(i) * s * s * s;
            pt.y = a_y(i) + b_y(i) * s + c_y(i) * s * s + d_y(i) * s * s * s;
            path.points.push_back(pt);
        }
    }
    path_pub_.publish(path);
}