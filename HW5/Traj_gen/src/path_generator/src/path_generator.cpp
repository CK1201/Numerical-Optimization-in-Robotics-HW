#include "path_generator/path_generator.h"

namespace plt = matplotlibcpp;

void PathGenerator::init(ros::NodeHandle &nh){
    need_path_ = false;
    
    start_goal_sub_ = nh.subscribe<geometry_msgs::PoseStamped>("/move_base_simple/goal", 1, &PathGenerator::StartGoalCallback, this);
    
    map_pub_ = nh.advertise<visualization_msgs::MarkerArray>("/path_generator/map", 1);
    cloud_pub_ = nh.advertise<sensor_msgs::PointCloud2>("/path_generator/global_cloud", 1);
    path_pub_ = nh.advertise<visualization_msgs::Marker>("/path_generator/path", 1);
    // path_pub_ = nh.advertise<nav_msgs::Path>("/path_generator/path", 1);
    
    exec_timer_ = nh.createTimer(ros::Duration(0.1), &PathGenerator::execCallback, this);

    // map
    srand(1);
    edge_num_ = 5;
    obs_type_ = Eigen::MatrixXd::Zero(10, 3);
    // FIXME obstacle
    obs_type_ << 0, 0, 5, // x,y,radius
                 9, 0, 3,
                 0, -8, 2,
                 3, 15, 5,
                 15, 10, 6,
                 -15, -2, 3,
                 -10, -12, 4,
                 6, -18, 3,
                 15, -5, 4,
                 -15, 10, 4;
    int delta_angle = 360 / 5;
    double angle, radius;
    // Eigen::Vector5d temp(0.2, 0.4, 0.6, 0.8, 1.0);
    obs_info_ = Eigen::MatrixXd::Zero(10, edge_num_ * 2);
    obs_half_space_ = Eigen::MatrixXd::Zero(10, edge_num_ * 3);
    for(int i = 0; i < obs_type_.rows(); ++i){
        // int temp = i % edge_num_;
        for(int j = 0; j < edge_num_; ++j){
            angle = (delta_angle * j + rand() % delta_angle) / 180.0 * M_PI;
            // radius = obs_type_(i, 2) * ((j + temp) % 5 * 0.2);
            radius = obs_type_(i, 2);
            obs_info_(i, 2 * j) = obs_type_(i, 0) + cos(angle) * radius;
            obs_info_(i, 2 * j + 1) = obs_type_(i, 1) + sin(angle) * radius;
        }

        for(int j = 0; j < edge_num_; ++j){
            int j_1 = (j + 1) % edge_num_;
            obs_half_space_(i, 3 * j) = -(obs_info_(i, 2 * j_1 + 1) - obs_info_(i, 2 * j + 1)) / (obs_info_(i, 2 * j_1) - obs_info_(i, 2 * j));
            obs_half_space_(i, 3 * j + 1) = 1;
            obs_half_space_(i, 3 * j + 2) = obs_info_(i, 2 * j + 1) + obs_half_space_(i, 3 * j) * obs_info_(i, 2 * j);
            if(obs_half_space_(i, 3 * j) * obs_type_(i, 0) + obs_half_space_(i, 3 * j + 1) * obs_type_(i, 1) > obs_half_space_(i, 3 * j + 2)){
                // std::cout << obs_half_space_ << std::endl;
                obs_half_space_(i, 3 * j) *= -1;
                obs_half_space_(i, 3 * j + 1) *= -1;
                obs_half_space_(i, 3 * j + 2) *= -1;
                // obs_half_space_.row(i) = obs_half_space_.row(i) * -1;
                
            }
        }
    }
    // std::cout << obs_half_space_ << std::endl;
    VisObsHpre(obs_info_);
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
    // map_pub_.publish(map_markers_);
    cloud_pub_.publish(globalMap_pcd_);
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

    TOPP(x);

    /* Report the result. */
    // std::cout << std::setprecision(4)
    //             << "================================" << std::endl
    //             << "L-BFGS Optimization Returned: " << ret << std::endl
    //             << "Minimized Cost: " << finalCost << std::endl
    //             << "Optimal Variables: " << std::endl
    //             << x.transpose() << std::endl;
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

void PathGenerator::VisObsHpre(Eigen::MatrixXd obs_info){
    int obs_size = obs_info.rows();
    pcl::PointCloud<pcl::PointXYZ> cloudMap;
    pcl::PointXYZ pt;
    
    for(int i = 0; i < obs_size; ++i){
        for(int j = 0; j < edge_num_; ++j){
            int j_1 = (j + 1) % edge_num_;
            double temp1 = obs_info(i, j_1 * 2) - obs_info(i, j * 2);
            double temp2 = obs_info(i, j_1 * 2 + 1) - obs_info(i, j * 2 + 1);
            double dis = sqrt(temp1 * temp1 + temp2 * temp2);
            int num = dis * 10;
            for(int k = 0; k < num; ++k){
                pt.x = obs_info(i, j * 2) + (obs_info(i, j_1 * 2) - obs_info(i, j * 2)) * k / num;
                pt.y = obs_info(i, j * 2 + 1) + (obs_info(i, j_1 * 2 + 1) - obs_info(i, j * 2 + 1)) * k / num;
                for(double z = 0; z < 2; z += 0.1){
                    pt.z = z;
                    cloudMap.points.push_back(pt);
                }
            }
        }
        
    }

    cloudMap.width = cloudMap.points.size();
    cloudMap.height = 1;
    cloudMap.is_dense = true;
    
    pcl::toROSMsg(cloudMap, globalMap_pcd_);
    globalMap_pcd_.header.frame_id = "map";
    cloud_pub_.publish(globalMap_pcd_);
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

    // int resolution = 10;
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
    path.color.r = 1.0;
    path.color.g = 0.0;
    path.color.b = 0.0;


    geometry_msgs::Point pt;
    pt.z = 1;

    for(int i = 0; i < N; ++i){
        for(int j = 0; j < 2; ++j){
            // s = double(j) / resolution;
            s = j;
            pt.x = a_x(i) + b_x(i) * s + c_x(i) * s * s + d_x(i) * s * s * s;
            pt.y = a_y(i) + b_y(i) * s + c_y(i) * s * s + d_y(i) * s * s * s;
            path.points.push_back(pt);
        }
    }
    pt.x = all_x(N);
    pt.y = all_y(N);
    path.points.push_back(pt);
    path_pub_.publish(path);
}

void PathGenerator::TOPP(Eigen::VectorXd x){
    const int n = x.size();
    const int mid_pt_num = n / 2;
    const int N = mid_pt_num + 1;
    const int resolution = 10;
    const int K = N * resolution;

    s_ = Eigen::VectorXd::Zero(K + 1);
    px_ = Eigen::VectorXd::Zero(K + 1);
    py_ = Eigen::VectorXd::Zero(K + 1);
    px_dot_ = Eigen::VectorXd::Zero(K + 1);
    py_dot_ = Eigen::VectorXd::Zero(K + 1);
    px_dot_dot_ = Eigen::VectorXd::Zero(K + 1);
    py_dot_dot_ = Eigen::VectorXd::Zero(K + 1);

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

    double p = 0;
    int cnt = 1;
    // double last_x, last_y, now_x, now_y;
    px_(0) = a_x(0) + b_x(0) * p + c_x(0) * p * p + d_x(0) * p * p * p;
    py_(0) = a_y(0) + b_y(0) * p + c_y(0) * p * p + d_y(0) * p * p * p;
    s_(0) = 0;

    for(int i = 0; i < N; ++i){
        for(int j = 1; j <= resolution; ++j){
            p = double(j) / resolution;
            px_(cnt) = a_x(i) + b_x(i) * p + c_x(i) * p * p + d_x(i) * p * p * p;
            py_(cnt) = a_y(i) + b_y(i) * p + c_y(i) * p * p + d_y(i) * p * p * p;
            s_(cnt) = s_(cnt - 1) + sqrt((px_(cnt) - px_(cnt - 1)) * (px_(cnt) - px_(cnt - 1)) + (py_(cnt) - py_(cnt - 1)) * (py_(cnt) - py_(cnt - 1)));
            
            // std::cout << sqrt((px_(cnt) - px_(cnt - 1)) * (px_(cnt) - px_(cnt - 1)) + (py_(cnt) - py_(cnt - 1)) * (py_(cnt) - py_(cnt - 1))) << ", " << px_(cnt) - px_(cnt - 1) << ", " << py_(cnt) - py_(cnt - 1) << ", " << p << std::endl;
            ++cnt;
        }

    }
    
    // std::cout << K << std::endl;

    for(int i = 0; i < K + 1; ++i){
        if(i == 0){
            px_dot_(i) = (px_(i + 1) - px_(i)) / (s_(i + 1) - s_(i));
            py_dot_(i) = (py_(i + 1) - py_(i)) / (s_(i + 1) - s_(i));
        }else if (i == K){
            px_dot_(i) = (px_(i) - px_(i - 1)) / (s_(i) - s_(i - 1));
            py_dot_(i) = (py_(i) - py_(i - 1)) / (s_(i) - s_(i - 1));
        }else{
            px_dot_(i) = (px_(i + 1) - px_(i - 1)) / (s_(i + 1) - s_(i - 1));
            py_dot_(i) = (py_(i + 1) - py_(i - 1)) / (s_(i + 1) - s_(i - 1));
        }
    }

    for(int i = 0; i < K + 1; ++i){
        if(i == 0){
            px_dot_dot_(i) = (px_dot_(i + 1) - px_dot_(i)) / (s_(i + 1) - s_(i));
            py_dot_dot_(i) = (py_dot_(i + 1) - py_dot_(i)) / (s_(i + 1) - s_(i));
        }else if (i == K){
            px_dot_dot_(i) = (px_dot_(i) - px_dot_(i - 1)) / (s_(i) - s_(i - 1));
            py_dot_dot_(i) = (py_dot_(i) - py_dot_(i - 1)) / (s_(i) - s_(i - 1));
        }else{
            px_dot_dot_(i) = (px_dot_(i + 1) - px_dot_(i - 1)) / (s_(i + 1) - s_(i - 1));
            py_dot_dot_(i) = (py_dot_(i + 1) - py_dot_(i - 1)) / (s_(i + 1) - s_(i - 1));
        }
    }

    // Conic ALM
    Eigen::MatrixXd Ad = Eigen::MatrixXd::Zero(K, K + 1), temp1 = Eigen::MatrixXd::Zero(K + 1, 4 * K + 2), temp2 = Eigen::MatrixXd::Zero(K, 4 * K + 2), temp3 = Eigen::MatrixXd::Zero(K, 4 * K + 2);
    temp1.block(0, K, K + 1, K + 1) = Eigen::MatrixXd::Identity(K + 1, K + 1);
    temp2.block(0, 0, K, K) = Eigen::MatrixXd::Identity(K, K);
    temp3.block(0, 3 * K + 2, K, K) = Eigen::MatrixXd::Identity(K, K);
    for(int i = 0; i < K; ++i){
        Ad(i, i) = -1;
        Ad(i, i + 1) = 1;
    }
    G_ = Eigen::MatrixXd::Zero(K + 2, 4 * K + 2);
    h_ = Eigen::MatrixXd::Zero(K + 2, 1);
    lambda_ = Eigen::MatrixXd::Zero(K + 2, 1);
    G_.block(0, 0, K, 4 * K + 2) = Ad * temp1 - 2 * (Ad * s_).asDiagonal() * temp2;
    G_(K, K) = 1;
    G_(K + 1, 2 * K) = 1;
    f_ = 2 * temp3.transpose() * Ad * s_;
    // std::cout << "===================" << std::endl;
    // Conic constraint
    Eigen::MatrixXd tempA, tempb;
    A1_.clear();
    b1_.clear();
    A2_.clear();
    b2_.clear();
    A3_.clear();
    b3_.clear();
    A41_.clear();
    b41_.clear();
    A42_.clear();
    b42_.clear();
    A43_.clear();
    b43_.clear();
    A44_.clear();
    b44_.clear();
    A51_.clear();
    b51_.clear();
    A52_.clear();
    b52_.clear();
    A53_.clear();
    b53_.clear();
    A54_.clear();
    b54_.clear();
    mu1_.clear();
    mu2_.clear();
    mu3_.clear();
    mu41_.clear();
    mu42_.clear();
    mu43_.clear();
    mu44_.clear();
    mu51_.clear();
    mu52_.clear();
    mu53_.clear();
    mu54_.clear();
    for(int k = 0; k <= K; ++k){
        if(k <= K - 1){
            tempA = Eigen::MatrixXd::Zero(3, 4 * K + 2);
            tempb = Eigen::MatrixXd::Zero(3, 1);
            tempA(0, 2 * K + 1 + k) = 1;
            tempA(0, 2 * K + 1 + k + 1) = 1;
            tempA(0, 3 * K + 2 + k) = 1;
            tempA(2, 2 * K + 1 + k) = 1;
            tempA(2, 2 * K + 1 + k + 1) = 1;
            tempA(2, 3 * K + 2 + k) = 1;
            tempb(1, 0) = 2;
            A1_.push_back(tempA);
            b1_.push_back(tempb);
            mu1_.push_back(Eigen::MatrixXd::Zero(3, 1));
        }
        // std::cout << "1===================" << std::endl;
        tempA = Eigen::MatrixXd::Zero(3, 4 * K + 2);
        tempb = Eigen::MatrixXd::Zero(3, 1);
        tempA(0, K + k) = 1;
        tempA(1, 2 * K + k + 1) = 1;
        tempA(2, K + k) = 1;
        tempb(0, 0) = 1;
        tempb(2, 0) = -1;
        A2_.push_back(tempA);
        b2_.push_back(tempb);
        mu2_.push_back(Eigen::MatrixXd::Zero(3, 1));
        // std::cout << "2===================" << std::endl;
        // tempA = Eigen::MatrixXd::Zero(2, 4 * K + 2);
        // tempb = Eigen::MatrixXd::Zero(2, 1);
        // tempA(1, K + k) = -1;
        // A3_.push_back(tempA);
        // b3_.push_back(tempb);
        // mu3_.push_back(Eigen::MatrixXd::Zero(2, 1));
        tempA = Eigen::MatrixXd::Zero(1, 4 * K + 2);
        tempb = Eigen::MatrixXd::Zero(1, 1);
        tempA(0, K + k) = 1;
        A3_.push_back(tempA);
        b3_.push_back(tempb);
        mu3_.push_back(Eigen::MatrixXd::Zero(1, 1));
        // std::cout << "3===================" << std::endl;
        // tempb = Eigen::MatrixXd::Zero(2, 1);
        // tempb(0, 0) = vmax_ * vmax_;
        // b41_.push_back(tempb);
        // b42_.push_back(tempb);
        // b43_.push_back(tempb);
        // b44_.push_back(tempb);
        // mu41_.push_back(Eigen::MatrixXd::Zero(2, 1));
        // mu42_.push_back(Eigen::MatrixXd::Zero(2, 1));
        // mu43_.push_back(Eigen::MatrixXd::Zero(2, 1));
        // mu44_.push_back(Eigen::MatrixXd::Zero(2, 1));
        // tempA = Eigen::MatrixXd::Zero(2, 4 * K + 2);
        // tempA(1, K + k) = px_dot_(k) * px_dot_(k);
        // A41_.push_back(tempA);
        // tempA(1, K + k) = py_dot_(k) * py_dot_(k);
        // A42_.push_back(tempA);
        // tempA(1, K + k) = -px_dot_(k) * px_dot_(k);
        // A43_.push_back(tempA);
        // tempA(1, K + k) = -py_dot_(k) * py_dot_(k);
        // A44_.push_back(tempA);

        tempb = Eigen::MatrixXd::Zero(1, 1);
        tempb(0, 0) = vmax_ * vmax_;
        b41_.push_back(tempb);
        b42_.push_back(tempb);
        b43_.push_back(tempb);
        b44_.push_back(tempb);
        mu41_.push_back(Eigen::MatrixXd::Zero(1, 1));
        mu42_.push_back(Eigen::MatrixXd::Zero(1, 1));
        mu43_.push_back(Eigen::MatrixXd::Zero(1, 1));
        mu44_.push_back(Eigen::MatrixXd::Zero(1, 1));
        tempA = Eigen::MatrixXd::Zero(1, 4 * K + 2);
        tempA(0, K + k) = -px_dot_(k) * px_dot_(k);
        A41_.push_back(tempA);
        tempA(0, K + k) = -py_dot_(k) * py_dot_(k);
        A42_.push_back(tempA);
        tempA(0, K + k) = px_dot_(k) * px_dot_(k);
        A43_.push_back(tempA);
        tempA(0, K + k) = py_dot_(k) * py_dot_(k);
        A44_.push_back(tempA);
        // std::cout << "5===================" << std::endl;
        if(k <= K - 1){
            // tempb = Eigen::MatrixXd::Zero(2, 1);
            // tempb(0, 0) = amax_;
            // b51_.push_back(tempb);
            // b52_.push_back(tempb);
            // b53_.push_back(tempb);
            // b54_.push_back(tempb);
            // mu51_.push_back(Eigen::MatrixXd::Zero(2, 1));
            // mu52_.push_back(Eigen::MatrixXd::Zero(2, 1));
            // mu53_.push_back(Eigen::MatrixXd::Zero(2, 1));
            // mu54_.push_back(Eigen::MatrixXd::Zero(2, 1));
            // tempA = Eigen::MatrixXd::Zero(2, 4 * K + 2);
            // tempA(1, k) = px_dot_(k);
            // tempA(1, K + k) = px_dot_dot_(k);
            // A51_.push_back(tempA);
            // tempA(1, k) = py_dot_(k);
            // tempA(1, K + k) = py_dot_dot_(k);
            // A52_.push_back(tempA);
            // tempA(1, k) = -px_dot_(k);
            // tempA(1, K + k) = -px_dot_dot_(k);
            // A53_.push_back(tempA);
            // tempA(1, k) = -py_dot_(k);
            // tempA(1, K + k) = -py_dot_dot_(k);
            // A54_.push_back(tempA);

            tempb = Eigen::MatrixXd::Zero(1, 1);
            tempb(0, 0) = amax_;
            b51_.push_back(tempb);
            b52_.push_back(tempb);
            b53_.push_back(tempb);
            b54_.push_back(tempb);
            mu51_.push_back(Eigen::MatrixXd::Zero(1, 1));
            mu52_.push_back(Eigen::MatrixXd::Zero(1, 1));
            mu53_.push_back(Eigen::MatrixXd::Zero(1, 1));
            mu54_.push_back(Eigen::MatrixXd::Zero(1, 1));
            tempA = Eigen::MatrixXd::Zero(1, 4 * K + 2);
            tempA(0, k) = -px_dot_(k);
            tempA(0, K + k) = -px_dot_dot_(k);
            A51_.push_back(tempA);
            tempA(0, k) = -py_dot_(k);
            tempA(0, K + k) = -py_dot_dot_(k);
            A52_.push_back(tempA);
            tempA(0, k) = px_dot_(k);
            tempA(0, K + k) = px_dot_dot_(k);
            A53_.push_back(tempA);
            tempA(0, k) = py_dot_(k);
            tempA(0, K + k) = py_dot_dot_(k);
            A54_.push_back(tempA);
        }
    }

    // std::cout << "6===================" << std::endl;
    // std::cout << (Ad * s_) << std::endl;
    // std::cout << temp1.block(0, K, K + 1, K + 1) << std::endl;

    int times = 0;
    double finalCost, e_cons = 1e-4, e_prec = 1e-4, res_cons = 1e5, res_prec = 1e5;
    Eigen::VectorXd x_opt(4 * K + 2);
    // int sign_x = px_[K] - px_[0] > 0 ? 1 : -1;
    // int sign_y = py_[K] - py_[0] > 0 ? 1 : -1;
    // for(int i = 0; i < K; ++i){
    //     if(i < K / 2)
    //         x_opt(i) = sign_x * amax_ / px_dot_(i);
    //     else
    //         x_opt(i) = -sign_x * amax_ / px_dot_(i);
    // }
    
    /* Set the minimization parameters */
    lbfgs::lbfgs_parameter_t params;
    params.mem_size = 32;
    params.past = 0;
    params.g_epsilon = 1e-8;
    params.min_step = 1e-32;
    params.max_linesearch = 128;
    params.delta = 1e-4;

    /* Start minimization */
    while((res_cons > e_cons || res_prec > e_prec) && times < 20){
        ++times;
        /* Start minimization */
        int ret = lbfgs::lbfgs_optimize(x_opt,
                                finalCost,
                                costFunction_TOPP,
                                monitorProgress,
                                this,
                                params);

        res_cons = (G_ * x_opt - h_).lpNorm<Eigen::Infinity>();
        Eigen::VectorXd g = f_ + G_.transpose() * (lambda_ + rho_ * (G_ * x_opt - h_));
        lambda_ = lambda_ + rho_ * (G_ * x_opt - h_);
        double temp;
        for(int k = 0; k <= K; ++k){
            // std::cout << "1===================" << std::endl;
            if(k <= K - 1){
                // cost = cost + (proj2SOC(mu1_[k] / rho_ - A1_[k] * x_opt - b1_[k], 3)).norm();
                // cost = cost + (proj2SOC(mu51_[k] / rho_ - A51_[k] * x_opt - b51_[k], 2)).norm();
                // cost = cost + (proj2SOC(mu52_[k] / rho_ - A52_[k] * x_opt - b52_[k], 2)).norm();
                // cost = cost + (proj2SOC(mu53_[k] / rho_ - A53_[k] * x_opt - b53_[k], 2)).norm();
                // cost = cost + (proj2SOC(mu54_[k] / rho_ - A54_[k] * x_opt - b54_[k], 2)).norm();
                g -= A1_[k].transpose() * proj2SOC(mu1_[k] - rho_ * (A1_[k] * x_opt + b1_[k]));
                g -= A51_[k].transpose() * proj2SOC(mu51_[k] - rho_ * (A51_[k] * x_opt + b51_[k]));
                g -= A51_[k].transpose() * proj2SOC(mu52_[k] - rho_ * (A52_[k] * x_opt + b52_[k]));
                g -= A51_[k].transpose() * proj2SOC(mu53_[k] - rho_ * (A53_[k] * x_opt + b53_[k]));
                g -= A51_[k].transpose() * proj2SOC(mu54_[k] - rho_ * (A54_[k] * x_opt + b54_[k]));

                
                temp = std::abs((mu1_[k] / rho_ - proj2SOC(mu1_[k] / rho_ - (A1_[k] * x_opt + b1_[k]))).lpNorm<Eigen::Infinity>());
                if(temp > res_cons)
                    res_cons = temp;
                temp = std::abs((mu51_[k] / rho_ - proj2SOC(mu51_[k] / rho_ - (A51_[k] * x_opt + b51_[k]))).lpNorm<Eigen::Infinity>());
                if(temp > res_cons)
                    res_cons = temp;
                temp = std::abs((mu52_[k] / rho_ - proj2SOC(mu52_[k] / rho_ - (A52_[k] * x_opt + b52_[k]))).lpNorm<Eigen::Infinity>());
                if(temp > res_cons)
                    res_cons = temp;
                temp = std::abs((mu53_[k] / rho_ - proj2SOC(mu53_[k] / rho_ - (A53_[k] * x_opt + b53_[k]))).lpNorm<Eigen::Infinity>());
                if(temp > res_cons)
                    res_cons = temp;
                temp = std::abs((mu54_[k] / rho_ - proj2SOC(mu54_[k] / rho_ - (A54_[k] * x_opt + b54_[k]))).lpNorm<Eigen::Infinity>());
                if(temp > res_cons)
                    res_cons = temp;

                mu1_[k] = proj2SOC(mu1_[k] - rho_ * (A1_[k] * x_opt + b1_[k]));
                mu51_[k] = proj2SOC(mu51_[k] - rho_ * (A51_[k] * x_opt + b51_[k]));
                mu52_[k] = proj2SOC(mu52_[k] - rho_ * (A52_[k] * x_opt + b52_[k]));
                mu53_[k] = proj2SOC(mu53_[k] - rho_ * (A53_[k] * x_opt + b53_[k]));
                mu54_[k] = proj2SOC(mu54_[k] - rho_ * (A54_[k] * x_opt + b54_[k]));


            }
            // std::cout << "2===================" << std::endl;
            // cost = cost + (proj2SOC(mu2_[k] / rho_ - A2_[k] * x_opt - b2_[k])).norm();
            // cost = cost + (proj2SOC(mu3_[k] / rho_ - A3_[k] * x_opt - b3_[k])).norm();
            // cost = cost + (proj2SOC(mu41_[k] / rho_ - A41_[k] * x_opt - b41_[k])).norm();
            // cost = cost + (proj2SOC(mu42_[k] / rho_ - A42_[k] * x_opt - b42_[k])).norm();
            // cost = cost + (proj2SOC(mu43_[k] / rho_ - A43_[k] * x_opt - b43_[k])).norm();
            // cost = cost + (proj2SOC(mu44_[k] / rho_ - A44_[k] * x_opt - b44_[k])).norm();
            // std::cout << "3===================" << std::endl;
            g -= A2_[k].transpose() * proj2SOC(mu2_[k] - rho_ * (A2_[k] * x_opt + b2_[k]));
            g -= A3_[k].transpose() * proj2SOC(mu3_[k] - rho_ * (A3_[k] * x_opt + b3_[k]));
            g -= A41_[k].transpose() * proj2SOC(mu41_[k] - rho_ * (A41_[k] * x_opt + b41_[k]));
            g -= A42_[k].transpose() * proj2SOC(mu42_[k] - rho_ * (A42_[k] * x_opt + b42_[k]));
            g -= A43_[k].transpose() * proj2SOC(mu43_[k] - rho_ * (A43_[k] * x_opt + b43_[k]));
            g -= A44_[k].transpose() * proj2SOC(mu44_[k] - rho_ * (A44_[k] * x_opt + b44_[k]));

            temp = std::abs((mu2_[k] / rho_ - proj2SOC(mu2_[k] / rho_ - (A2_[k] * x_opt + b2_[k]))).lpNorm<Eigen::Infinity>());
            if(temp > res_cons)
                res_cons = temp;
            temp = std::abs((mu3_[k] / rho_ - proj2SOC(mu3_[k] / rho_ - (A3_[k] * x_opt + b3_[k]))).lpNorm<Eigen::Infinity>());
            if(temp > res_cons)
                res_cons = temp;
            temp = std::abs((mu41_[k] / rho_ - proj2SOC(mu41_[k] / rho_ - (A41_[k] * x_opt + b41_[k]))).lpNorm<Eigen::Infinity>());
            if(temp > res_cons)
                res_cons = temp;
            temp = std::abs((mu42_[k] / rho_ - proj2SOC(mu42_[k] / rho_ - (A42_[k] * x_opt + b42_[k]))).lpNorm<Eigen::Infinity>());
            if(temp > res_cons)
                res_cons = temp;
            temp = std::abs((mu43_[k] / rho_ - proj2SOC(mu43_[k] / rho_ - (A43_[k] * x_opt + b43_[k]))).lpNorm<Eigen::Infinity>());
            if(temp > res_cons)
                res_cons = temp;
            temp = std::abs((mu44_[k] / rho_ - proj2SOC(mu44_[k] / rho_ - (A44_[k] * x_opt + b44_[k]))).lpNorm<Eigen::Infinity>());
            if(temp > res_cons)
                res_cons = temp;

            mu2_[k] = proj2SOC(mu2_[k] - rho_ * (A2_[k] * x_opt + b2_[k]));
            mu3_[k] = proj2SOC(mu3_[k] - rho_ * (A3_[k] * x_opt + b3_[k]));
            mu41_[k] = proj2SOC(mu41_[k] - rho_ * (A41_[k] * x_opt + b41_[k]));
            mu42_[k] = proj2SOC(mu42_[k] - rho_ * (A42_[k] * x_opt + b42_[k]));
            mu43_[k] = proj2SOC(mu43_[k] - rho_ * (A43_[k] * x_opt + b43_[k]));
            mu44_[k] = proj2SOC(mu44_[k] - rho_ * (A44_[k] * x_opt + b44_[k]));
        }
        rho_ = std::min(beta_, (1 + gamma_) * rho_);
        // Eigen::VectorXd v = mu_ / rho_ - A_ba_ * x - b_ba_;
        res_prec = g.lpNorm<Eigen::Infinity>();
        
        // Eigen::VectorXd proj = mu_ / rho_ - proj2SOC(v);
        

        // // params.g_epsilon = std::max(params.g_epsilon * 0.1, 1e-8);

        // /* Update */
        // mu_ = proj2SOC(mu_ - rho_ * (A_ba_ * x + b_ba_));
        // rho_ = std::min(beta_, (1 + gamma_) * rho_);

        std::cout << "================================" << std::endl
                << "iter times: " << times << std::endl
                << "res_cons: " << res_cons << std::endl
                << "res_prec: " << res_prec << std::endl
                << "L-BFGS Optimization Returned: " << ret << std::endl
                << "Minimized Cost: " << finalCost << std::endl;
                // << "\033[32m" << "      Minimized Cost: " << f_.dot(x) << "\033[0m" << std::endl
                // << "given Minimized Cost: " << -156.9589 << std::endl
                // << "\033[32m" << "      Optimal Variables: " << x.transpose() << "\033[0m" << std::endl
                // << "given Optimal Variables: " << x_opt.transpose() << std::endl << std::endl;
    }
    // int ret = lbfgs::lbfgs_optimize(x_opt,
    //                                 finalCost,
    //                                 costFunction_TOPP,
    //                                 monitorProgress,
    //                                 this,
    //                                 params);

    /* Report the result. */
    // std::cout << std::setprecision(4)
    //             << "================================" << std::endl
    //             << "L-BFGS Optimization Returned: " << ret << std::endl
    //             << "Minimized Cost: " << finalCost << std::endl
    //             << "Optimal Variables: " << std::endl;
    //             << x_opt.transpose() << std::endl;

    std::vector<double> v_x(K + 1), v_y(K + 1), acc_x(K), acc_y(K), vmax(K + 1), _vmax(K + 1), amax(K), _amax(K), x_0_K(K + 1), x_0_K_1(K), b_k(K + 1), s_K(K + 1), s_K_1(K);
    for(int k = 0; k <= K; ++k){
        // v_x[k] = px_dot_(k) * px_dot_(k) * x_opt(K + k);
        // v_y[k] = py_dot_(k) * py_dot_(k) * x_opt(K + k);
        v_x[k] = px_dot_(k) * sqrt(x_opt(K + k));
        v_y[k] = py_dot_(k) * sqrt(x_opt(K + k));
        x_0_K[k] = k;
        s_K[k] = s_(k);
        b_k[k] = x_opt(K + k);
        vmax[k] = vmax_;
        _vmax[k] = -vmax_;
        if(k <= K - 1){
            acc_x[k] = px_dot_dot_(k) * x_opt(K + k) + px_dot_(k) * x_opt(k);
            acc_y[k] = py_dot_dot_(k) * x_opt(K + k) + py_dot_(k) * x_opt(k);
            x_0_K_1[k] = k;
            s_K_1[k] = s_(k);
            amax[k] = amax_;
            _amax[k] = -amax_;
        }
    }
    plt::figure();
    plt::subplot(3,1,1);
    plt::plot(s_K, v_x, {{"label", "v_x"}});
    plt::plot(s_K, v_y, {{"label", "v_y"}});
    plt::plot(s_K, vmax, {{"label", "v_max"}, {"ls", "--"}});
    plt::plot(s_K, _vmax, {{"label", "-v_max"}, {"ls", "--"}});
    plt::legend();
    plt::grid(true);

    plt::subplot(3,1,2);
    plt::plot(s_K_1, acc_x, {{"label", "a_x"}});
    plt::plot(s_K_1, acc_y, {{"label", "a_y"}});
    plt::plot(s_K_1, amax, {{"label", "a_max"}, {"ls", "--"}});
    plt::plot(s_K_1, _amax, {{"label", "-a_max"}, {"ls", "--"}});
    plt::legend();
    plt::grid(true);

    plt::subplot(3,1,3);
    plt::plot(s_K, b_k, {{"label", "b_k"}});
    plt::legend();
    plt::grid(true);
    plt::show();
}