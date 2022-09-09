#include "path_generator/path_generator.h"

int main(int argc, char **argv){
    ros::init(argc, argv, "path_generator_node");
    ros::NodeHandle nh("~");
    PathGenerator path_generator;
    path_generator.init(nh);
    ros::spin();
    return 0;
}