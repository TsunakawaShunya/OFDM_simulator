#include </Volumes/USB1/eigen-3.4.0/Eigen/Core>
#include </Volumes/USB1/eigen-3.4.0/Eigen/Eigen>
#include <iostream>

int main() {
    Eigen::Vector3d vec1(1, 2, 3);
    Eigen::Vector3d vec2(4, 5, 6);

    std::cout << "vec1" << std::endl;
    std::cout << vec1 << std::endl;
    std::cout << "vec2" << std::endl;
    std::cout << vec2 << std::endl;
    std::cout << "sum" << std::endl;
    std::cout << vec1 + vec2 << std::endl;
    
    return 0;
}