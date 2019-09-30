#include <iostream>
#include <fstream>
#include <string>
#include <sstream>
#include <cmath>
#include <vector>
#include <Eigen/Core>
#include <fftw3.h>
#include <typeinfo>

using namespace std;
using namespace Eigen;

typedef Eigen::Matrix<double,6,1> Vec6;
typedef Eigen::Vector3d          Vec3;
typedef Eigen::Vector2d          Vec2;

const double pi = M_PI;