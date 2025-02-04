#ifndef CAMERA_INFO_H
#define CAMERA_INFO_H

#include <algorithm>
#include <Eigen/Dense>
#include <vector>
#include <fstream>
#include <filesystem>

enum class CAMERA_MODEL {
    SIMPLE_PINHOLE = 0,
    PINHOLE = 1,
    SIMPLE_RADIAL = 2,
    RADIAL = 3,
    OPENCV = 4,
    OPENCV_FISHEYE = 5,
    FULL_OPENCV = 6,
    FOV = 7,
    SIMPLE_RADIAL_FISHEYE = 8,
    RADIAL_FISHEYE = 9,
    THIN_PRISM_FISHEYE = 10,
    UNDEFINED = 11
};

// This class stores all information about a camera at loading time
// To me this seems to be double work, since we already have a Camera class
// I guess this can be removed later on
// TODO: Check and remove this struct if possible
struct CameraInfo {
    uint32_t _camera_ID;
    Eigen::Matrix3f _R; // rotation  matrix
    Eigen::Vector3f _T; // translation vector
    float _fov_x;
    float _fov_y;
    std::string _image_name;
    std::string _image_path;
    CAMERA_MODEL _camera_model;
    int _width;
    int _height;
    int _img_w;
    int _img_h;
    int _channels;
    std::vector<double> _params;
    unsigned char* _img_data; // shallow copy is fine here. No ownership
};


#endif