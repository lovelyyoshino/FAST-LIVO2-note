#ifndef TYPES_H
#define TYPES_H

#include <Eigen/Eigen>
#include <pcl/point_cloud.h>
#include <pcl/point_types.h>

// 定义点类型，包含位置、强度和法线信息
typedef pcl::PointXYZINormal PointType;
// 定义带有RGB颜色信息的点类型
typedef pcl::PointXYZRGB PointTypeRGB;
// 定义带有RGBA颜色信息的点类型
typedef pcl::PointXYZRGBA PointTypeRGBA;
// 定义包含PointType的点云类型
typedef pcl::PointCloud<PointType> PointCloudXYZI;
// 定义使用Eigen对齐分配器的点向量类型
typedef std::vector<PointType, Eigen::aligned_allocator<PointType>> PointVector;
// 定义包含PointTypeRGB的点云类型
typedef pcl::PointCloud<PointTypeRGB> PointCloudXYZRGB;
// 定义包含PointTypeRGBA的点云类型
typedef pcl::PointCloud<PointTypeRGBA> PointCloudXYZRGBA;

// 定义二维和三维向量及矩阵类型
typedef Eigen::Vector2f V2F; // 单精度二维向量
typedef Eigen::Vector2d V2D; // 双精度二维向量
typedef Eigen::Vector3d V3D; // 双精度三维向量
typedef Eigen::Matrix3d M3D; // 双精度3x3矩阵
typedef Eigen::Vector3f V3F; // 单精度三维向量
typedef Eigen::Matrix3f M3F; // 单精度3x3矩阵

// 定义宏以便于创建特定大小的双精度和单精度矩阵
#define MD(a, b) Eigen::Matrix<double, (a), (b)> // 双精度矩阵
#define VD(a) Eigen::Matrix<double, (a), 1> // 双精度列向量
#define MF(a, b) Eigen::Matrix<float, (a), (b)> // 单精度矩阵
#define VF(a) Eigen::Matrix<float, (a), 1> // 单精度列向量

// 定义一个结构体，用于表示6D姿态
struct Pose6D
{
  /*** 在IMU测量时刻的预积分Lidar状态 ***/
  double offset_time; // IMU测量相对于第一个Lidar点的偏移时间
  double acc[3];      // 在Lidar原点的预积分总加速度（全局坐标系）
  double gyr[3];      // 在Lidar原点的无偏角速度（机体坐标系）
  double vel[3];      // 在Lidar原点的预积分速度（全局坐标系）
  double pos[3];      // 在Lidar原点的预积分位置（全局坐标系）
  double rot[9];      // 在Lidar原点的预积分旋转（全局坐标系），以列优先的方式存储的3x3旋转矩阵
};

#endif // TYPES_H
