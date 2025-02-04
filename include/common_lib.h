/* 
This file is part of FAST-LIVO2: Fast, Direct LiDAR-Inertial-Visual Odometry.

Developer: Chunran Zheng <zhengcr@connect.hku.hk>

For commercial use, please contact me at <zhengcr@connect.hku.hk> or
Prof. Fu Zhang at <fuzhang@hku.hk>.

This file is subject to the terms and conditions outlined in the 'LICENSE' file,
which is included as part of this source code package.
*/

#ifndef COMMON_LIB_H
#define COMMON_LIB_H

#include <utils/so3_math.h> // 导入SO(3)数学工具
#include <utils/types.h>    // 导入自定义类型
#include <utils/color.h>    // 导入颜色处理工具
#include <opencv2/opencv.hpp> // 导入OpenCV库
#include <sensor_msgs/Imu.h> // 导入IMU传感器消息
#include <sophus/se3.h>      // 导入Sophus库中的SE(3)变换
#include <tf/transform_broadcaster.h> // 导入TF变换广播器

using namespace std; // 使用标准命名空间
using namespace Eigen; // 使用Eigen库命名空间
using namespace Sophus; // 使用Sophus库命名空间

#define print_line std::cout << __FILE__ << ", " << __LINE__ << std::endl; // 打印当前文件和行号
#define G_m_s2 (9.81)   // 广州/中国的重力常数
#define DIM_STATE (19)  // 状态维度 (设定SO(3)维度为3)
#define INIT_COV (0.01) // 初始化协方差
#define SIZE_LARGE (500) // 大尺寸常量
#define SIZE_SMALL (100) // 小尺寸常量
#define VEC_FROM_ARRAY(v) v[0], v[1], v[2] // 从数组创建3维向量
#define MAT_FROM_ARRAY(v) v[0], v[1], v[2], v[3], v[4], v[5], v[6], v[7], v[8] // 从数组创建3x3矩阵
#define DEBUG_FILE_DIR(name) (string(string(ROOT_DIR) + "Log/" + name)) // 设置调试文件目录

// 激光雷达类型枚举
enum LID_TYPE
{
  AVIA = 1,
  VELO16 = 2,
  OUST64 = 3,
  L515 = 4,
  XT32 = 5,
  PANDAR128 = 6
};

// SLAM模式枚举
enum SLAM_MODE
{
  ONLY_LO = 0, // 仅使用定位
  ONLY_LIO = 1, // 仅使用激光与惯性测量
  LIVO = 2 // 同时使用激光、惯性与视觉
};

// EKF状态枚举
enum EKF_STATE
{
  WAIT = 0, // 等待状态
  VIO = 1, // 视觉惯性测量状态
  LIO = 2 // 激光惯性测量状态
};

// 测量组结构体
struct MeasureGroup
{
  double vio_time; // 视觉惯性测量时间
  double lio_time; // 激光惯性测量时间
  deque<sensor_msgs::Imu::ConstPtr> imu; // 存储IMU数据的队列
  cv::Mat img; // 存储图像数据

  MeasureGroup() // 构造函数
  {
    vio_time = 0.0; // 初始化视觉惯性测量时间
    lio_time = 0.0; // 初始化激光惯性测量时间
  };
};

// 激光测量组结构体
struct LidarMeasureGroup
{
  double lidar_frame_beg_time; // 激光帧开始时间
  double lidar_frame_end_time; // 激光帧结束时间
  double last_lio_update_time; // 上一次激光惯性更新的时间
  PointCloudXYZI::Ptr lidar; // 激光点云数据
  PointCloudXYZI::Ptr pcl_proc_cur; // 当前处理的点云
  PointCloudXYZI::Ptr pcl_proc_next; // 下一步处理的点云
  deque<struct MeasureGroup> measures; // 存储测量组的队列
  EKF_STATE lio_vio_flg; // 激光与视觉的状态标志
  int lidar_scan_index_now; // 当前激光扫描索引

  LidarMeasureGroup() // 构造函数
  {
    lidar_frame_beg_time = -0.0; // 初始化激光帧开始时间
    lidar_frame_end_time = 0.0; // 初始化激光帧结束时间
    last_lio_update_time = -1.0; // 初始化最后一次激光惯性更新的时间
    lio_vio_flg = WAIT; // 初始化状态标志为等待
    this->lidar.reset(new PointCloudXYZI()); // 初始化激光点云
    this->pcl_proc_cur.reset(new PointCloudXYZI()); // 初始化当前点云
    this->pcl_proc_next.reset(new PointCloudXYZI()); // 初始化下一步点云
    this->measures.clear(); // 清空测量组队列
    lidar_scan_index_now = 0; // 初始化激光扫描索引
    last_lio_update_time = -1.0; // 初始化最后一次更新的时间
  };
};

// 带有协方差的点结构体
typedef struct pointWithVar
{
  Eigen::Vector3d point_b;     // 激光坐标系下的点
  Eigen::Vector3d point_i;     // IMU坐标系下的点
  Eigen::Vector3d point_w;     // 世界坐标系下的点
  Eigen::Matrix3d var_nostate; // 去除状态协方差后的方差
  Eigen::Matrix3d body_var;    // 机体方差
  Eigen::Matrix3d var;          // 方差
  Eigen::Matrix3d point_crossmat; // 点的叉乘矩阵
  Eigen::Vector3d normal;      // 法向量

  pointWithVar() // 构造函数
  {
    var_nostate = Eigen::Matrix3d::Zero(); // 初始化去状态协方差的方差为零
    var = Eigen::Matrix3d::Zero(); // 初始化方差为零
    body_var = Eigen::Matrix3d::Zero(); // 初始化机体方差为零
    point_crossmat = Eigen::Matrix3d::Zero(); // 初始化点的叉乘矩阵为零
    point_b = Eigen::Vector3d::Zero(); // 初始化激光坐标系下的点为零
    point_i = Eigen::Vector3d::Zero(); // 初始化IMU坐标系下的点为零
    point_w = Eigen::Vector3d::Zero(); // 初始化世界坐标系下的点为零
    normal = Eigen::Vector3d::Zero(); // 初始化法向量为零
  };
} pointWithVar;

// 状态组结构体
struct StatesGroup
{
  StatesGroup() // 默认构造函数
  {
    this->rot_end = M3D::Identity(); // 初始化旋转矩阵为单位矩阵
    this->pos_end = V3D::Zero(); // 初始化位置为零
    this->vel_end = V3D::Zero(); // 初始化速度为零
    this->bias_g = V3D::Zero(); // 初始化陀螺仪偏差为零
    this->bias_a = V3D::Zero(); // 初始化加速度计偏差为零
    this->gravity = V3D::Zero(); // 初始化重力加速度为零
    this->inv_expo_time = 1.0; // 初始化逆曝光时间为1.0
    this->cov = MD(DIM_STATE, DIM_STATE)::Identity() * INIT_COV; // 初始化协方差矩阵
    this->cov(6, 6) = 0.00001; // 设置协方差矩阵的特定值
    this->cov.block<9, 9>(10, 10) = MD(9, 9)::Identity() * 0.00001; // 设置协方差矩阵的特定块
  };

  StatesGroup(const StatesGroup &b) // 拷贝构造函数
  {
    this->rot_end = b.rot_end; // 拷贝旋转矩阵
    this->pos_end = b.pos_end; // 拷贝位置
    this->vel_end = b.vel_end; // 拷贝速度
    this->bias_g = b.bias_g; // 拷贝陀螺仪偏差
    this->bias_a = b.bias_a; // 拷贝加速度计偏差
    this->gravity = b.gravity; // 拷贝重力加速度
    this->inv_expo_time = b.inv_expo_time; // 拷贝逆曝光时间
    this->cov = b.cov; // 拷贝协方差矩阵
  };

  StatesGroup &operator=(const StatesGroup &b) // 赋值操作符重载
  {
    this->rot_end = b.rot_end; // 赋值旋转矩阵
    this->pos_end = b.pos_end; // 赋值位置
    this->vel_end = b.vel_end; // 赋值速度
    this->bias_g = b.bias_g; // 赋值陀螺仪偏差
    this->bias_a = b.bias_a; // 赋值加速度计偏差
    this->gravity = b.gravity; // 赋值重力加速度
    this->inv_expo_time = b.inv_expo_time; // 赋值逆曝光时间
    this->cov = b.cov; // 赋值协方差矩阵
    return *this; // 返回当前对象
  };

  StatesGroup operator+(const Matrix<double, DIM_STATE, 1> &state_add) // 加法操作符重载
  {
    StatesGroup a; // 创建新的状态组
    a.rot_end = this->rot_end * Exp(state_add(0, 0), state_add(1, 0), state_add(2, 0)); // 更新旋转
    a.pos_end = this->pos_end + state_add.block<3, 1>(3, 0); // 更新位置
    a.inv_expo_time = this->inv_expo_time + state_add(6, 0); // 更新逆曝光时间
    a.vel_end = this->vel_end + state_add.block<3, 1>(7, 0); // 更新速度
    a.bias_g = this->bias_g + state_add.block<3, 1>(10, 0); // 更新陀螺仪偏差
    a.bias_a = this->bias_a + state_add.block<3, 1>(13, 0); // 更新加速度计偏差
    a.gravity = this->gravity + state_add.block<3, 1>(16, 0); // 更新重力加速度

    a.cov = this->cov; // 赋值协方差矩阵
    return a; // 返回新的状态组
  };

  StatesGroup &operator+=(const Matrix<double, DIM_STATE, 1> &state_add) // 加法赋值操作符重载
  {
    this->rot_end = this->rot_end * Exp(state_add(0, 0), state_add(1, 0), state_add(2, 0)); // 更新旋转
    this->pos_end += state_add.block<3, 1>(3, 0); // 更新位置
    this->inv_expo_time += state_add(6, 0); // 更新逆曝光时间
    this->vel_end += state_add.block<3, 1>(7, 0); // 更新速度
    this->bias_g += state_add.block<3, 1>(10, 0); // 更新陀螺仪偏差
    this->bias_a += state_add.block<3, 1>(13, 0); // 更新加速度计偏差
    this->gravity += state_add.block<3, 1>(16, 0); // 更新重力加速度
    return *this; // 返回当前对象
  };

  Matrix<double, DIM_STATE, 1> operator-(const StatesGroup &b) // 减法操作符重载
  {
    Matrix<double, DIM_STATE, 1> a; // 创建新的向量
    M3D rotd(b.rot_end.transpose() * this->rot_end); // 计算旋转差异
    a.block<3, 1>(0, 0) = Log(rotd); // 计算旋转差异的对数
    a.block<3, 1>(3, 0) = this->pos_end - b.pos_end; // 计算位置差异
    a(6, 0) = this->inv_expo_time - b.inv_expo_time; // 计算逆曝光时间差异
    a.block<3, 1>(7, 0) = this->vel_end - b.vel_end; // 计算速度差异
    a.block<3, 1>(10, 0) = this->bias_g - b.bias_g; // 计算陀螺仪偏差差异
    a.block<3, 1>(13, 0) = this->bias_a - b.bias_a; // 计算加速度计偏差差异
    a.block<3, 1>(16, 0) = this->gravity - b.gravity; // 计算重力加速度差异
    return a; // 返回差异向量
  };

  void resetpose() // 重置姿态函数
  {
    this->rot_end = M3D::Identity(); // 重置旋转矩阵为单位矩阵
    this->pos_end = V3D::Zero(); // 重置位置为零
    this->vel_end = V3D::Zero(); // 重置速度为零
  }

  M3D rot_end;                              // 激光点结束时的估计姿态（旋转矩阵）
  V3D pos_end;                              // 激光点结束时的估计位置（世界坐标系）
  V3D vel_end;                              // 激光点结束时的估计速度（世界坐标系）
  double inv_expo_time;                     // 估计的逆曝光时间（无尺度）
  V3D bias_g;                               // 陀螺仪偏差
  V3D bias_a;                               // 加速度计偏差
  V3D gravity;                              // 估计的重力加速度
  Matrix<double, DIM_STATE, DIM_STATE> cov; // 状态协方差矩阵
};

// 设置6D姿态的模板函数
template <typename T>
auto set_pose6d(const double t, const Matrix<T, 3, 1> &a, const Matrix<T, 3, 1> &g, const Matrix<T, 3, 1> &v, const Matrix<T, 3, 1> &p,
                const Matrix<T, 3, 3> &R)
{
  Pose6D rot_kp; // 创建6D姿态对象
  rot_kp.offset_time = t; // 设置时间偏移
  for (int i = 0; i < 3; i++)
  {
    rot_kp.acc[i] = a(i); // 设置加速度
    rot_kp.gyr[i] = g(i); // 设置陀螺仪数据
    rot_kp.vel[i] = v(i); // 设置速度
    rot_kp.pos[i] = p(i); // 设置位置
    for (int j = 0; j < 3; j++)
      rot_kp.rot[i * 3 + j] = R(i, j); // 设置旋转矩阵
  }
  // Map<M3D>(rot_kp.rot, 3,3) = R; // 将旋转矩阵映射到rot_kp中
  return move(rot_kp); // 返回设置好的6D姿态
}

#endif
