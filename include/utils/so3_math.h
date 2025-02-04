#ifndef SO3_MATH_H
#define SO3_MATH_H

#include <Eigen/Core>
#include <math.h>

// 定义一个宏，用于生成反对称矩阵（斜对称矩阵）
#define SKEW_SYM_MATRX(v) 0.0, -v[2], v[1], v[2], 0.0, -v[0], -v[1], v[0], 0.0

// 计算旋转向量的指数映射，返回旋转矩阵
template <typename T> 
Eigen::Matrix<T, 3, 3> Exp(const Eigen::Matrix<T, 3, 1> &&ang)
{
  // 计算旋转向量的范数（大小）
  T ang_norm = ang.norm();
  // 创建单位矩阵
  Eigen::Matrix<T, 3, 3> Eye3 = Eigen::Matrix<T, 3, 3>::Identity();
  
  // 如果旋转向量的范数大于一个很小的值
  if (ang_norm > 0.0000001)
  {
    // 计算旋转轴
    Eigen::Matrix<T, 3, 1> r_axis = ang / ang_norm;
    Eigen::Matrix<T, 3, 3> K;
    K << SKEW_SYM_MATRX(r_axis); // 生成反对称矩阵

    // 罗德里格斯变换，计算旋转矩阵
    return Eye3 + std::sin(ang_norm) * K + (1.0 - std::cos(ang_norm)) * K * K;
  }
  else { return Eye3; } // 返回单位矩阵
}

// 计算角速度的指数映射，返回旋转矩阵
template <typename T, typename Ts> 
Eigen::Matrix<T, 3, 3> Exp(const Eigen::Matrix<T, 3, 1> &ang_vel, const Ts &dt)
{
  // 计算角速度的范数
  T ang_vel_norm = ang_vel.norm();
  Eigen::Matrix<T, 3, 3> Eye3 = Eigen::Matrix<T, 3, 3>::Identity();

  // 如果角速度的范数大于一个很小的值
  if (ang_vel_norm > 0.0000001)
  {
    // 计算旋转轴
    Eigen::Matrix<T, 3, 1> r_axis = ang_vel / ang_vel_norm;
    Eigen::Matrix<T, 3, 3> K;
    K << SKEW_SYM_MATRX(r_axis); // 生成反对称矩阵

    T r_ang = ang_vel_norm * dt; // 计算旋转角度

    // 罗德里格斯变换，计算旋转矩阵
    return Eye3 + std::sin(r_ang) * K + (1.0 - std::cos(r_ang)) * K * K;
  }
  else { return Eye3; } // 返回单位矩阵
}

// 计算给定三个分量的旋转向量的指数映射，返回旋转矩阵
template <typename T> 
Eigen::Matrix<T, 3, 3> Exp(const T &v1, const T &v2, const T &v3)
{
  // 计算旋转向量的范数
  T &&norm = sqrt(v1 * v1 + v2 * v2 + v3 * v3);
  Eigen::Matrix<T, 3, 3> Eye3 = Eigen::Matrix<T, 3, 3>::Identity();
  
  // 如果范数大于一个很小的值
  if (norm > 0.00001)
  {
    // 计算旋转轴
    T r_ang[3] = {v1 / norm, v2 / norm, v3 / norm};
    Eigen::Matrix<T, 3, 3> K;
    K << SKEW_SYM_MATRX(r_ang); // 生成反对称矩阵

    // 罗德里格斯变换，计算旋转矩阵
    return Eye3 + std::sin(norm) * K + (1.0 - std::cos(norm)) * K * K;
  }
  else { return Eye3; } // 返回单位矩阵
}

// 计算旋转矩阵的对数映射，返回旋转向量
template <typename T> 
Eigen::Matrix<T, 3, 1> Log(const Eigen::Matrix<T, 3, 3> &R)
{
  // 计算旋转角度
  T theta = (R.trace() > 3.0 - 1e-6) ? 0.0 : std::acos(0.5 * (R.trace() - 1));
  // 生成旋转向量
  Eigen::Matrix<T, 3, 1> K(R(2, 1) - R(1, 2), R(0, 2) - R(2, 0), R(1, 0) - R(0, 1));
  
  // 返回旋转向量
  return (std::abs(theta) < 0.001) ? (0.5 * K) : (0.5 * theta / std::sin(theta) * K);
}

// 将旋转矩阵转换为欧拉角，返回欧拉角向量
template <typename T> 
Eigen::Matrix<T, 3, 1> RotMtoEuler(const Eigen::Matrix<T, 3, 3> &rot)
{
  // 计算sy的值
  T sy = sqrt(rot(0, 0) * rot(0, 0) + rot(1, 0) * rot(1, 0));
  bool singular = sy < 1e-6; // 检测奇异性
  T x, y, z;
  
  // 如果不奇异
  if (!singular)
  {
    x = atan2(rot(2, 1), rot(2, 2)); // 计算x轴旋转
    y = atan2(-rot(2, 0), sy); // 计算y轴旋转
    z = atan2(rot(1, 0), rot(0, 0)); // 计算z轴旋转
  }
  else // 如果奇异
  {
    x = atan2(-rot(1, 2), rot(1, 1)); // 计算x轴旋转
    y = atan2(-rot(2, 0), sy); // 计算y轴旋转
    z = 0; // z轴旋转为0
  }
  
  // 返回欧拉角向量
  Eigen::Matrix<T, 3, 1> ang(x, y, z);
  return ang;
}

#endif // SO3_MATH_H
