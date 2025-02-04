/* 
This file is part of FAST-LIVO2: Fast, Direct LiDAR-Inertial-Visual Odometry.

Developer: Chunran Zheng <zhengcr@connect.hku.hk>

For commercial use, please contact me at <zhengcr@connect.hku.hk> or
Prof. Fu Zhang at <fuzhang@hku.hk>.

This file is subject to the terms and conditions outlined in the 'LICENSE' file,
which is included as part of this source code package.
*/
#include "preprocess.h"

// 定义返回值常量
#define RETURN0 0x00
#define RETURN0AND1 0x10

// 构造函数
Preprocess::Preprocess() : feature_enabled(0), lidar_type(AVIA), blind(0.01), point_filter_num(1)
{
  inf_bound = 10; // 设定无穷边界
  N_SCANS = 6; // 设定扫描线数量
  group_size = 8; // 设定分组大小
  disA = 0.01; // 设置距离A
  disA = 0.1; // B?
  p2l_ratio = 225; // 点到线的比率
  limit_maxmid = 6.25; // 最大中间限制
  limit_midmin = 6.25; // 中间最小限制
  limit_maxmin = 3.24; // 最大最小限制
  jump_up_limit = 170.0; // 跳跃上限
  jump_down_limit = 8.0; // 跳跃下限
  cos160 = 160.0; // 160度的余弦值
  edgea = 2; // 边缘A
  edgeb = 0.1; // 边缘B
  smallp_intersect = 172.5; // 小点交集
  smallp_ratio = 1.2; // 小点比率
  given_offset_time = false; // 是否给定偏移时间

  // 将角度转换为余弦值
  jump_up_limit = cos(jump_up_limit / 180 * M_PI);
  jump_down_limit = cos(jump_down_limit / 180 * M_PI);
  cos160 = cos(cos160 / 180 * M_PI);
  smallp_intersect = cos(smallp_intersect / 180 * M_PI);
}

// 析构函数
Preprocess::~Preprocess() {}

// 设置参数
void Preprocess::set(bool feat_en, int lid_type, double bld, int pfilt_num)
{
  feature_enabled = feat_en; // 设置特征启用状态
  lidar_type = lid_type; // 设置激光雷达类型
  blind = bld; // 设置盲区
  point_filter_num = pfilt_num; // 设置点过滤数量
}

// 处理来自livox_ros_driver的自定义消息
void Preprocess::process(const livox_ros_driver::CustomMsg::ConstPtr &msg, PointCloudXYZI::Ptr &pcl_out)
{
  avia_handler(msg); // 处理特定于AVIA的逻辑
  *pcl_out = pl_surf; // 将处理后的点云输出
}

// 处理来自sensor_msgs的点云消息
void Preprocess::process(const sensor_msgs::PointCloud2::ConstPtr &msg, PointCloudXYZI::Ptr &pcl_out)
{
  switch (lidar_type) // 根据激光雷达类型选择处理方法
  {
  case OUST64:
    oust64_handler(msg); // 处理OUST64类型
    break;

  case VELO16:
    velodyne_handler(msg); // 处理VELO16类型
    break;

  case L515:
    l515_handler(msg); // 处理L515类型
    break;

  case XT32:
    xt32_handler(msg); // 处理XT32类型
    break;

  case PANDAR128:
    Pandar128_handler(msg); // 处理Pandar128类型
    break;
  default:
    printf("Error LiDAR Type: %d \n", lidar_type); // 输出错误信息
    break;
  }
  *pcl_out = pl_surf; // 将处理后的点云输出
}

// 处理AVIA类型的点云数据
void Preprocess::avia_handler(const livox_ros_driver::CustomMsg::ConstPtr &msg)
{
  pl_surf.clear(); // 清空表面点云
  pl_corn.clear(); // 清空边缘点云
  pl_full.clear(); // 清空全量点云
  double t1 = omp_get_wtime(); // 记录开始时间
  int plsize = msg->point_num; // 获取输入点数
  printf("[ Preprocess ] Input point number: %d \n", plsize); // 输出输入点数

  pl_corn.reserve(plsize); // 预留空间
  pl_surf.reserve(plsize); // 预留空间
  pl_full.resize(plsize); // 调整全量点云大小

  // 初始化扫描缓冲区
  for (int i = 0; i < N_SCANS; i++)
  {
    pl_buff[i].clear();
    pl_buff[i].reserve(plsize); // 预留空间
  }
  uint valid_num = 0; // 有效点计数

  if (feature_enabled) // 如果启用特征提取
  {
    for (uint i = 1; i < plsize; i++)
    {
      if ((msg->points[i].line < N_SCANS) && ((msg->points[i].tag & 0x30) == 0x10)) // 检查点的线号和标签
      {
        pl_full[i].x = msg->points[i].x; // 复制X坐标
        pl_full[i].y = msg->points[i].y; // 复制Y坐标
        pl_full[i].z = msg->points[i].z; // 复制Z坐标
        pl_full[i].intensity = msg->points[i].reflectivity; // 复制反射强度
        pl_full[i].curvature = msg->points[i].offset_time / float(1000000); // 将偏移时间作为曲率

        bool is_new = false; // 新点标志
        // 检查当前点与前一个点的坐标差异
        if ((abs(pl_full[i].x - pl_full[i - 1].x) > 1e-7) || (abs(pl_full[i].y - pl_full[i - 1].y) > 1e-7) ||
            (abs(pl_full[i].z - pl_full[i - 1].z) > 1e-7))
        {
          pl_buff[msg->points[i].line].push_back(pl_full[i]); // 将当前点添加到相应的扫描缓冲区
        }
      }
    }
    static int count = 0; // 计数器
    static double time = 0.0; // 累计时间
    count++; // 增加计数
    double t0 = omp_get_wtime(); // 记录处理开始时间
    for (int j = 0; j < N_SCANS; j++) // 遍历每条扫描线
    {
      if (pl_buff[j].size() <= 5) continue; // 如果点数不足则跳过
      pcl::PointCloud<PointType> &pl = pl_buff[j]; // 获取当前扫描线的点云
      plsize = pl.size(); // 获取点云大小
      vector<orgtype> &types = typess[j]; // 获取当前扫描线的类型信息
      types.clear(); // 清空类型信息
      types.resize(plsize); // 调整类型信息大小
      plsize--; // 调整点云大小（最后一个点单独处理）
      for (uint i = 0; i < plsize; i++)
      {
        types[i].range = pl[i].x * pl[i].x + pl[i].y * pl[i].y; // 计算范围
        vx = pl[i].x - pl[i + 1].x; // 计算X差
        vy = pl[i].y - pl[i + 1].y; // 计算Y差
        vz = pl[i].z - pl[i + 1].z; // 计算Z差
        types[i].dista = vx * vx + vy * vy + vz * vz; // 计算距离
      }
      types[plsize].range = pl[plsize].x * pl[plsize].x + pl[plsize].y * pl[plsize].y; // 处理最后一个点的范围
      give_feature(pl, types); // 提取特征
      // pl_surf += pl; // 可以选择性添加到表面点云
    }
    time += omp_get_wtime() - t0; // 累计处理时间
    printf("Feature extraction time: %lf \n", time / count); // 输出特征提取时间
  }
  else // 如果不启用特征提取
  {
    for (uint i = 0; i < plsize; i++)
    {
      if ((msg->points[i].line < N_SCANS)) // 检查点的线号
      {
        valid_num++; // 有效点计数增加

        pl_full[i].x = msg->points[i].x; // 复制X坐标
        pl_full[i].y = msg->points[i].y; // 复制Y坐标
        pl_full[i].z = msg->points[i].z; // 复制Z坐标
        pl_full[i].intensity = msg->points[i].reflectivity; // 复制反射强度
        pl_full[i].curvature = msg->points[i].offset_time / float(1000000); // 将偏移时间作为曲率

        if (i == 0)
          pl_full[i].curvature = fabs(pl_full[i].curvature) < 1.0 ? pl_full[i].curvature : 0.0; // 处理第一个点的曲率
        else
        {
          // if(fabs(pl_full[i].curvature - pl_full[i - 1].curvature) > 1.0) ROS_ERROR("time jump: %f", fabs(pl_full[i].curvature - pl_full[i - 1].curvature));
          pl_full[i].curvature = fabs(pl_full[i].curvature - pl_full[i - 1].curvature) < 1.0
                                     ? pl_full[i].curvature
                                     : pl_full[i - 1].curvature + 0.004166667f; // float(100/24000)
        }

        if (valid_num % point_filter_num == 0) // 根据过滤数量选择性加入点云
        {
          if (pl_full[i].x * pl_full[i].x + pl_full[i].y * pl_full[i].y + pl_full[i].z * pl_full[i].z >= blind_sqr)
          {
            pl_surf.push_back(pl_full[i]); // 添加到表面点云
            // if (i % 100 == 0 || i == 0) printf("pl_full[i].curvature: %f \n",
            // pl_full[i].curvature);
          }
        }
      }
    }
  }
  printf("[ Preprocess ] Output point number: %zu \n", pl_surf.points.size()); // 输出处理后的点云数量
}

// 处理L515类型的点云数据
void Preprocess::l515_handler(const sensor_msgs::PointCloud2::ConstPtr &msg)
{
  pl_surf.clear(); // 清空表面点云
  pl_corn.clear(); // 清空边缘点云
  pl_full.clear(); // 清空全量点云
  pcl::PointCloud<pcl::PointXYZRGB> pl_orig; // 创建原始点云
  pcl::fromROSMsg(*msg, pl_orig); // 从ROS消息转换
  int plsize = pl_orig.size(); // 获取点云大小
  pl_corn.reserve(plsize); // 预留空间
  pl_surf.reserve(plsize); // 预留空间

  double time_stamp = msg->header.stamp.toSec(); // 获取时间戳
  // cout << "===================================" << endl;
  // printf("Pt size = %d, N_SCANS = %d\r\n", plsize, N_SCANS);
  for (int i = 0; i < pl_orig.points.size(); i++)
  {
    if (i % point_filter_num != 0) continue; // 根据过滤数量选择性处理

    double range = pl_orig.points[i].x * pl_orig.points[i].x + pl_orig.points[i].y * pl_orig.points[i].y + pl_orig.points[i].z * pl_orig.points[i].z;

    if (range < blind_sqr) continue; // 如果在盲区则跳过

    Eigen::Vector3d pt_vec; // 创建3D点向量
    PointType added_pt; // 创建添加的点
    added_pt.x = pl_orig.points[i].x; // 设置X坐标
    added_pt.y = pl_orig.points[i].y; // 设置Y坐标
    added_pt.z = pl_orig.points[i].z; // 设置Z坐标
    added_pt.normal_x = pl_orig.points[i].r; // 设置法向量X
    added_pt.normal_y = pl_orig.points[i].g; // 设置法向量Y
    added_pt.normal_z = pl_orig.points[i].b; // 设置法向量Z

    added_pt.curvature = 0.0; // 设置曲率
    pl_surf.points.push_back(added_pt); // 添加到表面点云
  }

  cout << "pl size:: " << pl_orig.points.size() << endl; // 输出原始点云大小
  // pub_func(pl_surf, pub_full, msg->header.stamp);
  // pub_func(pl_surf, pub_corn, msg->header.stamp);
}

// 处理OUST64类型的点云数据
void Preprocess::oust64_handler(const sensor_msgs::PointCloud2::ConstPtr &msg)
{
  pl_surf.clear(); // 清空表面点云
  pl_corn.clear(); // 清空边缘点云
  pl_full.clear(); // 清空全量点云
  pcl::PointCloud<ouster_ros::Point> pl_orig; // 创建原始点云
  pcl::fromROSMsg(*msg, pl_orig); // 从ROS消息转换
  int plsize = pl_orig.size(); // 获取点云大小
  pl_corn.reserve(plsize); // 预留空间
  pl_surf.reserve(plsize); // 预留空间
  if (feature_enabled) // 如果启用特征提取
  {
    for (int i = 0; i < N_SCANS; i++)
    {
      pl_buff[i].clear(); // 清空扫描缓冲区
      pl_buff[i].reserve(plsize); // 预留空间
    }

    for (uint i = 0; i < plsize; i++)
    {
      double range =
          pl_orig.points[i].x * pl_orig.points[i].x + pl_orig.points[i].y * pl_orig.points[i].y + pl_orig.points[i].z * pl_orig.points[i].z;
      if (range < blind_sqr) continue; // 如果在盲区则跳过
      Eigen::Vector3d pt_vec; // 创建3D点向量
      PointType added_pt; // 创建添加的点
      added_pt.x = pl_orig.points[i].x; // 设置X坐标
      added_pt.y = pl_orig.points[i].y; // 设置Y坐标
      added_pt.z = pl_orig.points[i].z; // 设置Z坐标
      added_pt.intensity = pl_orig.points[i].intensity; // 设置强度
      added_pt.normal_x = 0; // 设置法向量X
      added_pt.normal_y = 0; // 设置法向量Y
      added_pt.normal_z = 0; // 设置法向量Z
      double yaw_angle = atan2(added_pt.y, added_pt.x) * 57.3; // 计算偏航角
      if (yaw_angle >= 180.0) yaw_angle -= 360.0; // 处理角度范围
      if (yaw_angle <= -180.0) yaw_angle += 360.0; // 处理角度范围

      added_pt.curvature = pl_orig.points[i].t / 1e6; // 设置曲率
      if (pl_orig.points[i].ring < N_SCANS) { pl_buff[pl_orig.points[i].ring].push_back(added_pt); } // 添加到相应的扫描缓冲区
    }

    for (int j = 0; j < N_SCANS; j++)
    {
      PointCloudXYZI &pl = pl_buff[j]; // 获取当前扫描线的点云
      int linesize = pl.size(); // 获取点云大小
      vector<orgtype> &types = typess[j]; // 获取当前扫描线的类型信息
      types.clear(); // 清空类型信息
      types.resize(linesize); // 调整类型信息大小
      linesize--; // 调整点云大小（最后一个点单独处理）
      for (uint i = 0; i < linesize; i++)
      {
        types[i].range = sqrt(pl[i].x * pl[i].x + pl[i].y * pl[i].y); // 计算范围
        vx = pl[i].x - pl[i + 1].x; // 计算X差
        vy = pl[i].y - pl[i + 1].y; // 计算Y差
        vz = pl[i].z - pl[i + 1].z; // 计算Z差
        types[i].dista = vx * vx + vy * vy + vz * vz; // 计算距离
      }
      types[linesize].range = sqrt(pl[linesize].x * pl[linesize].x + pl[linesize].y * pl[linesize].y); // 处理最后一个点的范围
      give_feature(pl, types); // 提取特征
    }
  }
  else // 如果不启用特征提取
  {
    double time_stamp = msg->header.stamp.toSec(); // 获取时间戳
    // cout << "===================================" << endl;
    // printf("Pt size = %d, N_SCANS = %d\r\n", plsize, N_SCANS);
    for (int i = 0; i < pl_orig.points.size(); i++)
    {
      if (i % point_filter_num != 0) continue; // 根据过滤数量选择性处理

      double range =
          pl_orig.points[i].x * pl_orig.points[i].x + pl_orig.points[i].y * pl_orig.points[i].y + pl_orig.points[i].z * pl_orig.points[i].z;

      if (range < blind_sqr) continue; // 如果在盲区则跳过

      Eigen::Vector3d pt_vec; // 创建3D点向量
      PointType added_pt; // 创建添加的点
      added_pt.x = pl_orig.points[i].x; // 设置X坐标
      added_pt.y = pl_orig.points[i].y; // 设置Y坐标
      added_pt.z = pl_orig.points[i].z; // 设置Z坐标
      added_pt.intensity = pl_orig.points[i].intensity; // 设置强度
      added_pt.normal_x = 0; // 设置法向量X
      added_pt.normal_y = 0; // 设置法向量Y
      added_pt.normal_z = 0; // 设置法向量Z
      double yaw_angle = atan2(added_pt.y, added_pt.x) * 57.3; // 计算偏航角
      if (yaw_angle >= 180.0) yaw_angle -= 360.0; // 处理角度范围
      if (yaw_angle <= -180.0) yaw_angle += 360.0; // 处理角度范围

      added_pt.curvature = pl_orig.points[i].t / 1e6; // 设置曲率

      // cout<<added_pt.curvature<<endl;

      pl_surf.points.push_back(added_pt); // 添加到表面点云
    }
  }
  // pub_func(pl_surf, pub_full, msg->header.stamp);
  // pub_func(pl_surf, pub_corn, msg->header.stamp);
}

// 定义最大扫描行数
#define MAX_LINE_NUM 64

// 处理VELODYNE类型的点云数据
void Preprocess::velodyne_handler(const sensor_msgs::PointCloud2::ConstPtr &msg)
{
  pl_surf.clear(); // 清空表面点云
  pl_corn.clear(); // 清空边缘点云
  pl_full.clear(); // 清空全量点云

  pcl::PointCloud<velodyne_ros::Point> pl_orig; // 创建原始点云
  pcl::fromROSMsg(*msg, pl_orig); // 从ROS消息转换
  int plsize = pl_orig.points.size(); // 获取点云大小
  pl_surf.reserve(plsize); // 预留空间

  bool is_first[MAX_LINE_NUM]; // 是否为第一点标志
  double yaw_fp[MAX_LINE_NUM] = {0};     // 第一扫描点的偏航角
  double omega_l = 3.61;                 // 扫描角速度
  float yaw_last[MAX_LINE_NUM] = {0.0};  // 上一个扫描点的偏航角
  float time_last[MAX_LINE_NUM] = {0.0}; // 上一个偏移时间

  if (pl_orig.points[plsize - 1].t > 0) { given_offset_time = true; } // 根据最后一个点的时间戳判断是否给定偏移时间
  else
  {
    given_offset_time = false; // 未给定偏移时间
    memset(is_first, true, sizeof(is_first)); // 初始化第一点标志
    double yaw_first = atan2(pl_orig.points[0].y, pl_orig.points[0].x) * 57.29578; // 计算第一点的偏航角
    double yaw_end = yaw_first; // 记录结束偏航角
    int layer_first = pl_orig.points[0].ring; // 记录第一点的扫描层
    for (uint i = plsize - 1; i > 0; i--)
    {
      if (pl_orig.points[i].ring == layer_first) // 查找同一层的最后一个点
      {
        yaw_end = atan2(pl_orig.points[i].y, pl_orig.points[i].x) * 57.29578; // 计算结束偏航角
        break;
      }
    }
  }

  if (feature_enabled) // 如果启用特征提取
  {
    for (int i = 0; i < N_SCANS; i++)
    {
      pl_buff[i].clear(); // 清空扫描缓冲区
      pl_buff[i].reserve(plsize); // 预留空间
    }

    for (int i = 0; i < plsize; i++)
    {
      PointType added_pt; // 创建添加的点
      added_pt.normal_x = 0; // 设置法向量X
      added_pt.normal_y = 0; // 设置法向量Y
      added_pt.normal_z = 0; // 设置法向量Z
      int layer = pl_orig.points[i].ring; // 获取当前点的扫描层
      if (layer >= N_SCANS) continue; // 如果扫描层超出范围则跳过
      added_pt.x = pl_orig.points[i].x; // 设置X坐标
      added_pt.y = pl_orig.points[i].y; // 设置Y坐标
      added_pt.z = pl_orig.points[i].z; // 设置Z坐标
      added_pt.intensity = pl_orig.points[i].intensity; // 设置强度
      added_pt.curvature = pl_orig.points[i].t / 1000.0; // 设置曲率（单位：ms）

      if (!given_offset_time) // 如果未给定偏移时间
      {
        double yaw_angle = atan2(added_pt.y, added_pt.x) * 57.2957; // 计算偏航角
        if (is_first[layer]) // 如果是该层的第一点
        {
          // printf("layer: %d; is first: %d", layer, is_first[layer]);
          yaw_fp[layer] = yaw_angle; // 记录第一点的偏航角
          is_first[layer] = false; // 标记为非第一点
          added_pt.curvature = 0.0; // 设置曲率为0
          yaw_last[layer] = yaw_angle; // 记录上一个偏航角
          time_last[layer] = added_pt.curvature; // 记录上一个时间
          continue; // 继续处理下一个点
        }

        // 计算偏移时间
        if (yaw_angle <= yaw_fp[layer]) { added_pt.curvature = (yaw_fp[layer] - yaw_angle) / omega_l; }
        else { added_pt.curvature = (yaw_fp[layer] - yaw_angle + 360.0) / omega_l; }

        if (added_pt.curvature < time_last[layer]) added_pt.curvature += 360.0 / omega_l; // 确保曲率值递增

        yaw_last[layer] = yaw_angle; // 更新上一个偏航角
        time_last[layer] = added_pt.curvature; // 更新上一个时间
      }

      pl_buff[layer].points.push_back(added_pt); // 将当前点添加到相应的扫描缓冲区
    }

    for (int j = 0; j < N_SCANS; j++)
    {
      PointCloudXYZI &pl = pl_buff[j]; // 获取当前扫描线的点云
      int linesize = pl.size(); // 获取点云大小
      if (linesize < 2) continue; // 如果点数不足则跳过
      vector<orgtype> &types = typess[j]; // 获取当前扫描线的类型信息
      types.clear(); // 清空类型信息
      types.resize(linesize); // 调整类型信息大小
      linesize--; // 调整点云大小（最后一个点单独处理）
      for (uint i = 0; i < linesize; i++)
      {
        types[i].range = sqrt(pl[i].x * pl[i].x + pl[i].y * pl[i].y); // 计算范围
        vx = pl[i].x - pl[i + 1].x; // 计算X差
        vy = pl[i].y - pl[i + 1].y; // 计算Y差
        vz = pl[i].z - pl[i + 1].z; // 计算Z差
        types[i].dista = vx * vx + vy * vy + vz * vz; // 计算距离
      }
      types[linesize].range = sqrt(pl[linesize].x * pl[linesize].x + pl[linesize].y * pl[linesize].y); // 处理最后一个点的范围
      give_feature(pl, types); // 提取特征
    }
  }
  else // 如果不启用特征提取
  {
    for (int i = 0; i < plsize; i++)
    {
      PointType added_pt; // 创建添加的点
      // cout<<"!!!!!!"<<i<<" "<<plsize<<endl;

      added_pt.normal_x = 0; // 设置法向量X
      added_pt.normal_y = 0; // 设置法向量Y
      added_pt.normal_z = 0; // 设置法向量Z
      added_pt.x = pl_orig.points[i].x; // 设置X坐标
      added_pt.y = pl_orig.points[i].y; // 设置Y坐标
      added_pt.z = pl_orig.points[i].z; // 设置Z坐标
      added_pt.intensity = pl_orig.points[i].intensity; // 设置强度
      added_pt.curvature = pl_orig.points[i].t / 1000.0; // 设置曲率（单位：ms）

      if (!given_offset_time) // 如果未给定偏移时间
      {
        int layer = pl_orig.points[i].ring; // 获取当前点的扫描层
        double yaw_angle = atan2(added_pt.y, added_pt.x) * 57.2957; // 计算偏航角

        if (is_first[layer]) // 如果是该层的第一点
        {
          // printf("layer: %d; is first: %d", layer, is_first[layer]);
          yaw_fp[layer] = yaw_angle; // 记录第一点的偏航角
          is_first[layer] = false; // 标记为非第一点
          added_pt.curvature = 0.0; // 设置曲率为0
          yaw_last[layer] = yaw_angle; // 记录上一个偏航角
          time_last[layer] = added_pt.curvature; // 记录上一个时间
          continue; // 继续处理下一个点
        }

        // 计算偏移时间
        if (yaw_angle <= yaw_fp[layer]) { added_pt.curvature = (yaw_fp[layer] - yaw_angle) / omega_l; }
        else { added_pt.curvature = (yaw_fp[layer] - yaw_angle + 360.0) / omega_l; }

        if (added_pt.curvature < time_last[layer]) added_pt.curvature += 360.0 / omega_l; // 确保曲率值递增

        // added_pt.curvature = pl_orig.points[i].t;

        yaw_last[layer] = yaw_angle; // 更新上一个偏航角
        time_last[layer] = added_pt.curvature; // 更新上一个时间
      }

      // if(i==(plsize-1))  printf("index: %d layer: %d, yaw: %lf, offset-time:
      // %lf, condition: %d\n", i, layer, yaw_angle, added_pt.curvature,
      // prints);
      if (i % point_filter_num == 0) // 根据过滤数量选择性加入点云
      {
        if (added_pt.x * added_pt.x + added_pt.y * added_pt.y + added_pt.z * added_pt.z > blind_sqr)
        {
          pl_surf.points.push_back(added_pt); // 添加到表面点云
          // printf("time mode: %d time: %d \n", given_offset_time,
          // pl_orig.points[i].t);
        }
      }
    }
  }
  // pub_func(pl_surf, pub_full, msg->header.stamp);
  // pub_func(pl_surf, pub_surf, msg->header.stamp);
  // pub_func(pl_surf, pub_corn, msg->header.stamp);
}

// 处理Pandar128类型的点云数据
void Preprocess::Pandar128_handler(const sensor_msgs::PointCloud2::ConstPtr &msg)
{
  pl_surf.clear(); // 清空表面点云

  pcl::PointCloud<Pandar128_ros::Point> pl_orig; // 创建原始点云
  pcl::fromROSMsg(*msg, pl_orig); // 从ROS消息转换
  int plsize = pl_orig.points.size(); // 获取点云大小
  pl_surf.reserve(plsize); // 预留空间

  // double time_head = pl_orig.points[0].timestamp; // 获取时间头
  for (int i = 0; i < plsize; i++)
  {
    PointType added_pt; // 创建添加的点

    added_pt.normal_x = 0; // 设置法向量X
    added_pt.normal_y = 0; // 设置法向量Y
    added_pt.normal_z = 0; // 设置法向量Z
    added_pt.x = pl_orig.points[i].x; // 设置X坐标
    added_pt.y = pl_orig.points[i].y; // 设置Y坐标
    added_pt.z = pl_orig.points[i].z; // 设置Z坐标
    added_pt.curvature = pl_orig.points[i].timestamp * 1000.f; // 设置曲率

    if (i % point_filter_num == 0) // 根据过滤数量选择性加入点云
    {
      if (added_pt.x * added_pt.x + added_pt.y * added_pt.y + added_pt.z * added_pt.z > blind_sqr)
      {
        pl_surf.points.push_back(added_pt); // 添加到表面点云
        // printf("time mode: %d time: %d \n", given_offset_time,
        // pl_orig.points[i].t);
      }
    }
  }

  // 定义比较函数
  auto comparePoints = [](const PointType& a, const PointType& b) -> bool
  {
    return a.curvature < b.curvature; // 按照曲率排序
  };
  
  // 使用比较函数排序点云
  std::sort(pl_surf.points.begin(), pl_surf.points.end(), comparePoints);
  
  // cout << GREEN << "pl_surf.points[0].timestamp: " << pl_surf.points[0].curvature << RESET << endl;
  // cout << GREEN << "pl_surf.points[1000].timestamp: " << pl_surf.points[1000].curvature << RESET << endl;
  // cout << GREEN << "pl_surf.points[5000].timestamp: " << pl_surf.points[5000].curvature << RESET << endl;
  // cout << GREEN << "pl_surf.points[10000].timestamp: " << pl_surf.points[10000].curvature << RESET << endl;
  // cout << GREEN << "pl_surf.points[20000].timestamp: " << pl_surf.points[20000].curvature << RESET << endl;
  // cout << GREEN << "pl_surf.points[30000].timestamp: " << pl_surf.points[30000].curvature << RESET << endl;
  // cout << GREEN << "pl_surf.points[31000].timestamp: " << pl_surf.points[31000].curvature << RESET << endl;
}
// 处理XT32类型的点云数据
void Preprocess::xt32_handler(const sensor_msgs::PointCloud2::ConstPtr &msg)
{
  pl_surf.clear(); // 清空表面点云
  pl_corn.clear(); // 清空边缘点云
  pl_full.clear(); // 清空全量点云

  pcl::PointCloud<xt32_ros::Point> pl_orig; // 创建原始点云
  pcl::fromROSMsg(*msg, pl_orig); // 从ROS消息转换
  int plsize = pl_orig.points.size(); // 获取点云大小
  pl_surf.reserve(plsize); // 预留空间

  bool is_first[MAX_LINE_NUM]; // 第一点标志数组
  double yaw_fp[MAX_LINE_NUM] = {0};     // 第一扫描点的偏航角
  double omega_l = 3.61;                 // 扫描角速度
  float yaw_last[MAX_LINE_NUM] = {0.0};  // 上一个扫描点的偏航角
  float time_last[MAX_LINE_NUM] = {0.0}; // 上一个偏移时间

  // 判断是否给定偏移时间
  if (pl_orig.points[plsize - 1].timestamp > 0) { 
    given_offset_time = true; 
  } else {
    given_offset_time = false; // 未给定偏移时间
    memset(is_first, true, sizeof(is_first)); // 初始化第一点标志
    double yaw_first = atan2(pl_orig.points[0].y, pl_orig.points[0].x) * 57.29578; // 计算第一点的偏航角
    double yaw_end = yaw_first; // 记录结束偏航角
    int layer_first = pl_orig.points[0].ring; // 记录第一点的扫描层
    for (uint i = plsize - 1; i > 0; i--) {
      if (pl_orig.points[i].ring == layer_first) {
        yaw_end = atan2(pl_orig.points[i].y, pl_orig.points[i].x) * 57.29578; // 计算结束偏航角
        break;
      }
    }
  }

  double time_head = pl_orig.points[0].timestamp; // 获取时间头

  if (feature_enabled) // 如果启用特征提取
  {
    for (int i = 0; i < N_SCANS; i++) {
      pl_buff[i].clear(); // 清空扫描缓冲区
      pl_buff[i].reserve(plsize); // 预留空间
    }

    for (int i = 0; i < plsize; i++) // 遍历每个点
    {
      PointType added_pt; // 创建添加的点
      added_pt.normal_x = 0; // 设置法向量X
      added_pt.normal_y = 0; // 设置法向量Y
      added_pt.normal_z = 0; // 设置法向量Z
      int layer = pl_orig.points[i].ring; // 获取当前点的扫描层
      if (layer >= N_SCANS) continue; // 如果扫描层超出范围则跳过
      added_pt.x = pl_orig.points[i].x; // 设置X坐标
      added_pt.y = pl_orig.points[i].y; // 设置Y坐标
      added_pt.z = pl_orig.points[i].z; // 设置Z坐标
      added_pt.intensity = pl_orig.points[i].intensity; // 设置强度
      added_pt.curvature = pl_orig.points[i].timestamp / 1000.0; // 设置曲率（单位：ms）

      if (!given_offset_time) // 如果未给定偏移时间
      {
        double yaw_angle = atan2(added_pt.y, added_pt.x) * 57.2957; // 计算偏航角
        if (is_first[layer]) // 如果是该层的第一点
        {
          // printf("layer: %d; is first: %d", layer, is_first[layer]);
          yaw_fp[layer] = yaw_angle; // 记录第一点的偏航角
          is_first[layer] = false; // 标记为非第一点
          added_pt.curvature = 0.0; // 设置曲率为0
          yaw_last[layer] = yaw_angle; // 记录上一个偏航角
          time_last[layer] = added_pt.curvature; // 记录上一个时间
          continue; // 继续处理下一个点
        }

        // 计算偏移时间
        if (yaw_angle <= yaw_fp[layer]) { 
          added_pt.curvature = (yaw_fp[layer] - yaw_angle) / omega_l; 
        } else { 
          added_pt.curvature = (yaw_fp[layer] - yaw_angle + 360.0) / omega_l; 
        }

        if (added_pt.curvature < time_last[layer]) 
          added_pt.curvature += 360.0 / omega_l; // 确保曲率值递增

        yaw_last[layer] = yaw_angle; // 更新上一个偏航角
        time_last[layer] = added_pt.curvature; // 更新上一个时间
      }

      pl_buff[layer].points.push_back(added_pt); // 将当前点添加到相应的扫描缓冲区
    }

    // 遍历每条扫描线
    for (int j = 0; j < N_SCANS; j++)
    {
      PointCloudXYZI &pl = pl_buff[j]; // 获取当前扫描线的点云
      int linesize = pl.size(); // 获取点云大小
      if (linesize < 2) continue; // 如果点数不足则跳过
      vector<orgtype> &types = typess[j]; // 获取当前扫描线的类型信息
      types.clear(); // 清空类型信息
      types.resize(linesize); // 调整类型信息大小
      linesize--; // 调整点云大小（最后一个点单独处理）

      for (uint i = 0; i < linesize; i++)
      {
        types[i].range = sqrt(pl[i].x * pl[i].x + pl[i].y * pl[i].y); // 计算范围
        vx = pl[i].x - pl[i + 1].x; // 计算X差
        vy = pl[i].y - pl[i + 1].y; // 计算Y差
        vz = pl[i].z - pl[i + 1].z; // 计算Z差
        types[i].dista = vx * vx + vy * vy + vz * vz; // 计算距离
      }
      types[linesize].range = sqrt(pl[linesize].x * pl[linesize].x + pl[linesize].y * pl[linesize].y); // 处理最后一个点的范围
      give_feature(pl, types); // 提取特征
    }
  }
  else // 如果不启用特征提取
  {
    for (int i = 0; i < plsize; i++)
    {
      PointType added_pt; // 创建添加的点
      // cout<<"!!!!!!"<<i<<" "<<plsize<<endl;

      added_pt.normal_x = 0; // 设置法向量X
      added_pt.normal_y = 0; // 设置法向量Y
      added_pt.normal_z = 0; // 设置法向量Z
      added_pt.x = pl_orig.points[i].x; // 设置X坐标
      added_pt.y = pl_orig.points[i].y; // 设置Y坐标
      added_pt.z = pl_orig.points[i].z; // 设置Z坐标
      added_pt.intensity = pl_orig.points[i].intensity; // 设置强度
      added_pt.curvature = (pl_orig.points[i].timestamp - time_head) * 1000.f; // 设置曲率

      // printf("added_pt.curvature: %lf %lf \n", added_pt.curvature,
      // pl_orig.points[i].timestamp);

      if (i % point_filter_num == 0) // 根据过滤数量选择性加入点云
      {
        if (added_pt.x * added_pt.x + added_pt.y * added_pt.y + added_pt.z * added_pt.z > blind_sqr)
        {
          pl_surf.points.push_back(added_pt); // 添加到表面点云
          // printf("time mode: %d time: %d \n", given_offset_time,
          // pl_orig.points[i].t);
        }
      }
    }
  }
  // pub_func(pl_surf, pub_full, msg->header.stamp);
  // pub_func(pl_surf, pub_surf, msg->header.stamp);
  // pub_func(pl_surf, pub_corn, msg->header.stamp);
}

// 提取特征
void Preprocess::give_feature(pcl::PointCloud<PointType> &pl, vector<orgtype> &types)
{
  int plsize = pl.size(); // 点云大小
  int plsize2; // 第二个点云大小
  if (plsize == 0)
  {
    printf("something wrong\n");
    return; // 如果点云为空，直接返回
  }
  uint head = 0; // 头部索引

  // 寻找第一个有效点
  while (types[head].range < blind_sqr)
  {
    head++; // 移动头部索引
  }

  // Surf
  plsize2 = (plsize > group_size) ? (plsize - group_size) : 0; // 计算第二个点云的大小

  Eigen::Vector3d curr_direct(Eigen::Vector3d::Zero()); // 当前方向
  Eigen::Vector3d last_direct(Eigen::Vector3d::Zero()); // 上一方向

  uint i_nex = 0, i2; // 下一点索引
  uint last_i = 0; // 上一个点索引
  uint last_i_nex = 0; // 上一个点的下一点索引
  int last_state = 0; // 上一个状态
  int plane_type; // 平面类型

  // 遍历点云
  for (uint i = head; i < plsize2; i++)
  {
    if (types[i].range < blind_sqr) { continue; } // 如果在盲区则跳过

    i2 = i; // 记录当前索引

    plane_type = plane_judge(pl, types, i, i_nex, curr_direct); // 判断平面类型

    if (plane_type == 1) // 如果是平面类型
    {
      for (uint j = i; j <= i_nex; j++)
      {
        if (j != i && j != i_nex) { types[j].ftype = Real_Plane; } // 设置真实平面
        else { types[j].ftype = Poss_Plane; } // 设置可能平面
      }

      // if(last_state==1 && fabs(last_direct.sum())>0.5)
      if (last_state == 1 && last_direct.norm() > 0.1) // 检查上一个状态和方向
      {
        double mod = last_direct.transpose() * curr_direct; // 计算方向的点积
        if (mod > -0.707 && mod < 0.707) { types[i].ftype = Edge_Plane; } // 如果方向相差较大，标记为边缘平面
        else { types[i].ftype = Real_Plane; } // 否则标记为真实平面
      }

      i = i_nex - 1; // 更新当前索引
      last_state = 1; // 更新状态
    }
    else // 如果不是平面类型
    {
      i = i_nex; // 更新当前索引
      last_state = 0; // 更新状态
    }

    last_i = i2; // 更新上一个点索引
    last_i_nex = i_nex; // 更新上一个点的下一点索引
    last_direct = curr_direct; // 更新上一个方向
  }

  plsize2 = plsize > 3 ? plsize - 3 : 0; // 计算第二个点云的大小
  for (uint i = head + 3; i < plsize2; i++)
  {
    if (types[i].range < blind_sqr || types[i].ftype >= Real_Plane) { continue; } // 如果在盲区或已标记为真实平面则跳过

    if (types[i - 1].dista < 1e-16 || types[i].dista < 1e-16) { continue; } // 如果距离无效则跳过

    Eigen::Vector3d vec_a(pl[i].x, pl[i].y, pl[i].z); // 当前点向量
    Eigen::Vector3d vecs[2]; // 存储前后点向量

    for (int j = 0; j < 2; j++)
    {
      int m = -1; // 初始化索引
      if (j == 1) { m = 1; } // 如果是第二个点，则索引为1

      if (types[i + m].range < blind_sqr) // 如果在盲区
      {
        if (types[i].range > inf_bound) { types[i].edj[j] = Nr_inf; } // 设置边缘状态
        else { types[i].edj[j] = Nr_blind; }
        continue; // 继续处理下一个点
      }

      vecs[j] = Eigen::Vector3d(pl[i + m].x, pl[i + m].y, pl[i + m].z); // 设置向量
      vecs[j] = vecs[j] - vec_a; // 计算向量差

      types[i].angle[j] = vec_a.dot(vecs[j]) / vec_a.norm() / vecs[j].norm(); // 计算夹角
      if (types[i].angle[j] < jump_up_limit) { types[i].edj[j] = Nr_180; } // 设置边缘状态
      else if (types[i].angle[j] > jump_down_limit) { types[i].edj[j] = Nr_zero; }
    }

    // 计算交点
    types[i].intersect = vecs[Prev].dot(vecs[Next]) / vecs[Prev].norm() / vecs[Next].norm();
    if (types[i].edj[Prev] == Nr_nor && types[i].edj[Next] == Nr_zero && types[i].dista > 0.0225 && types[i].dista > 4 * types[i - 1].dista)
    {
      if (types[i].intersect > cos160)
      {
        if (edge_jump_judge(pl, types, i, Prev)) { types[i].ftype = Edge_Jump; } // 判断边缘跳跃
      }
    }
    else if (types[i].edj[Prev] == Nr_zero && types[i].edj[Next] == Nr_nor && types[i - 1].dista > 0.0225 && types[i - 1].dista > 4 * types[i].dista)
    {
      if (types[i].intersect > cos160)
      {
        if (edge_jump_judge(pl, types, i, Next)) { types[i].ftype = Edge_Jump; } // 判断边缘跳跃
      }
    }
    else if (types[i].edj[Prev] == Nr_nor && types[i].edj[Next] == Nr_inf)
    {
      if (edge_jump_judge(pl, types, i, Prev)) { types[i].ftype = Edge_Jump; } // 判断边缘跳跃
    }
    else if (types[i].edj[Prev] == Nr_inf && types[i].edj[Next] == Nr_nor)
    {
      if (edge_jump_judge(pl, types, i, Next)) { types[i].ftype = Edge_Jump; } // 判断边缘跳跃
    }
    else if (types[i].edj[Prev] > Nr_nor && types[i].edj[Next] > Nr_nor)
    {
      if (types[i].ftype == Nor) { types[i].ftype = Wire; } // 标记为线
    }
  }

  plsize2 = plsize - 1; // 更新点云大小
  double ratio; // 比率
  for (uint i = head + 1; i < plsize2; i++)
  {
    if (types[i].range < blind_sqr || types[i - 1].range < blind_sqr || types[i + 1].range < blind_sqr) { continue; } // 如果在盲区则跳过

    if (types[i - 1].dista < 1e-8 || types[i].dista < 1e-8) { continue; } // 如果距离无效则跳过

    if (types[i].ftype == Nor) // 如果是普通类型
    {
      if (types[i - 1].dista > types[i].dista) { ratio = types[i - 1].dista / types[i].dista; } // 计算比率
      else { ratio = types[i].dista / types[i - 1].dista; }

      if (types[i].intersect < smallp_intersect && ratio < smallp_ratio) // 判断小点交集和比率
      {
        if (types[i - 1].ftype == Nor) { types[i - 1].ftype = Real_Plane; } // 标记为真实平面
        if (types[i + 1].ftype == Nor) { types[i + 1].ftype = Real_Plane; } // 标记为真实平面
        types[i].ftype = Real_Plane; // 标记为真实平面
      }
    }
  }

  int last_surface = -1; // 上一个表面索引
  for (uint j = head; j < plsize; j++)
  {
    if (types[j].ftype == Poss_Plane || types[j].ftype == Real_Plane) // 如果是可能平面或真实平面
    {
      if (last_surface == -1) { last_surface = j; } // 记录上一个表面索引

      if (j == uint(last_surface + point_filter_num - 1)) // 如果达到过滤数量
      {
        PointType ap; // 创建添加的点
        ap.x = pl[j].x; // 设置X坐标
        ap.y = pl[j].y; // 设置Y坐标
        ap.z = pl[j].z; // 设置Z坐标
        ap.curvature = pl[j].curvature; // 设置曲率
        pl_surf.push_back(ap); // 添加到表面点云

        last_surface = -1; // 重置上一个表面索引
      }
    }
    else // 如果不是平面
    {
      if (types[j].ftype == Edge_Jump || types[j].ftype == Edge_Plane) { pl_corn.push_back(pl[j]); } // 如果是边缘跳跃或边缘平面，添加到边缘点云
      if (last_surface != -1) // 如果存在上一个表面
      {
        PointType ap; // 创建添加的点
        for (uint k = last_surface; k < j; k++)
        {
          ap.x += pl[k].x; // 累加X坐标
          ap.y += pl[k].y; // 累加Y坐标
          ap.z += pl[k].z; // 累加Z坐标
          ap.curvature += pl[k].curvature; // 累加曲率
        }
        ap.x /= (j - last_surface); // 计算平均X坐标
        ap.y /= (j - last_surface); // 计算平均Y坐标
        ap.z /= (j - last_surface); // 计算平均Z坐标
        ap.curvature /= (j - last_surface); // 计算平均曲率
        pl_surf.push_back(ap); // 添加到表面点云
      }
      last_surface = -1; // 重置上一个表面索引
    }
  }
}

// 发布点云数据
void Preprocess::pub_func(PointCloudXYZI &pl, const ros::Time &ct)
{
  pl.height = 1; // 设置点云高度
  pl.width = pl.size(); // 设置点云宽度
  sensor_msgs::PointCloud2 output; // 创建ROS消息
  pcl::toROSMsg(pl, output); // 转换为ROS消息
  output.header.frame_id = "livox"; // 设置帧ID
  output.header.stamp = ct; // 设置时间戳
}

// 判断平面类型
int Preprocess::plane_judge(const PointCloudXYZI &pl, vector<orgtype> &types, uint i_cur, uint &i_nex, Eigen::Vector3d &curr_direct)
{
  double group_dis = disA * types[i_cur].range + disB; // 计算分组距离
  group_dis = group_dis * group_dis; // 平方

  double two_dis; // 两点间的距离
  vector<double> disarr; // 距离数组
  disarr.reserve(20); // 预留空间

  for (i_nex = i_cur; i_nex < i_cur + group_size; i_nex++) // 遍历分组大小
  {
    if (types[i_nex].range < blind_sqr)
    {
      curr_direct.setZero(); // 设置当前方向为零
      return 2; // 返回状态
    }
    disarr.push_back(types[i_nex].dista); // 添加距离到数组
  }

  for (;;)
  {
    if ((i_cur >= pl.size()) || (i_nex >= pl.size())) break; // 如果超出范围则跳出

    if (types[i_nex].range < blind_sqr)
    {
      curr_direct.setZero(); // 设置当前方向为零
      return 2; // 返回状态
    }
    vx = pl[i_nex].x - pl[i_cur].x; // 计算X差
    vy = pl[i_nex].y - pl[i_cur].y; // 计算Y差
    vz = pl[i_nex].z - pl[i_cur].z; // 计算Z差
    two_dis = vx * vx + vy * vy + vz * vz; // 计算距离的平方
    if (two_dis >= group_dis) { break; } // 如果超出分组距离则跳出
    disarr.push_back(types[i_nex].dista); // 添加距离到数组
    i_nex++; // 移动到下一个点
  }

  double leng_wid = 0; // 长宽
  double v1[3], v2[3]; // 向量
  for (uint j = i_cur + 1; j < i_nex; j++) // 遍历当前到下一点之间的点
  {
    if ((j >= pl.size()) || (i_cur >= pl.size())) break; // 如果超出范围则跳出
    v1[0] = pl[j].x - pl[i_cur].x; // 计算向量X
    v1[1] = pl[j].y - pl[i_cur].y; // 计算向量Y
    v1[2] = pl[j].z - pl[i_cur].z; // 计算向量Z

    v2[0] = v1[1] * vz - vy * v1[2]; // 计算叉积
    v2[1] = v1[2] * vx - v1[0] * vz;
    v2[2] = v1[0] * vy - vx * v1[1];

    double lw = v2[0] * v2[0] + v2[1] * v2[1] + v2[2] * v2[2]; // 计算长度平方
    if (lw > leng_wid) { leng_wid = lw; } // 更新最大长度
  }

  if ((two_dis * two_dis / leng_wid) < p2l_ratio) // 判断是否为平面
  {
    curr_direct.setZero(); // 设置当前方向为零
    return 0; // 返回状态
  }

  uint disarrsize = disarr.size(); // 距离数组大小
  for (uint j = 0; j < disarrsize - 1; j++) // 冒泡排序
  {
    for (uint k = j + 1; k < disarrsize; k++)
    {
      if (disarr[j] < disarr[k])
      {
        leng_wid = disarr[j]; // 交换
        disarr[j] = disarr[k];
        disarr[k] = leng_wid;
      }
    }
  }

  if (disarr[disarr.size() - 2] < 1e-16) // 判断是否有效
  {
    curr_direct.setZero(); // 设置当前方向为零
    return 0; // 返回状态
  }

  // 根据激光雷达类型进行判断
  if (lidar_type == AVIA)
  {
    double dismax_mid = disarr[0] / disarr[disarrsize / 2]; // 计算最大中间比率
    double dismid_min = disarr[disarrsize / 2] / disarr[disarrsize - 2]; // 计算中间最小比率

    if (dismax_mid >= limit_maxmid || dismid_min >= limit_midmin) // 判断是否超出限制
    {
      curr_direct.setZero(); // 设置当前方向为零
      return 0; // 返回状态
    }
  }
  else
  {
    double dismax_min = disarr[0] / disarr[disarrsize - 2]; // 计算最大最小比率
    if (dismax_min >= limit_maxmin) // 判断是否超出限制
    {
      curr_direct.setZero(); // 设置当前方向为零
      return 0; // 返回状态
    }
  }

  curr_direct << vx, vy, vz; // 设置当前方向
  curr_direct.normalize(); // 归一化
  return 1; // 返回状态
}

// 判断边缘跳跃
bool Preprocess::edge_jump_judge(const PointCloudXYZI &pl, vector<orgtype> &types, uint i, Surround nor_dir)
{
  if (nor_dir == 0) // 前一个方向
  {
    if (types[i - 1].range < blind_sqr || types[i - 2].range < blind_sqr) { return false; } // 如果在盲区则返回false
  }
  else if (nor_dir == 1) // 后一个方向
  {
    if (types[i + 1].range < blind_sqr || types[i + 2].range < blind_sqr) { return false; } // 如果在盲区则返回false
  }
  double d1 = types[i + nor_dir - 1].dista; // 获取距离
  double d2 = types[i + 3 * nor_dir - 2].dista; // 获取距离
  double d;

  if (d1 < d2) // 如果d1小于d2，则交换
  {
    d = d1;
    d1 = d2;
    d2 = d;
  }

  d1 = sqrt(d1); // 计算平方根
  d2 = sqrt(d2); // 计算平方根

  if (d1 > edgea * d2 || (d1 - d2) > edgeb) { return false; } // 判断是否为边缘跳跃

  return true; // 返回true
}