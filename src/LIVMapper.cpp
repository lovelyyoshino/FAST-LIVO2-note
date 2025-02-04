/* 
This file is part of FAST-LIVO2: Fast, Direct LiDAR-Inertial-Visual Odometry.

Developer: Chunran Zheng <zhengcr@connect.hku.hk>

For commercial use, please contact me at <zhengcr@connect.hku.hk> or
Prof. Fu Zhang at <fuzhang@hku.hk>.

This file is subject to the terms and conditions outlined in the 'LICENSE' file,
which is included as part of this source code package.
*/
#include "LIVMapper.h"

// LIVMapper类的构造函数
LIVMapper::LIVMapper(ros::NodeHandle &nh)
    : extT(0, 0, 0), // 初始化外部平移向量
      extR(M3D::Identity()) // 初始化外部旋转矩阵为单位矩阵
{
  // 初始化外部参数向量
  extrinT.assign(3, 0.0);
  extrinR.assign(9, 0.0);
  cameraextrinT.assign(3, 0.0);
  cameraextrinR.assign(9, 0.0);

  // 创建预处理和IMU处理对象
  p_pre.reset(new Preprocess());
  p_imu.reset(new ImuProcess());

  readParameters(nh); // 读取参数
  VoxelMapConfig voxel_config; // 初始化体素地图配置
  loadVoxelConfig(nh, voxel_config); // 从参数加载体素配置

  // 初始化点云对象
  visual_sub_map.reset(new PointCloudXYZI());
  feats_undistort.reset(new PointCloudXYZI());
  feats_down_body.reset(new PointCloudXYZI());
  feats_down_world.reset(new PointCloudXYZI());
  pcl_w_wait_pub.reset(new PointCloudXYZI());
  pcl_wait_pub.reset(new PointCloudXYZI());
  pcl_wait_save.reset(new PointCloudXYZRGB());

  // 创建体素地图管理器和VIO管理器
  voxelmap_manager.reset(new VoxelMapManager(voxel_config, voxel_map));
  vio_manager.reset(new VIOManager());

  root_dir = ROOT_DIR; // 设置根目录
  initializeFiles(); // 初始化文件
  initializeComponents(); // 初始化组件
  path.header.stamp = ros::Time::now(); // 设置路径时间戳
  path.header.frame_id = "camera_init"; // 设置路径帧ID
}

// LIVMapper类的析构函数
LIVMapper::~LIVMapper() {}

// 读取参数函数
void LIVMapper::readParameters(ros::NodeHandle &nh)
{
  // 读取通用参数
  nh.param<string>("common/lid_topic", lid_topic, "/livox/lidar"); // 激光雷达话题
  nh.param<string>("common/imu_topic", imu_topic, "/livox/imu"); // IMU话题
  nh.param<bool>("common/ros_driver_bug_fix", ros_driver_fix_en, false); // ROS驱动修复标志
  nh.param<int>("common/img_en", img_en, 1); // 图像启用标志
  nh.param<int>("common/lidar_en", lidar_en, 1); // 激光雷达启用标志
  nh.param<string>("common/img_topic", img_topic, "/left_camera/image"); // 图像话题

  // 读取VIO参数
  nh.param<bool>("vio/normal_en", normal_en, true); // 正常模式启用标志
  nh.param<bool>("vio/inverse_composition_en", inverse_composition_en, false); // 反向组合启用标志
  nh.param<int>("vio/max_iterations", max_iterations, 5); // 最大迭代次数
  nh.param<double>("vio/img_point_cov", IMG_POINT_COV, 100); // 图像点协方差
  nh.param<bool>("vio/raycast_en", raycast_en, false); // 射线追踪启用标志
  nh.param<bool>("vio/exposure_estimate_en", exposure_estimate_en, true); // 曝光估计启用标志
  nh.param<double>("vio/inv_expo_cov", inv_expo_cov, 0.2); // 逆曝光协方差
  nh.param<int>("vio/grid_size", grid_size, 5); // 网格大小
  nh.param<int>("vio/grid_n_height", grid_n_height, 17); // 网格高度
  nh.param<int>("vio/patch_pyrimid_level", patch_pyrimid_level, 3); // 补丁金字塔层级
  nh.param<int>("vio/patch_size", patch_size, 8); // 补丁大小
  nh.param<double>("vio/outlier_threshold", outlier_threshold, 1000); // 离群值阈值

  // 读取时间偏移参数
  nh.param<double>("time_offset/exposure_time_init", exposure_time_init, 0.0); // 曝光时间初始化
  nh.param<double>("time_offset/img_time_offset", img_time_offset, 0.0); // 图像时间偏移
  nh.param<double>("time_offset/imu_time_offset", imu_time_offset, 0.0); // IMU时间偏移
  nh.param<bool>("uav/imu_rate_odom", imu_prop_enable, false); // UAV IMU率里程计
  nh.param<bool>("uav/gravity_align_en", gravity_align_en, false); // 重力对齐启用标志

  // 读取EVO参数
  nh.param<string>("evo/seq_name", seq_name, "01"); // 序列名称
  nh.param<bool>("evo/pose_output_en", pose_output_en, false); // 位姿输出启用标志
  nh.param<double>("imu/gyr_cov", gyr_cov, 1.0); // 陀螺仪协方差
  nh.param<double>("imu/acc_cov", acc_cov, 1.0); // 加速度计协方差
  nh.param<int>("imu/imu_int_frame", imu_int_frame, 3); // IMU集成帧数
  nh.param<bool>("imu/imu_en", imu_en, false); // IMU启用标志
  nh.param<bool>("imu/gravity_est_en", gravity_est_en, true); // 重力估计启用标志
  nh.param<bool>("imu/ba_bg_est_en", ba_bg_est_en, true); // 偏置估计启用标志

  // 读取预处理参数
  nh.param<double>("preprocess/blind", p_pre->blind, 0.01); // 盲区参数
  nh.param<double>("preprocess/filter_size_surf", filter_size_surf_min, 0.5); // 表面滤波器大小
  nh.param<int>("preprocess/lidar_type", p_pre->lidar_type, AVIA); // 激光雷达类型
  nh.param<int>("preprocess/scan_line", p_pre->N_SCANS, 6); // 扫描线数
  nh.param<int>("preprocess/point_filter_num", p_pre->point_filter_num, 3); // 点过滤数量
  nh.param<bool>("preprocess/feature_extract_enabled", p_pre->feature_enabled, false); // 特征提取启用标志

  // 读取点云保存参数
  nh.param<int>("pcd_save/interval", pcd_save_interval, -1); // 点云保存间隔
  nh.param<bool>("pcd_save/pcd_save_en", pcd_save_en, false); // 点云保存启用标志
  nh.param<bool>("pcd_save/colmap_output_en", colmap_output_en, false); // COLMAP输出启用标志
  nh.param<double>("pcd_save/filter_size_pcd", filter_size_pcd, 0.5); // 点云保存过滤器大小
  nh.param<vector<double>>("extrin_calib/extrinsic_T", extrinT, vector<double>()); // 外部平移
  nh.param<vector<double>>("extrin_calib/extrinsic_R", extrinR, vector<double>()); // 外部旋转
  nh.param<vector<double>>("extrin_calib/Pcl", cameraextrinT, vector<double>()); // 相机平移
  nh.param<vector<double>>("extrin_calib/Rcl", cameraextrinR, vector<double>()); // 相机旋转
  nh.param<double>("debug/plot_time", plot_time, -10); // 调试绘图时间
  nh.param<int>("debug/frame_cnt", frame_cnt, 6); // 调试帧计数

  // 读取发布参数
  nh.param<double>("publish/blind_rgb_points", blind_rgb_points, 0.01); // RGB点盲区
  nh.param<int>("publish/pub_scan_num", pub_scan_num, 1); // 发布扫描数量
  nh.param<bool>("publish/pub_effect_point_en", pub_effect_point_en, false); // 发布有效点启用标志
  nh.param<bool>("publish/dense_map_en", dense_map_en, false); // 密集地图启用标志

  p_pre->blind_sqr = p_pre->blind * p_pre->blind; // 计算盲区平方
}

// 初始化组件函数
void LIVMapper::initializeComponents() 
{
  downSizeFilterSurf.setLeafSize(filter_size_surf_min, filter_size_surf_min, filter_size_surf_min); // 设置表面下采样滤波器的叶子大小
  extT << VEC_FROM_ARRAY(extrinT); // 将外部平移向量赋值
  extR << MAT_FROM_ARRAY(extrinR); // 将外部旋转矩阵赋值

  voxelmap_manager->extT_ << VEC_FROM_ARRAY(extrinT); // 将外部平移赋值给体素地图管理器
  voxelmap_manager->extR_ << MAT_FROM_ARRAY(extrinR); // 将外部旋转赋值给体素地图管理器

  // 从ROS命名空间加载相机模型
  if (!vk::camera_loader::loadFromRosNs("laserMapping", vio_manager->cam)) throw std::runtime_error("Camera model not correctly specified.");

  // 设置VIO管理器参数
  vio_manager->grid_size = grid_size;
  vio_manager->patch_size = patch_size;
  vio_manager->outlier_threshold = outlier_threshold;
  vio_manager->setImuToLidarExtrinsic(extT, extR); // 设置IMU到激光雷达的外部参数
  vio_manager->setLidarToCameraExtrinsic(cameraextrinR, cameraextrinT); // 设置激光雷达到相机的外部参数
  vio_manager->state = &_state; // 设置状态指针
  vio_manager->state_propagat = &state_propagat; // 设置传播状态指针
  vio_manager->max_iterations = max_iterations; // 设置最大迭代次数
  vio_manager->img_point_cov = IMG_POINT_COV; // 设置图像点协方差
  vio_manager->normal_en = normal_en; // 设置法线启用标志
  vio_manager->inverse_composition_en = inverse_composition_en; // 设置反向组合启用标志
  vio_manager->raycast_en = raycast_en; // 设置射线追踪启用标志
  vio_manager->grid_n_width = grid_n_width; // 设置网格宽度
  vio_manager->grid_n_height = grid_n_height; // 设置网格高度
  vio_manager->patch_pyrimid_level = patch_pyrimid_level; // 设置补丁金字塔层级
  vio_manager->exposure_estimate_en = exposure_estimate_en; // 设置曝光估计启用标志
  vio_manager->colmap_output_en = colmap_output_en; // 设置COLMAP输出启用标志
  vio_manager->initializeVIO(); // 初始化VIO管理器

  // 设置IMU处理参数
  p_imu->set_extrinsic(extT, extR); // 设置IMU外部参数
  p_imu->set_gyr_cov_scale(V3D(gyr_cov, gyr_cov, gyr_cov)); // 设置陀螺仪协方差比例
  p_imu->set_acc_cov_scale(V3D(acc_cov, acc_cov, acc_cov)); // 设置加速度计协方差比例
  p_imu->set_inv_expo_cov(inv_expo_cov); // 设置逆曝光协方差
  p_imu->set_gyr_bias_cov(V3D(0.0001, 0.0001, 0.0001)); // 设置陀螺仪偏置协方差
  p_imu->set_acc_bias_cov(V3D(0.0001, 0.0001, 0.0001)); // 设置加速度计偏置协方差
  p_imu->set_imu_init_frame_num(imu_int_frame); // 设置IMU初始化帧数

  // 根据标志启用或禁用IMU相关功能
  if (!imu_en) p_imu->disable_imu();
  if (!gravity_est_en) p_imu->disable_gravity_est();
  if (!ba_bg_est_en) p_imu->disable_bias_est();
  if (!exposure_estimate_en) p_imu->disable_exposure_est();

  // 设置SLAM模式
  slam_mode_ = (img_en && lidar_en) ? LIVO : imu_en ? ONLY_LIO : ONLY_LO;
}

// 初始化文件函数
void LIVMapper::initializeFiles() 
{
  // 如果启用点云保存和COLMAP输出
  if (pcd_save_en && colmap_output_en)
  {
      const std::string folderPath = std::string(ROOT_DIR) + "/scripts/colmap_output.sh"; // COLMAP输出脚本路径
      
      std::string chmodCommand = "chmod +x " + folderPath; // 设置脚本可执行权限
      
      int chmodRet = system(chmodCommand.c_str());  
      if (chmodRet != 0) {
          std::cerr << "Failed to set execute permissions for the script." << std::endl; // 输出错误信息
          return;
      }

      int executionRet = system(folderPath.c_str()); // 执行脚本
      if (executionRet != 0) {
          std::cerr << "Failed to execute the script." << std::endl; // 输出错误信息
          return;
      }
  }
  // 打开输出文件
  if(colmap_output_en) fout_points.open(std::string(ROOT_DIR) + "Log/Colmap/sparse/0/points3D.txt", std::ios::out);
  if(pcd_save_interval > 0) fout_pcd_pos.open(std::string(ROOT_DIR) + "Log/PCD/scans_pos.json", std::ios::out);
  fout_pre.open(DEBUG_FILE_DIR("mat_pre.txt"), std::ios::out);
  fout_out.open(DEBUG_FILE_DIR("mat_out.txt"), std::ios::out);
}

// 初始化订阅者和发布者函数
void LIVMapper::initializeSubscribersAndPublishers(ros::NodeHandle &nh, image_transport::ImageTransport &it) 
{
  // 根据激光雷达类型选择相应的订阅回调函数
  sub_pcl = p_pre->lidar_type == AVIA ? 
            nh.subscribe(lid_topic, 200000, &LIVMapper::livox_pcl_cbk, this): 
            nh.subscribe(lid_topic, 200000, &LIVMapper::standard_pcl_cbk, this);
  sub_imu = nh.subscribe(imu_topic, 200000, &LIVMapper::imu_cbk, this); // 订阅IMU话题
  sub_img = nh.subscribe(img_topic, 200000, &LIVMapper::img_cbk, this); // 订阅图像话题
  
  // 发布点云和其他消息
  pubLaserCloudFullRes = nh.advertise<sensor_msgs::PointCloud2>("/cloud_registered", 100);
  pubNormal = nh.advertise<visualization_msgs::MarkerArray>("visualization_marker", 100);
  pubSubVisualMap = nh.advertise<sensor_msgs::PointCloud2>("/cloud_visual_sub_map_before", 100);
  pubLaserCloudEffect = nh.advertise<sensor_msgs::PointCloud2>("/cloud_effected", 100);
  pubLaserCloudMap = nh.advertise<sensor_msgs::PointCloud2>("/Laser_map", 100);
  pubOdomAftMapped = nh.advertise<nav_msgs::Odometry>("/aft_mapped_to_init", 10);
  pubPath = nh.advertise<nav_msgs::Path>("/path", 10);
  plane_pub = nh.advertise<visualization_msgs::Marker>("/planner_normal", 1);
  voxel_pub = nh.advertise<visualization_msgs::MarkerArray>("/voxels", 1);
  pubLaserCloudDyn = nh.advertise<sensor_msgs::PointCloud2>("/dyn_obj", 100);
  pubLaserCloudDynRmed = nh.advertise<sensor_msgs::PointCloud2>("/dyn_obj_removed", 100);
  pubLaserCloudDynDbg = nh.advertise<sensor_msgs::PointCloud2>("/dyn_obj_dbg_hist", 100);
  mavros_pose_publisher = nh.advertise<geometry_msgs::PoseStamped>("/mavros/vision_pose/pose", 10);
  pubImage = it.advertise("/rgb_img", 1);
  pubImuPropOdom = nh.advertise<nav_msgs::Odometry>("/LIVO2/imu_propagate", 10000);
  imu_prop_timer = nh.createTimer(ros::Duration(0.004), &LIVMapper::imu_prop_callback, this); // 创建IMU传播定时器
}

// 处理第一帧的函数
void LIVMapper::handleFirstFrame() 
{
  if (!is_first_frame) // 如果不是第一帧
  {
    _first_lidar_time = LidarMeasures.last_lio_update_time; // 获取最后LIO更新时间
    p_imu->first_lidar_time = _first_lidar_time; // 仅用于IMU数据日志
    is_first_frame = true; // 设置为第一帧
    cout << "FIRST LIDAR FRAME!" << endl; // 输出信息
  }
}

// 重力对齐处理函数
void LIVMapper::gravityAlignment() 
{
  if (!p_imu->imu_need_init && !gravity_align_finished) // 如果IMU不需要初始化且重力对齐未完成
  {
    std::cout << "Gravity Alignment Starts" << std::endl; // 输出开始信息
    V3D ez(0, 0, -1), gz(_state.gravity); // 定义重力向量
    Quaterniond G_q_I0 = Quaterniond::FromTwoVectors(gz, ez); // 计算重力对齐的四元数
    M3D G_R_I0 = G_q_I0.toRotationMatrix(); // 计算重力对齐的旋转矩阵

    // 更新状态
    _state.pos_end = G_R_I0 * _state.pos_end;
    _state.rot_end = G_R_I0 * _state.rot_end;
    _state.vel_end = G_R_I0 * _state.vel_end;
    _state.gravity = G_R_I0 * _state.gravity;
    gravity_align_finished = true; // 设置重力对齐完成标志
    std::cout << "Gravity Alignment Finished" << std::endl; // 输出完成信息
  }
}

// 处理IMU数据的函数
void LIVMapper::processImu() 
{
  // double t0 = omp_get_wtime(); // 记录开始时间

  p_imu->Process2(LidarMeasures, _state, feats_undistort); // 处理IMU数据

  if (gravity_align_en) gravityAlignment(); // 如果启用重力对齐，进行重力对齐

  state_propagat = _state; // 更新传播状态
  voxelmap_manager->state_ = _state; // 更新体素地图管理器状态
  voxelmap_manager->feats_undistort_ = feats_undistort; // 更新体素地图管理器去畸变特征
}

// 状态估计和映射函数
void LIVMapper::stateEstimationAndMapping() 
{
  switch (LidarMeasures.lio_vio_flg) // 根据LIO/VIO标志进行不同的处理
  {
    case VIO:
      handleVIO(); // 处理VIO
      break;
    case LIO:
      handleLIO(); // 处理LIO
      break;
  }
}
// 处理VIO（视觉惯性测量）数据的函数
void LIVMapper::handleVIO() 
{
  // 将当前的旋转矩阵转换为欧拉角
  euler_cur = RotMtoEuler(_state.rot_end);
  
  // 将当前状态信息写入文件
  fout_pre << std::setw(20) << LidarMeasures.last_lio_update_time - _first_lidar_time << " " << euler_cur.transpose() * 57.3 << " "
            << _state.pos_end.transpose() << " " << _state.vel_end.transpose() << " " << _state.bias_g.transpose() << " "
            << _state.bias_a.transpose() << " " << V3D(_state.inv_expo_time, 0, 0).transpose() << std::endl;
    
  // 检查待发布的点云是否为空
  if (pcl_w_wait_pub->empty() || (pcl_w_wait_pub == nullptr)) 
  {
    std::cout << "[ VIO ] No point!!!" << std::endl; // 输出信息
    return; // 退出函数
  }
    
  // 输出待发布点云的数量
  std::cout << "[ VIO ] Raw feature num: " << pcl_w_wait_pub->points.size() << std::endl;

  // 根据时间决定是否启用绘图标志
  if (fabs((LidarMeasures.last_lio_update_time - _first_lidar_time) - plot_time) < (frame_cnt / 2 * 0.1)) 
  {
    vio_manager->plot_flag = true; // 启用绘图
  } 
  else 
  {
    vio_manager->plot_flag = false; // 禁用绘图
  }

  // 处理当前帧数据
  vio_manager->processFrame(LidarMeasures.measures.back().img, _pv_list, voxelmap_manager->voxel_map_, LidarMeasures.last_lio_update_time - _first_lidar_time);

  // 如果IMU传播启用
  if (imu_prop_enable) 
  {
    ekf_finish_once = true; // 设置EKF完成标志
    latest_ekf_state = _state; // 更新最新的EKF状态
    latest_ekf_time = LidarMeasures.last_lio_update_time; // 更新时间戳
    state_update_flg = true; // 设置状态更新标志
  }

  // 发布点云和图像
  publish_frame_world(pubLaserCloudFullRes, vio_manager);
  publish_img_rgb(pubImage, vio_manager);

  // 将当前状态信息写入文件
  euler_cur = RotMtoEuler(_state.rot_end);
  fout_out << std::setw(20) << LidarMeasures.last_lio_update_time - _first_lidar_time << " " << euler_cur.transpose() * 57.3 << " "
            << _state.pos_end.transpose() << " " << _state.vel_end.transpose() << " " << _state.bias_g.transpose() << " "
            << _state.bias_a.transpose() << " " << V3D(_state.inv_expo_time, 0, 0).transpose() << " " << feats_undistort->points.size() << std::endl;
}

// 处理LIO（激光惯性测量）数据的函数
void LIVMapper::handleLIO() 
{    
  // 将当前的旋转矩阵转换为欧拉角
  euler_cur = RotMtoEuler(_state.rot_end);
  
  // 将当前状态信息写入文件
  fout_pre << setw(20) << LidarMeasures.last_lio_update_time - _first_lidar_time << " " << euler_cur.transpose() * 57.3 << " "
           << _state.pos_end.transpose() << " " << _state.vel_end.transpose() << " " << _state.bias_g.transpose() << " "
           << _state.bias_a.transpose() << " " << V3D(_state.inv_expo_time, 0, 0).transpose() << endl;
           
  // 检查去畸变特征是否为空
  if (feats_undistort->empty() || (feats_undistort == nullptr)) 
  {
    std::cout << "[ LIO ]: No point!!!" << std::endl; // 输出信息
    return; // 退出函数
  }

  double t0 = omp_get_wtime(); // 记录开始时间

  // 下采样滤波
  downSizeFilterSurf.setInputCloud(feats_undistort);
  downSizeFilterSurf.filter(*feats_down_body);
  
  double t_down = omp_get_wtime(); // 记录下采样时间

  feats_down_size = feats_down_body->points.size(); // 获取下采样后点云的大小
  voxelmap_manager->feats_down_body_ = feats_down_body; // 更新体素地图管理器的下采样特征
  transformLidar(_state.rot_end, _state.pos_end, feats_down_body, feats_down_world); // 转换激光雷达坐标系
  voxelmap_manager->feats_down_world_ = feats_down_world; // 更新体素地图管理器的世界坐标系特征
  voxelmap_manager->feats_down_size_ = feats_down_size; // 更新下采样特征的大小
  
  // 初始化激光雷达地图
  if (!lidar_map_inited) 
  {
    lidar_map_inited = true; // 设置为已初始化
    voxelmap_manager->BuildVoxelMap(); // 构建体素地图
  }

  double t1 = omp_get_wtime(); // 记录状态估计时间

  // 状态估计
  voxelmap_manager->StateEstimation(state_propagat); 
  _state = voxelmap_manager->state_; // 更新状态
  _pv_list = voxelmap_manager->pv_list_; // 更新点云列表

  double t2 = omp_get_wtime(); // 记录状态更新时间

  // 如果IMU传播启用
  if (imu_prop_enable) 
  {
    ekf_finish_once = true; // 设置EKF完成标志
    latest_ekf_state = _state; // 更新最新的EKF状态
    latest_ekf_time = LidarMeasures.last_lio_update_time; // 更新时间戳
    state_update_flg = true; // 设置状态更新标志
  }

  // 如果启用位姿输出
  if (pose_output_en) 
  {
    static bool pos_opend = false; // 位置文件打开标志
    static int ocount = 0; // 计数器
    std::ofstream outFile, evoFile; // 输出文件流
    if (!pos_opend) 
    {
      // 打开新文件
      evoFile.open(std::string(ROOT_DIR) + "Log/result/" + seq_name + ".txt", std::ios::out);
      pos_opend = true; // 设置为已打开
      if (!evoFile.is_open()) ROS_ERROR("open fail\n"); // 输出错误信息
    } 
    else 
    {
      // 以追加方式打开文件
      evoFile.open(std::string(ROOT_DIR) + "Log/result/" + seq_name + ".txt", std::ios::app);
      if (!evoFile.is_open()) ROS_ERROR("open fail\n"); // 输出错误信息
    }
    Eigen::Matrix4d outT; // 4x4变换矩阵
    Eigen::Quaterniond q(_state.rot_end); // 将旋转状态转换为四元数
    // 将状态信息写入文件
    evoFile << std::fixed;
    evoFile << LidarMeasures.last_lio_update_time << " " << _state.pos_end[0] << " " << _state.pos_end[1] << " " << _state.pos_end[2] << " "
            << q.x() << " " << q.y() << " " << q.z() << " " << q.w() << std::endl;
  }
  
  // 将当前的旋转矩阵转换为欧拉角并生成四元数
  euler_cur = RotMtoEuler(_state.rot_end);
  geoQuat = tf::createQuaternionMsgFromRollPitchYaw(euler_cur(0), euler_cur(1), euler_cur(2));
  publish_odometry(pubOdomAftMapped); // 发布里程计数据

  double t3 = omp_get_wtime(); // 记录发布里程计时间

  // 处理激光雷达点云
  PointCloudXYZI::Ptr world_lidar(new PointCloudXYZI());
  transformLidar(_state.rot_end, _state.pos_end, feats_down_body, world_lidar); // 转换激光雷达点云
  for (size_t i = 0; i < world_lidar->points.size(); i++) 
  {
    voxelmap_manager->pv_list_[i].point_w << world_lidar->points[i].x, world_lidar->points[i].y, world_lidar->points[i].z; // 更新点云位置
    M3D point_crossmat = voxelmap_manager->cross_mat_list_[i]; // 获取交叉矩阵
    M3D var = voxelmap_manager->body_cov_list_[i]; // 获取点云协方差
    // 更新协方差
    var = (_state.rot_end * extR) * var * (_state.rot_end * extR).transpose() +
          (-point_crossmat) * _state.cov.block<3, 3>(0, 0) * (-point_crossmat).transpose() + _state.cov.block<3, 3>(3, 3);
    voxelmap_manager->pv_list_[i].var = var; // 更新点云变量
  }
  voxelmap_manager->UpdateVoxelMap(voxelmap_manager->pv_list_); // 更新体素地图
  std::cout << "[ LIO ] Update Voxel Map" << std::endl; // 输出更新信息
  _pv_list = voxelmap_manager->pv_list_; // 更新点云列表
  
  double t4 = omp_get_wtime(); // 记录更新后时间

  // 如果启用滑动地图
  if(voxelmap_manager->config_setting_.map_sliding_en)
  {
    if(voxelmap_manager->mapSliding()) 
    {
      // update_local_voxel_map();
    }
    // publish_local_voxelmap(local_voxel_clouds_publisher);
  }
  
  // 处理激光雷达点云
  PointCloudXYZI::Ptr laserCloudFullRes(dense_map_en ? feats_undistort : feats_down_body);
  int size = laserCloudFullRes->points.size(); // 获取点云大小
  PointCloudXYZI::Ptr laserCloudWorld(new PointCloudXYZI(size, 1)); // 创建新的点云对象

  // 转换点云到世界坐标系
  for (int i = 0; i < size; i++) 
  {
    RGBpointBodyToWorld(&laserCloudFullRes->points[i], &laserCloudWorld->points[i]);
  }
  *pcl_w_wait_pub = *laserCloudWorld; // 更新待发布的点云

  // 根据启用标志发布点云
  if (!img_en) publish_frame_world(pubLaserCloudFullRes, vio_manager);
  if (pub_effect_point_en) publish_effect_world(pubLaserCloudEffect, voxelmap_manager->ptpl_list_);
  publish_path(pubPath); // 发布路径
  publish_mavros(mavros_pose_publisher); // 发布MAVROS数据

  frame_num++; // 增加帧计数
  aver_time_consu = aver_time_consu * (frame_num - 1) / frame_num + (t4 - t0) / frame_num; // 更新平均时间

  // 输出时间统计信息
  printf("\033[1;34m+-------------------------------------------------------------+\033[0m\n");
  printf("\033[1;34m|                         LIO Mapping Time                    |\033[0m\n");
  printf("\033[1;34m+-------------------------------------------------------------+\033[0m\n");
  printf("\033[1;34m| %-29s | %-27s |\033[0m\n", "Algorithm Stage", "Time (secs)");
  printf("\033[1;34m+-------------------------------------------------------------+\033[0m\n");
  printf("\033[1;36m| %-29s | %-27f |\033[0m\n", "DownSample", t_down - t0);
  printf("\033[1;36m| %-29s | %-27f |\033[0m\n", "ICP", t2 - t1);
  printf("\033[1;36m| %-29s | %-27f |\033[0m\n", "updateVoxelMap", t4 - t3);
  printf("\033[1;34m+-------------------------------------------------------------+\033[0m\n");
  printf("\033[1;36m| %-29s | %-27f |\033[0m\n", "Current Total Time", t4 - t0);
  printf("\033[1;36m| %-29s | %-27f |\033[0m\n", "Average Total Time", aver_time_consu);
  printf("\033[1;34m+-------------------------------------------------------------+\033[0m\n");

  // 将当前状态信息写入文件
  euler_cur = RotMtoEuler(_state.rot_end);
  fout_out << std::setw(20) << LidarMeasures.last_lio_update_time - _first_lidar_time << " " << euler_cur.transpose() * 57.3 << " "
            << _state.pos_end.transpose() << " " << _state.vel_end.transpose() << " " << _state.bias_g.transpose() << " "
            << _state.bias_a.transpose() << " " << V3D(_state.inv_expo_time, 0, 0).transpose() << " " << feats_undistort->points.size() << std::endl;
}

// 保存PCD文件的函数
void LIVMapper::savePCD() 
{
  // 如果启用PCD保存且待保存点云数量大于0且保存间隔小于0
  if (pcd_save_en && pcl_wait_save->points.size() > 0 && pcd_save_interval < 0) 
  {
    pcl::PointCloud<pcl::PointXYZRGB>::Ptr downsampled_cloud(new pcl::PointCloud<pcl::PointXYZRGB>); // 创建下采样点云
    pcl::VoxelGrid<pcl::PointXYZRGB> voxel_filter; // 创建体素滤波器
    voxel_filter.setInputCloud(pcl_wait_save); // 设置输入点云
    voxel_filter.setLeafSize(filter_size_pcd, filter_size_pcd, filter_size_pcd); // 设置滤波器的叶子大小
    voxel_filter.filter(*downsampled_cloud); // 执行滤波

    std::string raw_points_dir = std::string(ROOT_DIR) + "Log/PCD/all_raw_points.pcd"; // 原始点云路径
    std::string downsampled_points_dir = std::string(ROOT_DIR) + "Log/PCD/all_downsampled_points.pcd"; // 下采样点云路径

    pcl::PCDWriter pcd_writer; // 创建PCD写入器

    // 保存原始点云数据
    pcd_writer.writeBinary(raw_points_dir, *pcl_wait_save);
    std::cout << GREEN << "Raw point cloud data saved to: " << raw_points_dir 
              << " with point count: " << pcl_wait_save->points.size() << RESET << std::endl;

    // 保存下采样点云数据
    pcd_writer.writeBinary(downsampled_points_dir, *downsampled_cloud);
    std::cout << GREEN << "Downsampled point cloud data saved to: " << downsampled_points_dir 
          << " with point count after filtering: " << downsampled_cloud->points.size() << RESET << std::endl;

    // 如果启用COLMAP输出
    if(colmap_output_en)
    {
      fout_points << "# 3D point list with one line of data per point\n";
      fout_points << "#  POINT_ID, X, Y, Z, R, G, B, ERROR\n";
      for (size_t i = 0; i < downsampled_cloud->size(); ++i) 
      {
          const auto& point = downsampled_cloud->points[i];
          fout_points << i << " "
                      << std::fixed << std::setprecision(6)
                      << point.x << " " << point.y << " " << point.z << " "
                      << static_cast<int>(point.r) << " "
                      << static_cast<int>(point.g) << " "
                      << static_cast<int>(point.b) << " "
                      << 0 << std::endl; // 输出点的ID和坐标
      }
    }
  }
}

// LIVMapper类的主运行函数
void LIVMapper::run() 
{
  ros::Rate rate(5000); // 设置循环频率为5000Hz
  while (ros::ok()) // 当ROS处于运行状态时
  {
    ros::spinOnce(); // 处理ROS回调
    if (!sync_packages(LidarMeasures)) // 同步数据包
    {
      rate.sleep(); // 休眠以保持循环频率
      continue; // 继续下一次循环
    }   
    handleFirstFrame(); // 处理第一帧数据

    processImu(); // 处理IMU数据

    // if (!p_imu->imu_time_init) continue; // 如果IMU时间未初始化则跳过

    stateEstimationAndMapping(); // 状态估计与映射
  }
  savePCD(); // 保存点云数据
}

// 单次IMU传播的函数
void LIVMapper::prop_imu_once(StatesGroup &imu_prop_state, const double dt, V3D acc_avr, V3D angvel_avr)
{
  double mean_acc_norm = p_imu->IMU_mean_acc_norm; // 获取IMU的平均加速度规范
  acc_avr = acc_avr * G_m_s2 / mean_acc_norm - imu_prop_state.bias_a; // 计算加速度
  angvel_avr -= imu_prop_state.bias_g; // 减去角速度偏差

  M3D Exp_f = Exp(angvel_avr, dt); // 计算旋转的指数映射
  /* IMU姿态传播 */
  imu_prop_state.rot_end = imu_prop_state.rot_end * Exp_f; // 更新IMU的旋转状态

  /* IMU的全局特定加速度 */
  V3D acc_imu = imu_prop_state.rot_end * acc_avr + V3D(imu_prop_state.gravity[0], imu_prop_state.gravity[1], imu_prop_state.gravity[2]);

  /* IMU传播位置 */
  imu_prop_state.pos_end = imu_prop_state.pos_end + imu_prop_state.vel_end * dt + 0.5 * acc_imu * dt * dt; // 更新位置

  /* IMU的速度 */
  imu_prop_state.vel_end = imu_prop_state.vel_end + acc_imu * dt; // 更新速度
}

// IMU传播回调函数
void LIVMapper::imu_prop_callback(const ros::TimerEvent &e)
{
  if (p_imu->imu_need_init || !new_imu || !ekf_finish_once) { return; } // 如果IMU需要初始化或没有新IMU数据则返回
  mtx_buffer_imu_prop.lock(); // 锁定IMU传播缓冲区
  new_imu = false; // 控制传播频率与IMU频率一致
  if (imu_prop_enable && !prop_imu_buffer.empty())
  {
    static double last_t_from_lidar_end_time = 0; // 保存上一次激光雷达结束时间
    if (state_update_flg) // 如果状态更新标志为真
    {
      imu_propagate = latest_ekf_state; // 设置当前IMU传播状态为最新EKF状态
      // 丢弃所有无用的IMU包
      while ((!prop_imu_buffer.empty() && prop_imu_buffer.front().header.stamp.toSec() < latest_ekf_time))
      {
        prop_imu_buffer.pop_front(); // 从缓冲区中移除过时的IMU数据
      }
      last_t_from_lidar_end_time = 0; // 重置激光雷达结束时间
      for (int i = 0; i < prop_imu_buffer.size(); i++)
      {
        double t_from_lidar_end_time = prop_imu_buffer[i].header.stamp.toSec() - latest_ekf_time; // 计算IMU时间与激光雷达结束时间差
        double dt = t_from_lidar_end_time - last_t_from_lidar_end_time; // 计算时间间隔
        // cout << "prop dt" << dt << ", " << t_from_lidar_end_time << ", " << last_t_from_lidar_end_time << endl;
        V3D acc_imu(prop_imu_buffer[i].linear_acceleration.x, prop_imu_buffer[i].linear_acceleration.y, prop_imu_buffer[i].linear_acceleration.z); // 获取IMU线性加速度
        V3D omg_imu(prop_imu_buffer[i].angular_velocity.x, prop_imu_buffer[i].angular_velocity.y, prop_imu_buffer[i].angular_velocity.z); // 获取IMU角速度
        prop_imu_once(imu_propagate, dt, acc_imu, omg_imu); // 调用IMU传播函数
        last_t_from_lidar_end_time = t_from_lidar_end_time; // 更新激光雷达结束时间
      }
      state_update_flg = false; // 重置状态更新标志
    }
    else // 如果状态更新标志为假
    {
      V3D acc_imu(newest_imu.linear_acceleration.x, newest_imu.linear_acceleration.y, newest_imu.linear_acceleration.z); // 获取最新IMU的线性加速度
      V3D omg_imu(newest_imu.angular_velocity.x, newest_imu.angular_velocity.y, newest_imu.angular_velocity.z); // 获取最新IMU的角速度
      double t_from_lidar_end_time = newest_imu.header.stamp.toSec() - latest_ekf_time; // 计算IMU时间与激光雷达结束时间差
      double dt = t_from_lidar_end_time - last_t_from_lidar_end_time; // 计算时间间隔
      prop_imu_once(imu_propagate, dt, acc_imu, omg_imu); // 调用IMU传播函数
      last_t_from_lidar_end_time = t_from_lidar_end_time; // 更新激光雷达结束时间
    }

    V3D posi, vel_i; // 声明位置和速度变量
    Eigen::Quaterniond q; // 声明四元数变量
    posi = imu_propagate.pos_end; // 获取IMU传播后的末位置
    vel_i = imu_propagate.vel_end; // 获取IMU传播后的末速度
    q = Eigen::Quaterniond(imu_propagate.rot_end); // 获取IMU传播后的旋转状态
    imu_prop_odom.header.frame_id = "world"; // 设置坐标系为世界坐标系
    imu_prop_odom.header.stamp = newest_imu.header.stamp; // 设置时间戳
    imu_prop_odom.pose.pose.position.x = posi.x(); // 设置位置x坐标
    imu_prop_odom.pose.pose.position.y = posi.y(); // 设置位置y坐标
    imu_prop_odom.pose.pose.position.z = posi.z(); // 设置位置z坐标
    imu_prop_odom.pose.pose.orientation.w = q.w(); // 设置四元数w分量
    imu_prop_odom.pose.pose.orientation.x = q.x(); // 设置四元数x分量
    imu_prop_odom.pose.pose.orientation.y = q.y(); // 设置四元数y分量
    imu_prop_odom.pose.pose.orientation.z = q.z(); // 设置四元数z分量
    imu_prop_odom.twist.twist.linear.x = vel_i.x(); // 设置速度x分量
    imu_prop_odom.twist.twist.linear.y = vel_i.y(); // 设置速度y分量
    imu_prop_odom.twist.twist.linear.z = vel_i.z(); // 设置速度z分量
    pubImuPropOdom.publish(imu_prop_odom); // 发布IMU传播的里程计信息
  }
  mtx_buffer_imu_prop.unlock(); // 解锁IMU传播缓冲区
}

// 激光雷达坐标系转换函数
void LIVMapper::transformLidar(const Eigen::Matrix3d rot, const Eigen::Vector3d t, const PointCloudXYZI::Ptr &input_cloud, PointCloudXYZI::Ptr &trans_cloud)
{
  PointCloudXYZI().swap(*trans_cloud); // 交换点云内容
  trans_cloud->reserve(input_cloud->size()); // 预留空间
  for (size_t i = 0; i < input_cloud->size(); i++) // 遍历输入点云
  {
    pcl::PointXYZINormal p_c = input_cloud->points[i]; // 获取当前点
    Eigen::Vector3d p(p_c.x, p_c.y, p_c.z); // 创建点的3D向量
    p = (rot * (extR * p + extT) + t); // 进行坐标变换
    PointType pi; // 声明输出点
    pi.x = p(0); // 设置输出点x坐标
    pi.y = p(1); // 设置输出点y坐标
    pi.z = p(2); // 设置输出点z坐标
    pi.intensity = p_c.intensity; // 设置输出点强度
    trans_cloud->points.push_back(pi); // 将输出点添加到点云中
  }
}

// 从体坐标系转换到世界坐标系的函数
void LIVMapper::pointBodyToWorld(const PointType &pi, PointType &po)
{
  V3D p_body(pi.x, pi.y, pi.z); // 创建体坐标系点
  V3D p_global(_state.rot_end * (extR * p_body + extT) + _state.pos_end); // 转换到全局坐标系
  po.x = p_global(0); // 设置输出点x坐标
  po.y = p_global(1); // 设置输出点y坐标
  po.z = p_global(2); // 设置输出点z坐标
  po.intensity = pi.intensity; // 设置输出点强度
}

// 从体坐标系转换到世界坐标系的模板函数
template <typename T> void LIVMapper::pointBodyToWorld(const Matrix<T, 3, 1> &pi, Matrix<T, 3, 1> &po)
{
  V3D p_body(pi[0], pi[1], pi[2]); // 创建体坐标系点
  V3D p_global(_state.rot_end * (extR * p_body + extT) + _state.pos_end); // 转换到全局坐标系
  po[0] = p_global(0); // 设置输出点x坐标
  po[1] = p_global(1); // 设置输出点y坐标
  po[2] = p_global(2); // 设置输出点z坐标
}

// 从体坐标系转换到世界坐标系并返回的模板函数
template <typename T> Matrix<T, 3, 1> LIVMapper::pointBodyToWorld(const Matrix<T, 3, 1> &pi)
{
  V3D p(pi[0], pi[1], pi[2]); // 创建体坐标系点
  p = (_state.rot_end * (extR * p + extT) + _state.pos_end); // 转换到全局坐标系
  Matrix<T, 3, 1> po(p[0], p[1], p[2]); // 创建输出点
  return po; // 返回输出点
}

// 从体坐标系转换到世界坐标系的RGB点转换函数
void LIVMapper::RGBpointBodyToWorld(PointType const *const pi, PointType *const po)
{
  V3D p_body(pi->x, pi->y, pi->z); // 创建体坐标系点
  V3D p_global(_state.rot_end * (extR * p_body + extT) + _state.pos_end); // 转换到全局坐标系
  po->x = p_global(0); // 设置输出点x坐标
  po->y = p_global(1); // 设置输出点y坐标
  po->z = p_global(2); // 设置输出点z坐标
  po->intensity = pi->intensity; // 设置输出点强度
}

// 标准点云回调函数
void LIVMapper::standard_pcl_cbk(const sensor_msgs::PointCloud2::ConstPtr &msg)
{
  if (!lidar_en) return; // 如果激光雷达未启用则返回
  mtx_buffer.lock(); // 锁定缓冲区
  // cout<<"got feature"<<endl;
  if (msg->header.stamp.toSec() < last_timestamp_lidar) // 检查时间戳是否回绕
  {
    ROS_ERROR("lidar loop back, clear buffer"); // 输出错误信息
    lid_raw_data_buffer.clear(); // 清空缓冲区
  }
  // ROS_INFO("get point cloud at time: %.6f", msg->header.stamp.toSec());
  PointCloudXYZI::Ptr ptr(new PointCloudXYZI()); // 创建新的点云
  p_pre->process(msg, ptr); // 处理点云数据
  lid_raw_data_buffer.push_back(ptr); // 将处理后的点云添加到缓冲区
  lid_header_time_buffer.push_back(msg->header.stamp.toSec()); // 添加时间戳
  last_timestamp_lidar = msg->header.stamp.toSec(); // 更新最后的激光雷达时间戳

  mtx_buffer.unlock(); // 解锁缓冲区
  sig_buffer.notify_all(); // 通知所有等待的线程
}

// Livox激光雷达回调函数
void LIVMapper::livox_pcl_cbk(const livox_ros_driver::CustomMsg::ConstPtr &msg_in)
{
  if (!lidar_en) return; // 如果激光雷达未启用则返回
  mtx_buffer.lock(); // 锁定缓冲区
  livox_ros_driver::CustomMsg::Ptr msg(new livox_ros_driver::CustomMsg(*msg_in)); // 创建新的Livox消息
  // if ((abs(msg->header.stamp.toSec() - last_timestamp_lidar) > 0.2 && last_timestamp_lidar > 0) || sync_jump_flag)
  // {
  //   ROS_WARN("lidar jumps %.3f\n", msg->header.stamp.toSec() - last_timestamp_lidar);
  //   sync_jump_flag = true;
  //   msg->header.stamp = ros::Time().fromSec(last_timestamp_lidar + 0.1);
  // }
  if (abs(last_timestamp_imu - msg->header.stamp.toSec()) > 1.0 && !imu_buffer.empty()) // 检查IMU与激光雷达时间戳差异
  {
    double timediff_imu_wrt_lidar = last_timestamp_imu - msg->header.stamp.toSec(); // 计算IMU时间与激光雷达时间的差
    printf("\033[95mSelf sync IMU and LiDAR, HARD time lag is %.10lf \n\033[0m", timediff_imu_wrt_lidar - 0.100); // 输出时间延迟
    // imu_time_offset = timediff_imu_wrt_lidar;
  }

  double cur_head_time = msg->header.stamp.toSec(); // 获取当前头部时间
  ROS_INFO("Get LiDAR, its header time: %.6f", cur_head_time); // 输出激光雷达时间
  if (cur_head_time < last_timestamp_lidar) // 检查时间戳是否回绕
  {
    ROS_ERROR("lidar loop back, clear buffer"); // 输出错误信息
    lid_raw_data_buffer.clear(); // 清空缓冲区
  }
  // ROS_INFO("get point cloud at time: %.6f", msg->header.stamp.toSec());
  PointCloudXYZI::Ptr ptr(new PointCloudXYZI()); // 创建新的点云
  p_pre->process(msg, ptr); // 处理点云数据
  lid_raw_data_buffer.push_back(ptr); // 将处理后的点云添加到缓冲区
  lid_header_time_buffer.push_back(cur_head_time); // 添加时间戳
  last_timestamp_lidar = cur_head_time; // 更新最后的激光雷达时间戳

  mtx_buffer.unlock(); // 解锁缓冲区
  sig_buffer.notify_all(); // 通知所有等待的线程
}

// IMU回调函数
void LIVMapper::imu_cbk(const sensor_msgs::Imu::ConstPtr &msg_in)
{
  if (!imu_en) return; // 如果IMU未启用则返回

  if (last_timestamp_lidar < 0.0) return; // 如果激光雷达时间戳无效则返回
  // ROS_INFO("get imu at time: %.6f", msg_in->header.stamp.toSec());
  sensor_msgs::Imu::Ptr msg(new sensor_msgs::Imu(*msg_in)); // 创建新的IMU消息
  msg->header.stamp = ros::Time().fromSec(msg->header.stamp.toSec() - imu_time_offset); // 校正时间戳
  double timestamp = msg->header.stamp.toSec(); // 获取IMU时间戳

  if (fabs(last_timestamp_lidar - timestamp) > 0.5 && (!ros_driver_fix_en)) // 检查IMU与激光雷达时间戳差异
  {
    ROS_WARN("IMU and LiDAR not synced! delta time: %lf .\n", last_timestamp_lidar - timestamp); // 输出警告信息
  }

  if (ros_driver_fix_en) timestamp += std::round(last_timestamp_lidar - timestamp); // 修正时间戳
  msg->header.stamp = ros::Time().fromSec(timestamp); // 设置校正后的时间戳

  mtx_buffer.lock(); // 锁定缓冲区

  if (last_timestamp_imu > 0.0 && timestamp < last_timestamp_imu) // 检查IMU时间戳是否回绕
  {
    mtx_buffer.unlock(); // 解锁缓冲区
    sig_buffer.notify_all(); // 通知所有等待的线程
    ROS_ERROR("imu loop back, offset: %lf \n", last_timestamp_imu - timestamp); // 输出错误信息
    return; // 返回
  }

  if (last_timestamp_imu > 0.0 && timestamp > last_timestamp_imu + 0.2) // 检查IMU时间戳是否跳变
  {
    ROS_WARN("imu time stamp Jumps %0.4lf seconds \n", timestamp - last_timestamp_imu); // 输出警告信息
    mtx_buffer.unlock(); // 解锁缓冲区
    sig_buffer.notify_all(); // 通知所有等待的线程
    return; // 返回
  }

  last_timestamp_imu = timestamp; // 更新最后的IMU时间戳

  imu_buffer.push_back(msg); // 将IMU消息添加到缓冲区
  // cout<<"got imu: "<<timestamp<<" imu size "<<imu_buffer.size()<<endl;
  mtx_buffer.unlock(); // 解锁缓冲区
  if (imu_prop_enable) // 如果IMU传播启用
  {
    mtx_buffer_imu_prop.lock(); // 锁定IMU传播缓冲区
    if (imu_prop_enable && !p_imu->imu_need_init) { prop_imu_buffer.push_back(*msg); } // 将IMU消息添加到传播缓冲区
    newest_imu = *msg; // 更新最新IMU消息
    new_imu = true; // 设置新IMU消息标志
    mtx_buffer_imu_prop.unlock(); // 解锁IMU传播缓冲区
  }
  sig_buffer.notify_all(); // 通知所有等待的线程
}

// 从ROS消息获取图像的函数
cv::Mat LIVMapper::getImageFromMsg(const sensor_msgs::ImageConstPtr &img_msg)
{
  cv::Mat img; // 声明图像变量
  img = cv_bridge::toCvCopy(img_msg, "bgr8")->image; // 将ROS图像消息转换为OpenCV格式
  return img; // 返回图像
}

// 图像回调函数
void LIVMapper::img_cbk(const sensor_msgs::ImageConstPtr &msg_in)
{
  if (!img_en) return; // 如果图像未启用则返回
  sensor_msgs::Image::Ptr msg(new sensor_msgs::Image(*msg_in)); // 创建新的图像消息
  // if ((abs(msg->header.stamp.toSec() - last_timestamp_img) > 0.2 && last_timestamp_img > 0) || sync_jump_flag)
  // {
  //   ROS_WARN("img jumps %.3f\n", msg->header.stamp.toSec() - last_timestamp_img);
  //   sync_jump_flag = true;
  //   msg->header.stamp = ros::Time().fromSec(last_timestamp_img + 0.1);
  // }

  // Hiliti2022 40Hz
  // if (hilti_en)
  // {
  //   i++;
  //   if (i % 4 != 0) return;
  // }
  // double msg_header_time =  msg->header.stamp.toSec();
  double msg_header_time = msg->header.stamp.toSec() + img_time_offset; // 获取图像时间戳并校正
  if (abs(msg_header_time - last_timestamp_img) < 0.001) return; // 如果时间戳几乎相同则返回
  ROS_INFO("Get image, its header time: %.6f", msg_header_time); // 输出图像时间
  if (last_timestamp_lidar < 0) return; // 如果激光雷达时间戳无效则返回

  if (msg_header_time < last_timestamp_img) // 检查时间戳是否回绕
  {
    ROS_ERROR("image loop back. \n"); // 输出错误信息
    return; // 返回
  }

  mtx_buffer.lock(); // 锁定缓冲区

  double img_time_correct = msg_header_time; // 校正图像时间
  // last_timestamp_lidar + 0.105;

  if (img_time_correct - last_timestamp_img < 0.02) // 检查图像时间差
  {
    ROS_WARN("Image need Jumps: %.6f", img_time_correct); // 输出警告信息
    mtx_buffer.unlock(); // 解锁缓冲区
    sig_buffer.notify_all(); // 通知所有等待的线程
    return; // 返回
  }

  cv::Mat img_cur = getImageFromMsg(msg); // 获取当前图像
  img_buffer.push_back(img_cur); // 将图像添加到缓冲区
  img_time_buffer.push_back(img_time_correct); // 添加图像时间戳

  // ROS_INFO("Correct Image time: %.6f", img_time_correct);

  last_timestamp_img = img_time_correct; // 更新最后的图像时间戳
  // cv::imshow("img", img);
  // cv::waitKey(1);
  // cout<<"last_timestamp_img:::"<<last_timestamp_img<<endl;
  mtx_buffer.unlock(); // 解锁缓冲区
  sig_buffer.notify_all(); // 通知所有等待的线程
}

// 同步数据包的函数
bool LIVMapper::sync_packages(LidarMeasureGroup &meas)
{
  if (lid_raw_data_buffer.empty() && lidar_en) return false; // 如果激光雷达缓冲区为空则返回假
  if (img_buffer.empty() && img_en) return false; // 如果图像缓冲区为空则返回假
  if (imu_buffer.empty() && imu_en) return false; // 如果IMU缓冲区为空则返回假

  switch (slam_mode_) // 根据SLAM模式进行处理
  {
  case ONLY_LIO: // 仅LIO模式
  {
    if (meas.last_lio_update_time < 0.0) meas.last_lio_update_time = lid_header_time_buffer.front(); // 更新LIO时间
    if (!lidar_pushed) // 如果未推送激光雷达数据
    {
      // 如果未推送激光雷达到测量数据缓冲区
      meas.lidar = lid_raw_data_buffer.front(); // 推送第一个激光雷达主题
      if (meas.lidar->points.size() <= 1) return false; // 如果点云数量小于等于1则返回假

      meas.lidar_frame_beg_time = lid_header_time_buffer.front(); // 生成激光雷达帧开始时间
      meas.lidar_frame_end_time = meas.lidar_frame_beg_time + meas.lidar->points.back().curvature / double(1000); // 计算激光雷达扫描结束时间
      meas.pcl_proc_cur = meas.lidar; // 当前处理点云
      lidar_pushed = true; // 设置标志
    }

    if (imu_en && last_timestamp_imu < meas.lidar_frame_end_time) // 检查IMU时间戳
    { // 等待IMU消息需要大于激光雷达帧结束时间
      // 确保完成传播
      // ROS_ERROR("out sync");
      return false; // 返回假
    }

    struct MeasureGroup m; // 标准方法保存IMU消息

    m.imu.clear(); // 清空IMU数据
    m.lio_time = meas.lidar_frame_end_time; // 设置LIO时间
    mtx_buffer.lock(); // 锁定缓冲区
    while (!imu_buffer.empty()) // 遍历IMU缓冲区
    {
      if (imu_buffer.front()->header.stamp.toSec() > meas.lidar_frame_end_time) break; // 如果IMU时间戳大于激光雷达结束时间则退出
      m.imu.push_back(imu_buffer.front()); // 将IMU数据添加到测量组
      imu_buffer.pop_front(); // 从缓冲区中移除IMU数据
    }
    lid_raw_data_buffer.pop_front(); // 从激光雷达缓冲区中移除第一个数据
    lid_header_time_buffer.pop_front(); // 从时间戳缓冲区中移除第一个数据
    mtx_buffer.unlock(); // 解锁缓冲区
    sig_buffer.notify_all(); // 通知所有等待的线程

    meas.lio_vio_flg = LIO; // 处理激光雷达主题，因此时间戳应为激光雷达扫描结束
    meas.measures.push_back(m); // 将测量组添加到测量数据中
    // ROS_INFO("ONlY HAS LiDAR and IMU, NO IMAGE!");
    lidar_pushed = false; // 同步一个完整的激光雷达扫描
    return true; // 返回真

    break;
  }

  case LIVO: // LIVO模式
  {
    /*** 对于LIVO模式，LIO更新的时间设置为与VIO相同，LIO
     * 首先在VIO之前更新 ***/
    EKF_STATE last_lio_vio_flg = meas.lio_vio_flg; // 获取最后的LIO-VIO标志
    // double t0 = omp_get_wtime();
    switch (last_lio_vio_flg)
    {
    // double img_capture_time = meas.lidar_frame_beg_time + exposure_time_init;
    case WAIT:
    case VIO:
    {
      // printf("!!! meas.lio_vio_flg: %d \n", meas.lio_vio_flg);
      double img_capture_time = img_time_buffer.front() + exposure_time_init; // 获取图像捕获时间
      /*** 有图像主题，但图像主题时间戳大于激光雷达结束时间，
       * 处理激光雷达主题。更新LIO后，meas.lidar_frame_end_time
       * 将被刷新。 ***/
      if (meas.last_lio_update_time < 0.0) meas.last_lio_update_time = lid_header_time_buffer.front(); // 更新LIO时间

      double lid_newest_time = lid_header_time_buffer.back() + lid_raw_data_buffer.back()->points.back().curvature / double(1000); // 获取最新激光雷达时间
      double imu_newest_time = imu_buffer.back()->header.stamp.toSec(); // 获取最新IMU时间

      if (img_capture_time < meas.last_lio_update_time + 0.00001) // 如果图像捕获时间小于LIO更新的时间
      {
        img_buffer.pop_front(); // 从图像缓冲区中移除第一个数据
        img_time_buffer.pop_front(); // 从时间戳缓冲区中移除第一个数据
        ROS_ERROR("[ Data Cut ] Throw one image frame! \n"); // 输出错误信息
        return false; // 返回假
      }

      if (img_capture_time > lid_newest_time || img_capture_time > imu_newest_time) // 检查时间戳是否超出范围
      {
        // ROS_ERROR("lost first camera frame");
        // printf("img_capture_time, lid_newest_time, imu_newest_time: %lf , %lf
        // , %lf \n", img_capture_time, lid_newest_time, imu_newest_time);
        return false; // 返回假
      }

      struct MeasureGroup m; // 创建测量组

      // printf("[ Data Cut ] LIO \n");
      // printf("[ Data Cut ] img_capture_time: %lf \n", img_capture_time);
      m.imu.clear(); // 清空IMU数据
      m.lio_time = img_capture_time; // 设置LIO时间
      mtx_buffer.lock(); // 锁定缓冲区
      while (!imu_buffer.empty()) // 遍历IMU缓冲区
      {
        if (imu_buffer.front()->header.stamp.toSec() > m.lio_time) break; // 如果IMU时间戳大于LIO时间则退出

        if (imu_buffer.front()->header.stamp.toSec() > meas.last_lio_update_time) m.imu.push_back(imu_buffer.front()); // 如果IMU时间戳大于LIO更新时间则添加到测量组

        imu_buffer.pop_front(); // 从IMU缓冲区中移除IMU数据
        // printf("[ Data Cut ] imu time: %lf \n",
        // imu_buffer.front()->header.stamp.toSec());
      }
      mtx_buffer.unlock(); // 解锁缓冲区
      sig_buffer.notify_all(); // 通知所有等待的线程

      *(meas.pcl_proc_cur) = *(meas.pcl_proc_next); // 交换当前和下一点云
      PointCloudXYZI().swap(*meas.pcl_proc_next); // 清空下一点云

      int lid_frame_num = lid_raw_data_buffer.size(); // 获取激光雷达帧数量
      int max_size = meas.pcl_proc_cur->size() + 24000 * lid_frame_num; // 计算最大大小
      meas.pcl_proc_cur->reserve(max_size); // 预留空间
      meas.pcl_proc_next->reserve(max_size); // 预留空间

      while (!lid_raw_data_buffer.empty()) // 遍历激光雷达缓冲区
      {
        if (lid_header_time_buffer.front() > img_capture_time) break; // 如果激光雷达时间戳大于图像捕获时间则退出
        auto pcl(lid_raw_data_buffer.front()->points); // 获取激光雷达点云
        double frame_header_time(lid_header_time_buffer.front()); // 获取激光雷达帧时间
        float max_offs_time_ms = (m.lio_time - frame_header_time) * 1000.0f; // 计算最大偏移时间

        for (int i = 0; i < pcl.size(); i++) // 遍历点云
        {
          auto pt = pcl[i]; // 获取当前点
          if (pcl[i].curvature < max_offs_time_ms) // 如果点的曲率小于最大偏移时间
          {
            pt.curvature += (frame_header_time - meas.last_lio_update_time) * 1000.0f; // 更新点的曲率
            meas.pcl_proc_cur->points.push_back(pt); // 将点添加到当前处理点云
          }
          else // 否则
          {
            pt.curvature += (frame_header_time - m.lio_time) * 1000.0f; // 更新点的曲率
            meas.pcl_proc_next->points.push_back(pt); // 将点添加到下一处理点云
          }
        }
        lid_raw_data_buffer.pop_front(); // 从激光雷达缓冲区中移除第一个数据
        lid_header_time_buffer.pop_front(); // 从时间戳缓冲区中移除第一个数据
      }

      meas.measures.push_back(m); // 将测量组添加到测量数据中
      meas.lio_vio_flg = LIO; // 设置LIO-VIO标志
      // meas.last_lio_update_time = m.lio_time;
      // printf("!!! meas.lio_vio_flg: %d \n", meas.lio_vio_flg);
      // printf("[ Data Cut ] pcl_proc_cur number: %d \n", meas.pcl_proc_cur
      // ->points.size()); printf("[ Data Cut ] LIO process time: %lf \n",
      // omp_get_wtime() - t0);
      return true; // 返回真
    }

    case LIO: // LIO模式
    {
      double img_capture_time = img_time_buffer.front() + exposure_time_init; // 获取图像捕获时间
      meas.lio_vio_flg = VIO; // 设置LIO-VIO标志
      // printf("[ Data Cut ] VIO \n");
      meas.measures.clear(); // 清空测量数据
      double imu_time = imu_buffer.front()->header.stamp.toSec(); // 获取IMU时间

      struct MeasureGroup m; // 创建测量组
      m.vio_time = img_capture_time; // 设置VIO时间
      m.lio_time = meas.last_lio_update_time; // 设置LIO时间
      m.img = img_buffer.front(); // 设置图像
      mtx_buffer.lock(); // 锁定缓冲区
      // while ((!imu_buffer.empty() && (imu_time < img_capture_time)))
      // {
      //   imu_time = imu_buffer.front()->header.stamp.toSec();
      //   if (imu_time > img_capture_time) break;
      //   m.imu.push_back(imu_buffer.front());
      //   imu_buffer.pop_front();
      //   printf("[ Data Cut ] imu time: %lf \n",
      //   imu_buffer.front()->header.stamp.toSec());
      // }
      img_buffer.pop_front(); // 从图像缓冲区中移除图像数据
      img_time_buffer.pop_front(); // 从时间戳缓冲区中移除时间戳
      mtx_buffer.unlock(); // 解锁缓冲区
      sig_buffer.notify_all(); // 通知所有等待的线程
      meas.measures.push_back(m); // 将测量组添加到测量数据中
      lidar_pushed = false; // 在VIO更新后，激光雷达帧结束时间将被刷新
      // printf("[ Data Cut ] VIO process time: %lf \n", omp_get_wtime() - t0);
      return true; // 返回真
    }

    default: // 默认情况
    {
      // printf("!! WRONG EKF STATE !!");
      return false; // 返回假
    }
      // return false;
    }
    break;
  }

  default: // 默认SLAM类型
  {
    printf("!! WRONG SLAM TYPE !!");
    return false; // 返回假
  }
  }
  ROS_ERROR("out sync"); // 输出错误信息
}
// 发布RGB图像的函数
void LIVMapper::publish_img_rgb(const image_transport::Publisher &pubImage, VIOManagerPtr vio_manager)
{
  cv::Mat img_rgb = vio_manager->img_cp; // 从VIO管理器获取RGB图像
  cv_bridge::CvImage out_msg; // 创建cv_bridge图像消息
  out_msg.header.stamp = ros::Time::now(); // 设置时间戳
  // out_msg.header.frame_id = "camera_init"; // 设置坐标系（可以选择性使用）
  out_msg.encoding = sensor_msgs::image_encodings::BGR8; // 设置图像编码格式
  out_msg.image = img_rgb; // 将RGB图像赋值给消息
  pubImage.publish(out_msg.toImageMsg()); // 发布图像消息
}

// 发布全局坐标系下的激光雷达帧
void LIVMapper::publish_frame_world(const ros::Publisher &pubLaserCloudFullRes, VIOManagerPtr vio_manager)
{
  if (pcl_w_wait_pub->empty()) return; // 如果待发布的点云为空则返回
  PointCloudXYZRGB::Ptr laserCloudWorldRGB(new PointCloudXYZRGB()); // 创建新的点云对象
  if (img_en) // 如果启用图像
  {
    static int pub_num = 1; // 发布计数器
    *pcl_wait_pub += *pcl_w_wait_pub; // 将待发布的点云添加到发布缓冲区
    if(pub_num == pub_scan_num) // 如果达到发布数量
    {
      pub_num = 1; // 重置计数器
      size_t size = pcl_wait_pub->points.size(); // 获取点云大小
      laserCloudWorldRGB->reserve(size); // 预留点云空间
      cv::Mat img_rgb = vio_manager->img_rgb; // 获取RGB图像
      for (size_t i = 0; i < size; i++) // 遍历点云
      {
        PointTypeRGB pointRGB; // 创建RGB点
        pointRGB.x = pcl_wait_pub->points[i].x; // 设置X坐标
        pointRGB.y = pcl_wait_pub->points[i].y; // 设置Y坐标
        pointRGB.z = pcl_wait_pub->points[i].z; // 设置Z坐标

        V3D p_w(pcl_wait_pub->points[i].x, pcl_wait_pub->points[i].y, pcl_wait_pub->points[i].z); // 创建点的世界坐标
        V3D pf(vio_manager->new_frame_->w2f(p_w)); if (pf[2] < 0) continue; // 转换为相机坐标并检查深度
        V2D pc(vio_manager->new_frame_->w2c(p_w)); // 将世界坐标转换为图像坐标

        if (vio_manager->new_frame_->cam_->isInFrame(pc.cast<int>(), 3)) // 检查像素是否在图像帧内
        {
          V3F pixel = vio_manager->getInterpolatedPixel(img_rgb, pc); // 获取插值像素值
          pointRGB.r = pixel[2]; // 设置红色通道
          pointRGB.g = pixel[1]; // 设置绿色通道
          pointRGB.b = pixel[0]; // 设置蓝色通道
          // 以下注释部分可以用于强度归一化
          // pointRGB.r = pixel[2] * inv_expo; pointRGB.g = pixel[1] * inv_expo; pointRGB.b = pixel[0] * inv_expo;
          // if (pointRGB.r > 255) pointRGB.r = 255;
          // else if (pointRGB.r < 0) pointRGB.r = 0;
          // if (pointRGB.g > 255) pointRGB.g = 255;
          // else if (pointRGB.g < 0) pointRGB.g = 0;
          // if (pointRGB.b > 255) pointRGB.b = 255;
          // else if (pointRGB.b < 0) pointRGB.b = 0;
          if (pf.norm() > blind_rgb_points) laserCloudWorldRGB->push_back(pointRGB); // 如果点距离大于阈值则添加到点云
        }
      }
    }
    else
    {
      pub_num++; // 增加发布计数
    }
  }

  /*** 发布帧 ***/
  sensor_msgs::PointCloud2 laserCloudmsg; // 创建ROS点云消息
  if (img_en) // 如果启用图像
  {
    // cout << "RGB pointcloud size: " << laserCloudWorldRGB->size() << endl;
    pcl::toROSMsg(*laserCloudWorldRGB, laserCloudmsg); // 将RGB点云转换为ROS消息
  }
  else 
  { 
    pcl::toROSMsg(*pcl_w_wait_pub, laserCloudmsg); // 如果不启用图像则发布待发布的点云
  }
  laserCloudmsg.header.stamp = ros::Time::now(); //.fromSec(last_timestamp_lidar); // 设置时间戳
  laserCloudmsg.header.frame_id = "camera_init"; // 设置坐标系
  pubLaserCloudFullRes.publish(laserCloudmsg); // 发布激光雷达点云消息

  /**************** 保存地图 ****************/
  /* 1. 确保有足够的内存
  /* 2. 注意PCD保存将影响实时性能 **/
  if (pcd_save_en) // 如果启用PCD保存
  {
    int size = feats_undistort->points.size(); // 获取特征点数量
    PointCloudXYZI::Ptr laserCloudWorld(new PointCloudXYZI(size, 1)); // 创建新的点云对象
    static int scan_wait_num = 0; // 扫描等待计数器

    *pcl_wait_save += *laserCloudWorldRGB; // 将RGB点云添加到保存缓冲区
    scan_wait_num++; // 增加扫描等待计数

    if (pcl_wait_save->size() > 0 && pcd_save_interval > 0 && scan_wait_num >= pcd_save_interval) // 如果满足保存条件
    {
      pcd_index++; // 增加PCD索引
      string all_points_dir(string(string(ROOT_DIR) + "Log/PCD/") + to_string(pcd_index) + string(".pcd")); // 构建保存路径
      pcl::PCDWriter pcd_writer; // 创建PCD写入器
      if (pcd_save_en) // 如果启用PCD保存
      {
        cout << "current scan saved to /PCD/" << all_points_dir << endl; // 输出保存信息
        pcd_writer.writeBinary(all_points_dir, *pcl_wait_save); // 保存点云为二进制格式
        PointCloudXYZRGB().swap(*pcl_wait_save); // 清空保存缓冲区
        Eigen::Quaterniond q(_state.rot_end); // 获取当前旋转状态
        fout_pcd_pos << _state.pos_end[0] << " " << _state.pos_end[1] << " " << _state.pos_end[2] << " " << q.w() << " " << q.x() << " " << q.y()
                     << " " << q.z() << " " << endl; // 保存位置和姿态信息
        scan_wait_num = 0; // 重置扫描等待计数
      }
    }
  }
  if(laserCloudWorldRGB->size() > 0) // 如果RGB点云有数据
  {
    PointCloudXYZI().swap(*pcl_wait_pub); // 清空待发布点云
    PointCloudXYZRGB().swap(*laserCloudWorldRGB); // 清空RGB点云
  }
  PointCloudXYZI().swap(*pcl_w_wait_pub); // 清空待发布的点云
}

// 发布视觉子地图
void LIVMapper::publish_visual_sub_map(const ros::Publisher &pubSubVisualMap)
{
  PointCloudXYZI::Ptr laserCloudFullRes(visual_sub_map); // 获取视觉子地图
  int size = laserCloudFullRes->points.size(); if (size == 0) return; // 如果子地图为空则返回
  PointCloudXYZI::Ptr sub_pcl_visual_map_pub(new PointCloudXYZI()); // 创建新的点云对象
  *sub_pcl_visual_map_pub = *laserCloudFullRes; // 复制点云数据
  if (1) // 发送条件
  {
    sensor_msgs::PointCloud2 laserCloudmsg; // 创建ROS点云消息
    pcl::toROSMsg(*sub_pcl_visual_map_pub, laserCloudmsg); // 转换为ROS消息
    laserCloudmsg.header.stamp = ros::Time::now(); // 设置时间戳
    laserCloudmsg.header.frame_id = "camera_init"; // 设置坐标系
    pubSubVisualMap.publish(laserCloudmsg); // 发布视觉子地图
  }
}

// 发布效果点云
void LIVMapper::publish_effect_world(const ros::Publisher &pubLaserCloudEffect, const std::vector<PointToPlane> &ptpl_list)
{
  int effect_feat_num = ptpl_list.size(); // 获取效果特征数量
  PointCloudXYZI::Ptr laserCloudWorld(new PointCloudXYZI(effect_feat_num, 1)); // 创建新的点云对象
  for (int i = 0; i < effect_feat_num; i++) // 遍历特征列表
  {
    laserCloudWorld->points[i].x = ptpl_list[i].point_w_[0]; // 设置X坐标
    laserCloudWorld->points[i].y = ptpl_list[i].point_w_[1]; // 设置Y坐标
    laserCloudWorld->points[i].z = ptpl_list[i].point_w_[2]; // 设置Z坐标
  }
  sensor_msgs::PointCloud2 laserCloudFullRes3; // 创建ROS点云消息
  pcl::toROSMsg(*laserCloudWorld, laserCloudFullRes3); // 转换为ROS消息
  laserCloudFullRes3.header.stamp = ros::Time::now(); // 设置时间戳
  laserCloudFullRes3.header.frame_id = "camera_init"; // 设置坐标系
  pubLaserCloudEffect.publish(laserCloudFullRes3); // 发布效果点云
}

// 设置位姿和时间戳的模板函数
template <typename T> void LIVMapper::set_posestamp(T &out)
{
  out.position.x = _state.pos_end(0); // 设置位置X
  out.position.y = _state.pos_end(1); // 设置位置Y
  out.position.z = _state.pos_end(2); // 设置位置Z
  out.orientation.x = geoQuat.x; // 设置姿态四元数X
  out.orientation.y = geoQuat.y; // 设置姿态四元数Y
  out.orientation.z = geoQuat.z; // 设置姿态四元数Z
  out.orientation.w = geoQuat.w; // 设置姿态四元数W
}

// 发布里程计信息
void LIVMapper::publish_odometry(const ros::Publisher &pubOdomAftMapped)
{
  odomAftMapped.header.frame_id = "camera_init"; // 设置父坐标系
  odomAftMapped.child_frame_id = "aft_mapped"; // 设置子坐标系
  odomAftMapped.header.stamp = ros::Time::now(); //.ros::Time()fromSec(last_timestamp_lidar); // 设置时间戳
  set_posestamp(odomAftMapped.pose.pose); // 设置位姿和时间戳

  static tf::TransformBroadcaster br; // 创建TF广播器
  tf::Transform transform; // 创建TF变换
  tf::Quaternion q; // 创建四元数
  transform.setOrigin(tf::Vector3(_state.pos_end(0), _state.pos_end(1), _state.pos_end(2))); // 设置位置
  q.setW(geoQuat.w); // 设置四元数W
  q.setX(geoQuat.x); // 设置四元数X
  q.setY(geoQuat.y); // 设置四元数Y
  q.setZ(geoQuat.z); // 设置四元数Z
  transform.setRotation(q); // 设置旋转
  br.sendTransform( tf::StampedTransform(transform, odomAftMapped.header.stamp, "camera_init", "aft_mapped") ); // 发送TF变换
  pubOdomAftMapped.publish(odomAftMapped); // 发布里程计信息
}

// 发布Mavros位姿信息
void LIVMapper::publish_mavros(const ros::Publisher &mavros_pose_publisher)
{
  msg_body_pose.header.stamp = ros::Time::now(); // 设置时间戳
  msg_body_pose.header.frame_id = "camera_init"; // 设置坐标系
  set_posestamp(msg_body_pose.pose); // 设置位姿和时间戳
  mavros_pose_publisher.publish(msg_body_pose); // 发布Mavros位姿信息
}

// 发布路径信息
void LIVMapper::publish_path(const ros::Publisher pubPath)
{
  set_posestamp(msg_body_pose.pose); // 设置位姿
  msg_body_pose.header.stamp = ros::Time::now(); // 设置时间戳
  msg_body_pose.header.frame_id = "camera_init"; // 设置坐标系
  path.poses.push_back(msg_body_pose); // 将当前位姿添加到路径中
  pubPath.publish(path); // 发布路径信息
}
