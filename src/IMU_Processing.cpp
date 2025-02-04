/* 
This file is part of FAST-LIVO2: Fast, Direct LiDAR-Inertial-Visual Odometry.

Developer: Chunran Zheng <zhengcr@connect.hku.hk>

For commercial use, please contact me at <zhengcr@connect.hku.hk> or
Prof. Fu Zhang at <fuzhang@hku.hk>.

This file is subject to the terms and conditions outlined in the 'LICENSE' file,
which is included as part of this source code package.
*/
#include "IMU_Processing.h"

// 比较两个点的曲率，返回值用于排序
const bool time_list(PointType &x, PointType &y) { return (x.curvature < y.curvature); }

// IMU处理类的构造函数
ImuProcess::ImuProcess() : Eye3d(M3D::Identity()), // 初始化3D单位矩阵
                           Zero3d(0, 0, 0), // 初始化零向量
                           b_first_frame(true), // 标记是否为第一帧
                           imu_need_init(true) // 标记IMU是否需要初始化
{
  init_iter_num = 1; // 初始化迭代次数
  // 初始化协方差
  cov_acc = V3D(0.1, 0.1, 0.1);
  cov_gyr = V3D(0.1, 0.1, 0.1);
  cov_bias_gyr = V3D(0.1, 0.1, 0.1);
  cov_bias_acc = V3D(0.1, 0.1, 0.1);
  cov_inv_expo = 0.2; // 初始化曝光逆协方差
  mean_acc = V3D(0, 0, -1.0); // 初始化加速度均值
  mean_gyr = V3D(0, 0, 0); // 初始化角速度均值
  angvel_last = Zero3d; // 上一角速度
  acc_s_last = Zero3d; // 上一加速度
  Lid_offset_to_IMU = Zero3d; // 激光传感器到IMU的偏移
  Lid_rot_to_IMU = Eye3d; // 激光传感器到IMU的旋转
  last_imu.reset(new sensor_msgs::Imu()); // 初始化IMU消息
  cur_pcl_un_.reset(new PointCloudXYZI()); // 初始化点云
}

// IMU处理类的析构函数
ImuProcess::~ImuProcess() {}

// 重置IMU处理器
void ImuProcess::Reset()
{
  ROS_WARN("Reset ImuProcess"); // 输出警告信息
  mean_acc = V3D(0, 0, -1.0); // 重置加速度均值
  mean_gyr = V3D(0, 0, 0); // 重置角速度均值
  angvel_last = Zero3d; // 重置上一个角速度
  imu_need_init = true; // 设置IMU需要初始化
  init_iter_num = 1; // 重置初始化迭代次数
  IMUpose.clear(); // 清空IMU位姿
  last_imu.reset(new sensor_msgs::Imu()); // 重置IMU消息
  cur_pcl_un_.reset(new PointCloudXYZI()); // 重置点云
}

// 禁用IMU
void ImuProcess::disable_imu()
{
  cout << "IMU Disabled !!!!!" << endl; // 输出禁用信息
  imu_en = false; // 设置IMU启用状态为假
  imu_need_init = false; // 设置IMU不需要初始化
}

// 禁用重力估计
void ImuProcess::disable_gravity_est()
{
  cout << "Online Gravity Estimation Disabled !!!!!" << endl; // 输出禁用信息
  gravity_est_en = false; // 设置重力估计启用状态为假
}

// 禁用偏置估计
void ImuProcess::disable_bias_est()
{
  cout << "Bias Estimation Disabled !!!!!" << endl; // 输出禁用信息
  ba_bg_est_en = false; // 设置偏置估计启用状态为假
}

// 禁用曝光估计
void ImuProcess::disable_exposure_est()
{
  cout << "Online Time Offset Estimation Disabled !!!!!" << endl; // 输出禁用信息
  exposure_estimate_en = false; // 设置曝光估计启用状态为假
}

// 设置外参
void ImuProcess::set_extrinsic(const MD(4, 4) & T)
{
  Lid_offset_to_IMU = T.block<3, 1>(0, 3); // 设置激光传感器到IMU的平移
  Lid_rot_to_IMU = T.block<3, 3>(0, 0); // 设置激光传感器到IMU的旋转
}

// 设置外参（平移）
void ImuProcess::set_extrinsic(const V3D &transl)
{
  Lid_offset_to_IMU = transl; // 设置激光传感器到IMU的平移
  Lid_rot_to_IMU.setIdentity(); // 旋转设置为单位矩阵
}

// 设置外参（平移和旋转）
void ImuProcess::set_extrinsic(const V3D &transl, const M3D &rot)
{
  Lid_offset_to_IMU = transl; // 设置激光传感器到IMU的平移
  Lid_rot_to_IMU = rot; // 设置激光传感器到IMU的旋转
}

// 设置陀螺仪协方差缩放
void ImuProcess::set_gyr_cov_scale(const V3D &scaler) { cov_gyr = scaler; }

// 设置加速度协方差缩放
void ImuProcess::set_acc_cov_scale(const V3D &scaler) { cov_acc = scaler; }

// 设置陀螺仪偏置协方差
void ImuProcess::set_gyr_bias_cov(const V3D &b_g) { cov_bias_gyr = b_g; }

// 设置曝光逆协方差
void ImuProcess::set_inv_expo_cov(const double &inv_expo) { cov_inv_expo = inv_expo; }

// 设置加速度偏置协方差
void ImuProcess::set_acc_bias_cov(const V3D &b_a) { cov_bias_acc = b_a; }

// 设置IMU初始化帧数
void ImuProcess::set_imu_init_frame_num(const int &num) { MAX_INI_COUNT = num; }

// IMU初始化函数
void ImuProcess::IMU_init(const MeasureGroup &meas, StatesGroup &state_inout, int &N)
{
  /** 1. 初始化重力、陀螺仪偏置、加速度和陀螺仪协方差
   ** 2. 将加速度测量值归一化为单位重力 **/
  ROS_INFO("IMU Initializing: %.1f %%", double(N) / MAX_INI_COUNT * 100); // 输出初始化进度
  V3D cur_acc, cur_gyr; // 当前加速度和角速度

  if (b_first_frame) // 如果是第一帧
  {
    Reset(); // 重置IMU处理器
    N = 1; // 初始化计数
    b_first_frame = false; // 设置为非第一帧
    const auto &imu_acc = meas.imu.front()->linear_acceleration; // 获取第一个IMU的线性加速度
    const auto &gyr_acc = meas.imu.front()->angular_velocity; // 获取第一个IMU的角速度
    mean_acc << imu_acc.x, imu_acc.y, imu_acc.z; // 初始化加速度均值
    mean_gyr << gyr_acc.x, gyr_acc.y, gyr_acc.z; // 初始化角速度均值
  }

  // 遍历所有IMU测量数据
  for (const auto &imu : meas.imu)
  {
    const auto &imu_acc = imu->linear_acceleration; // 获取当前IMU的线性加速度
    const auto &gyr_acc = imu->angular_velocity; // 获取当前IMU的角速度
    cur_acc << imu_acc.x, imu_acc.y, imu_acc.z; // 更新当前加速度
    cur_gyr << gyr_acc.x, gyr_acc.y, gyr_acc.z; // 更新当前角速度

    // 更新加速度均值
    mean_acc += (cur_acc - mean_acc) / N;
    // 更新角速度均值
    mean_gyr += (cur_gyr - mean_gyr) / N;

    // N++; // 计数增加
    N++; // 计数增加
  }
  
  // 计算IMU加速度均值的范数
  IMU_mean_acc_norm = mean_acc.norm();
  // 更新状态中的重力
  state_inout.gravity = -mean_acc / mean_acc.norm() * G_m_s2; // 归一化加速度并标定重力
  state_inout.rot_end = Eye3d; // 初始化旋转矩阵
  state_inout.bias_g = Zero3d; // 初始化偏置

  last_imu = meas.imu.back(); // 更新最后的IMU数据
}
// IMU处理类的前向传播函数（不使用IMU数据）
void ImuProcess::Forward_without_imu(LidarMeasureGroup &meas, StatesGroup &state_inout, PointCloudXYZI &pcl_out)
{
  const double &pcl_beg_time = meas.lidar_frame_beg_time; // 获取点云开始时间

  /*** 按偏移时间对点云进行排序 ***/
  pcl_out = *(meas.lidar); // 将当前点云数据赋值给输出点云
  // sort(pcl_out->points.begin(), pcl_out->points.end(), time_list);
  const double &pcl_end_time = pcl_beg_time + pcl_out.points.back().curvature / double(1000); // 计算点云结束时间
  meas.last_lio_update_time = pcl_end_time; // 更新最后的LIO更新时间
  MD(DIM_STATE, DIM_STATE) F_x, cov_w; // 状态转移矩阵和协方差矩阵
  double dt = 0; // 时间间隔初始化

  if (b_first_frame) // 如果是第一帧
  {
    dt = 0.1; // 设置初始时间间隔
    b_first_frame = false; // 设置为非第一帧
  }
  else { dt = pcl_beg_time - time_last_scan; } // 计算时间间隔

  time_last_scan = pcl_beg_time; // 更新最后扫描时间

  /* 协方差传播 */
  M3D Exp_f = Exp(state_inout.bias_g, dt); // 计算状态转移矩阵的指数形式

  F_x.setIdentity(); // 初始化状态转移矩阵为单位矩阵
  cov_w.setZero(); // 初始化协方差矩阵为零矩阵

  F_x.block<3, 3>(0, 0) = Exp(state_inout.bias_g, -dt); // 更新状态转移矩阵的旋转部分
  F_x.block<3, 3>(0, 9) = Eye3d * dt; // 更新状态转移矩阵的速度部分
  F_x.block<3, 3>(3, 6) = Eye3d * dt; // 更新状态转移矩阵的加速度部分

  cov_w.block<3, 3>(9, 9).diagonal() = cov_gyr * dt * dt; // 为角速度的协方差设置值
  cov_w.block<3, 3>(6, 6).diagonal() = cov_acc * dt * dt; // 为速度的协方差设置值

  state_inout.cov = F_x * state_inout.cov * F_x.transpose() + cov_w; // 协方差传播
  state_inout.rot_end = state_inout.rot_end * Exp_f; // 更新旋转状态
  state_inout.pos_end = state_inout.pos_end + state_inout.vel_end * dt; // 更新位置
}

// IMU点云去畸变处理
void ImuProcess::UndistortPcl(LidarMeasureGroup &lidar_meas, StatesGroup &state_inout, PointCloudXYZI &pcl_out)
{
  double t0 = omp_get_wtime(); // 记录开始时间
  pcl_out.clear(); // 清空输出点云
  /*** 将上一帧的IMU数据添加到当前帧的开头 ***/
  MeasureGroup &meas = lidar_meas.measures.back(); // 获取最新的测量数据
  auto v_imu = meas.imu; // 获取IMU数据
  v_imu.push_front(last_imu); // 将最后的IMU数据插入到IMU数据的开头
  const double &imu_beg_time = v_imu.front()->header.stamp.toSec(); // 获取IMU开始时间
  const double &imu_end_time = v_imu.back()->header.stamp.toSec(); // 获取IMU结束时间
  const double prop_beg_time = last_prop_end_time; // 获取上次传播结束时间

  const double prop_end_time = lidar_meas.lio_vio_flg == LIO ? meas.lio_time : meas.vio_time; // 计算传播结束时间

  /*** 根据传播开始时间和所需的传播结束时间裁剪激光点 ***/
  if (lidar_meas.lio_vio_flg == LIO)
  {
    pcl_wait_proc.resize(lidar_meas.pcl_proc_cur->points.size()); // 调整待处理点云的大小
    pcl_wait_proc = *(lidar_meas.pcl_proc_cur); // 复制当前点云到待处理点云
    lidar_meas.lidar_scan_index_now = 0; // 重置激光扫描索引
    IMUpose.push_back(set_pose6d(0.0, acc_s_last, angvel_last, state_inout.vel_end, state_inout.pos_end, state_inout.rot_end)); // 保存IMU位姿
  }

  /*** 初始化IMU位姿 ***/
  V3D acc_imu(acc_s_last), angvel_avr(angvel_last), acc_avr, vel_imu(state_inout.vel_end), pos_imu(state_inout.pos_end); // 初始化加速度、角速度、速度和位置
  M3D R_imu(state_inout.rot_end); // 初始化旋转矩阵
  MD(DIM_STATE, DIM_STATE) F_x, cov_w; // 状态转移矩阵和协方差矩阵
  double dt, dt_all = 0.0; // 时间间隔和总时间间隔
  double offs_t;
  double tau; // 时间常数

  if (!imu_time_init) // 如果IMU时间未初始化
  {
    tau = 1.0; // 设置时间常数
    imu_time_init = true; // 设置IMU时间初始化为真
  }
  else
  {
    tau = state_inout.inv_expo_time; // 获取逆曝光时间
  }

  switch (lidar_meas.lio_vio_flg) // 根据LIO/VIO标志进行处理
  {
  case LIO:
  case VIO:
    dt = 0; // 时间间隔初始化
    for (int i = 0; i < v_imu.size() - 1; i++) // 遍历IMU数据
    {
      auto head = v_imu[i]; // 当前IMU数据
      auto tail = v_imu[i + 1]; // 下一个IMU数据

      if (tail->header.stamp.toSec() < prop_beg_time) continue; // 如果下一个IMU时间小于传播开始时间，跳过

      angvel_avr << 0.5 * (head->angular_velocity.x + tail->angular_velocity.x), 0.5 * (head->angular_velocity.y + tail->angular_velocity.y),
          0.5 * (head->angular_velocity.z + tail->angular_velocity.z); // 计算平均角速度

      acc_avr << 0.5 * (head->linear_acceleration.x + tail->linear_acceleration.x), 0.5 * (head->linear_acceleration.y + tail->linear_acceleration.y),
          0.5 * (head->linear_acceleration.z + tail->linear_acceleration.z); // 计算平均加速度

      angvel_avr -= state_inout.bias_g; // 减去陀螺仪偏置
      acc_avr = acc_avr * G_m_s2 / mean_acc.norm() - state_inout.bias_a; // 归一化加速度并减去加速度偏置

      // 计算时间间隔
      if (head->header.stamp.toSec() < prop_beg_time)
      {
        dt = tail->header.stamp.toSec() - last_prop_end_time; // 计算时间间隔
        offs_t = tail->header.stamp.toSec() - prop_beg_time;
      }
      else if (i != v_imu.size() - 2)
      {
        dt = tail->header.stamp.toSec() - head->header.stamp.toSec(); // 计算时间间隔
        offs_t = tail->header.stamp.toSec() - prop_beg_time;
      }
      else
      {
        dt = prop_end_time - head->header.stamp.toSec(); // 计算时间间隔
        offs_t = prop_end_time - prop_beg_time;
      }

      dt_all += dt; // 累加总时间间隔

      /* 协方差传播 */
      M3D acc_avr_skew; // 加速度的反对称矩阵
      M3D Exp_f = Exp(angvel_avr, dt); // 计算状态转移矩阵的指数形式
      acc_avr_skew << SKEW_SYM_MATRX(acc_avr); // 构造加速度的反对称矩阵

      F_x.setIdentity(); // 初始化状态转移矩阵为单位矩阵
      cov_w.setZero(); // 初始化协方差矩阵为零矩阵

      F_x.block<3, 3>(0, 0) = Exp(angvel_avr, -dt); // 更新状态转移矩阵的旋转部分
      if (ba_bg_est_en) F_x.block<3, 3>(0, 10) = -Eye3d * dt; // 更新状态转移矩阵的偏置部分
      F_x.block<3, 3>(3, 7) = Eye3d * dt; // 更新状态转移矩阵的加速度部分
      F_x.block<3, 3>(7, 0) = -R_imu * acc_avr_skew * dt; // 更新状态转移矩阵的加速度部分
      if (gravity_est_en) F_x.block<3, 3>(7, 16) = Eye3d * dt; // 更新状态转移矩阵的重力部分

      if (exposure_estimate_en) cov_w(6, 6) = cov_inv_expo * dt * dt; // 更新曝光协方差
      cov_w.block<3, 3>(0, 0).diagonal() = cov_gyr * dt * dt; // 更新角速度的协方差
      cov_w.block<3, 3>(7, 7) = R_imu * cov_acc.asDiagonal() * R_imu.transpose() * dt * dt; // 更新加速度的协方差
      cov_w.block<3, 3>(10, 10).diagonal() = cov_bias_gyr * dt * dt; // 更新陀螺仪偏置协方差
      cov_w.block<3, 3>(13, 13).diagonal() = cov_bias_acc * dt * dt; // 更新加速度偏置协方差

      state_inout.cov = F_x * state_inout.cov * F_x.transpose() + cov_w; // 协方差传播

      /* IMU姿态传播 */
      R_imu = R_imu * Exp_f; // 更新旋转矩阵

      /* IMU的特定加速度（全局坐标系） */
      acc_imu = R_imu * acc_avr + state_inout.gravity; // 计算IMU的全局加速度

      /* IMU的传播 */
      pos_imu = pos_imu + vel_imu * dt + 0.5 * acc_imu * dt * dt; // 更新位置
      vel_imu = vel_imu + acc_imu * dt; // 更新速度

      /* 保存每个IMU测量的位姿 */
      angvel_last = angvel_avr; // 更新上一个角速度
      acc_s_last = acc_imu; // 更新上一个加速度

      IMUpose.push_back(set_pose6d(offs_t, acc_imu, angvel_avr, vel_imu, pos_imu, R_imu)); // 保存IMU位姿
    }

    lidar_meas.last_lio_update_time = prop_end_time; // 更新最后的LIO更新时间
    break;
  }

  state_inout.vel_end = vel_imu; // 更新状态中的速度
  state_inout.rot_end = R_imu; // 更新状态中的旋转
  state_inout.pos_end = pos_imu; // 更新状态中的位置
  state_inout.inv_expo_time = tau; // 更新状态中的逆曝光时间

  last_imu = v_imu.back(); // 更新最后的IMU数据
  last_prop_end_time = prop_end_time; // 更新最后的传播结束时间

  double t1 = omp_get_wtime(); // 记录结束时间

  /*** 反向传播每个激光点（去畸变处理），仅适用于LIO更新 ***/
  if (lidar_meas.lio_vio_flg == LIO)
  {
    auto it_pcl = pcl_wait_proc.points.end() - 1; // 获取待处理点云的最后一个点
    M3D extR_Ri(Lid_rot_to_IMU.transpose() * state_inout.rot_end.transpose()); // 计算外部旋转
    V3D exrR_extT(Lid_rot_to_IMU.transpose() * Lid_offset_to_IMU); // 计算外部平移
    for (auto it_kp = IMUpose.end() - 1; it_kp != IMUpose.begin(); it_kp--) // 反向遍历IMU位姿
    {
      auto head = it_kp - 1; // 当前IMU位姿
      auto tail = it_kp; // 下一IMU位姿
      R_imu << MAT_FROM_ARRAY(head->rot); // 获取当前IMU的旋转矩阵
      acc_imu << VEC_FROM_ARRAY(head->acc); // 获取当前IMU的加速度
      vel_imu << VEC_FROM_ARRAY(head->vel); // 获取当前IMU的速度
      pos_imu << VEC_FROM_ARRAY(head->pos); // 获取当前IMU的位置
      angvel_avr << VEC_FROM_ARRAY(head->gyr); // 获取当前IMU的角速度

      for (; it_pcl->curvature / double(1000) > head->offset_time; it_pcl--) // 遍历待处理点云
      {
        dt = it_pcl->curvature / double(1000) - head->offset_time; // 计算时间间隔

        /* 转换到“结束”坐标系 */
        M3D R_i(R_imu * Exp(angvel_avr, dt)); // 计算当前点的旋转
        V3D T_ei(pos_imu + vel_imu * dt + 0.5 * acc_imu * dt * dt - state_inout.pos_end); // 计算当前点的平移

        V3D P_i(it_pcl->x, it_pcl->y, it_pcl->z); // 获取当前点的坐标
        V3D P_compensate = (extR_Ri * (R_i * (Lid_rot_to_IMU * P_i + Lid_offset_to_IMU) + T_ei) - exrR_extT); // 计算补偿后的点

        // 保存去畸变后的点及其旋转
        it_pcl->x = P_compensate(0); // 更新x坐标
        it_pcl->y = P_compensate(1); // 更新y坐标
        it_pcl->z = P_compensate(2); // 更新z坐标

        if (it_pcl == pcl_wait_proc.points.begin()) break; // 如果到达待处理点云的开头，退出循环
      }
    }
    pcl_out = pcl_wait_proc; // 输出去畸变后的点云
    pcl_wait_proc.clear(); // 清空待处理点云
    IMUpose.clear(); // 清空IMU位姿
  }
}

// 处理激光测量数据
void ImuProcess::Process2(LidarMeasureGroup &lidar_meas, StatesGroup &stat, PointCloudXYZI::Ptr cur_pcl_un_)
{
  double t1, t2, t3;
  t1 = omp_get_wtime(); // 记录开始时间
  ROS_ASSERT(lidar_meas.lidar != nullptr); // 确保激光数据不为空
  if (!imu_en) // 如果IMU未启用
  {
    Forward_without_imu(lidar_meas, stat, *cur_pcl_un_); // 调用不使用IMU的前向传播函数
    return;
  }

  MeasureGroup meas = lidar_meas.measures.back(); // 获取最新的测量数据

  if (imu_need_init) // 如果需要初始化IMU
  {
    double pcl_end_time = lidar_meas.lio_vio_flg == LIO ? meas.lio_time : meas.vio_time; // 获取点云结束时间

    if (meas.imu.empty()) { return; }; // 如果IMU数据为空，返回

    // 第一个激光帧
    IMU_init(meas, stat, init_iter_num); // 初始化IMU

    imu_need_init = true; // 设置IMU需要初始化为真

    last_imu = meas.imu.back(); // 更新最后的IMU数据

    if (init_iter_num > MAX_INI_COUNT) // 如果初始化迭代次数超过最大值
    {
      imu_need_init = false; // 设置IMU不需要初始化
      ROS_INFO("IMU Initials: Gravity: %.4f %.4f %.4f %.4f; acc covarience: "
               "%.8f %.8f %.8f; gry covarience: %.8f %.8f %.8f \n",
               stat.gravity[0], stat.gravity[1], stat.gravity[2], mean_acc.norm(), cov_acc[0], cov_acc[1], cov_acc[2], cov_gyr[0], cov_gyr[1],
               cov_gyr[2]); // 输出IMU初始化信息
      ROS_INFO("IMU Initials: ba covarience: %.8f %.8f %.8f; bg covarience: "
               "%.8f %.8f %.8f",
               cov_bias_acc[0], cov_bias_acc[1], cov_bias_acc[2], cov_bias_gyr[0], cov_bias_gyr[1], cov_bias_gyr[2]); // 输出IMU偏置协方差信息
      fout_imu.open(DEBUG_FILE_DIR("imu.txt"), ios::out); // 打开IMU调试文件
    }

    return; // 返回
  }

  UndistortPcl(lidar_meas, stat, *cur_pcl_un_); // 调用去畸变处理函数
}
