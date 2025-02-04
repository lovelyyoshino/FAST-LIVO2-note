/* 
This file is part of FAST-LIVO2: Fast, Direct LiDAR-Inertial-Visual Odometry.

Developer: Chunran Zheng <zhengcr@connect.hku.hk>

For commercial use, please contact me at <zhengcr@connect.hku.hk> or
Prof. Fu Zhang at <fuzhang@hku.hk>.

This file is subject to the terms and conditions outlined in the 'LICENSE' file,
which is included as part of this source code package.
*/
#include "vio.h"

// VIOManager类的构造函数
VIOManager::VIOManager()
{
  // downSizeFilter.setLeafSize(0.2, 0.2, 0.2); // 初始化滤波器的叶子大小（注释掉了）
}

// VIOManager类的析构函数
VIOManager::~VIOManager()
{
  delete visual_submap; // 删除视觉子图
  for (auto& pair : warp_map) delete pair.second; // 删除变换地图中的每个元素
  warp_map.clear(); // 清空变换地图
  for (auto& pair : feat_map) delete pair.second; // 删除特征地图中的每个元素
  feat_map.clear(); // 清空特征地图
}

// 设置IMU到激光雷达的外部参数
void VIOManager::setImuToLidarExtrinsic(const V3D &transl, const M3D &rot)
{
  Pli = -rot.transpose() * transl; // 计算IMU到激光雷达的平移
  Rli = rot.transpose(); // 计算IMU到激光雷达的旋转
}

// 设置激光雷达到相机的外部参数
void VIOManager::setLidarToCameraExtrinsic(vector<double> &R, vector<double> &P)
{
  Rcl << MAT_FROM_ARRAY(R); // 从数组中设置旋转矩阵
  Pcl << VEC_FROM_ARRAY(P); // 从数组中设置平移向量
}

// 初始化VIO系统
void VIOManager::initializeVIO()
{
  visual_submap = new SubSparseMap; // 创建一个新的子稀疏地图

  // 获取相机的内参数
  fx = cam->fx();
  fy = cam->fy();
  cx = cam->cx();
  cy = cam->cy();
  image_resize_factor = cam->scale(); // 获取图像缩放因子

  printf("intrinsic: %.6lf, %.6lf, %.6lf, %.6lf\n", fx, fy, cx, cy); // 打印内参数

  // 获取图像的宽度和高度
  width = cam->width();
  height = cam->height();

  printf("width: %d, height: %d, scale: %f\n", width, height, image_resize_factor); // 打印图像尺寸和缩放因子
  Rci = Rcl * Rli; // 计算相机到IMU的旋转矩阵
  Pci = Rcl * Pli + Pcl; // 计算相机到IMU的平移向量

  V3D Pic; // 创建相机到IMU的平移向量
  M3D tmp; // 创建临时矩阵
  Jdphi_dR = Rci; // 更新雅可比矩阵
  Pic = -Rci.transpose() * Pci; // 计算相机的平移
  tmp << SKEW_SYM_MATRX(Pic); // 计算平移的反对称矩阵
  Jdp_dR = -Rci * tmp; // 更新雅可比矩阵

  // 网格初始化
  if (grid_size > 10)
  {
    grid_n_width = ceil(static_cast<double>(width / grid_size)); // 计算网格宽度
    grid_n_height = ceil(static_cast<double>(height / grid_size)); // 计算网格高度
  }
  else
  {
    grid_size = static_cast<int>(height / grid_n_height); // 更新网格大小
    grid_n_height = ceil(static_cast<double>(height / grid_size)); // 计算网格高度
    grid_n_width = ceil(static_cast<double>(width / grid_size)); // 计算网格宽度
  }
  length = grid_n_width * grid_n_height; // 计算网格总长度

  // 如果启用光线投射
  if(raycast_en)
  {
    border_flag.resize(length, 0); // 初始化边界标志

    std::vector<std::vector<V3D>>().swap(rays_with_sample_points); // 清空样本点向量
    rays_with_sample_points.reserve(length); // 预留空间
    printf("grid_size: %d, grid_n_height: %d, grid_n_width: %d, length: %d\n", grid_size, grid_n_height, grid_n_width, length); // 打印网格信息

    float d_min = 0.1; // 最小深度
    float d_max = 3.0; // 最大深度
    float step = 0.2; // 步长
    for (int grid_row = 1; grid_row <= grid_n_height; grid_row++)
    {
      for (int grid_col = 1; grid_col <= grid_n_width; grid_col++)
      {
        std::vector<V3D> SamplePointsEachGrid; // 每个网格的样本点
        int index = (grid_row - 1) * grid_n_width + grid_col - 1; // 计算索引

        // 判断是否为边界
        if (grid_row == 1 || grid_col == 1 || grid_row == grid_n_height || grid_col == grid_n_width) border_flag[index] = 1;

        // 计算网格中心
        int u = grid_size / 2 + (grid_col - 1) * grid_size;
        int v = grid_size / 2 + (grid_row - 1) * grid_size;
        for (float d_temp = d_min; d_temp <= d_max; d_temp += step) // 遍历深度
        {
          V3D xyz;
          xyz = cam->cam2world(u, v); // 从相机坐标转换到世界坐标
          xyz *= d_temp / xyz[2]; // 计算实际坐标
          SamplePointsEachGrid.push_back(xyz); // 添加样本点
        }
        rays_with_sample_points.push_back(SamplePointsEachGrid); // 将样本点添加到光线样本点集合
      }
    }
  }

  // 如果启用COLMAP输出
  if(colmap_output_en)
  {
    pinhole_cam = dynamic_cast<vk::PinholeCamera*>(cam); // 动态转换为针孔相机
    fout_colmap.open(DEBUG_FILE_DIR("Colmap/sparse/0/images.txt"), ios::out); // 打开图像文件
    fout_colmap << "# Image list with two lines of data per image:\n";
    fout_colmap << "#   IMAGE_ID, QW, QX, QY, QZ, TX, TY, TZ, CAMERA_ID, NAME\n";
    fout_colmap << "#   POINTS2D[] as (X, Y, POINT3D_ID)\n";
    fout_camera.open(DEBUG_FILE_DIR("Colmap/sparse/0/cameras.txt"), ios::out); // 打开相机文件
    fout_camera << "# Camera list with one line of data per camera:\n";
    fout_camera << "#   CAMERA_ID, MODEL, WIDTH, HEIGHT, PARAMS[]\n";
    fout_camera << "1 PINHOLE " << width << " " << height << " "
        << std::fixed << std::setprecision(6)  // 控制浮点数精度为6位
        << fx << " " << fy << " "
        << cx << " " << cy << std::endl;
    fout_camera.close(); // 关闭相机文件
  }
  
  // 初始化地图相关数据
  grid_num.resize(length);
  map_index.resize(length);
  map_dist.resize(length);
  update_flag.resize(length);
  scan_value.resize(length);

  patch_size_total = patch_size * patch_size; // 计算补丁总大小
  patch_size_half = static_cast<int>(patch_size / 2); // 计算补丁的一半大小
  patch_buffer.resize(patch_size_total); // 初始化补丁缓冲区
  warp_len = patch_size_total * patch_pyrimid_level; // 计算变换长度
  border = (patch_size_half + 2) * 8; // 计算边界大小

  retrieve_voxel_points.reserve(length); // 预留体素点空间
  append_voxel_points.reserve(length); // 预留附加体素点空间

  sub_feat_map.clear(); // 清空子特征地图
}

// 重置网格
void VIOManager::resetGrid()
{
  fill(grid_num.begin(), grid_num.end(), TYPE_UNKNOWN); // 重置网格数量
  fill(map_index.begin(), map_index.end(), 0); // 重置地图索引
  fill(map_dist.begin(), map_dist.end(), 10000.0f); // 重置地图距离
  fill(update_flag.begin(), update_flag.end(), 0); // 重置更新标志
  fill(scan_value.begin(), scan_value.end(), 0.0f); // 重置扫描值

  retrieve_voxel_points.clear(); // 清空检索体素点
  retrieve_voxel_points.resize(length); // 重新调整大小

  append_voxel_points.clear(); // 清空附加体素点
  append_voxel_points.resize(length); // 重新调整大小

  total_points = 0; // 重置总点数
}

// 计算投影雅可比矩阵
void VIOManager::computeProjectionJacobian(V3D p, MD(2, 3) & J)
{
  const double x = p[0]; // 获取X坐标
  const double y = p[1]; // 获取Y坐标
  const double z_inv = 1. / p[2]; // 计算Z的倒数
  const double z_inv_2 = z_inv * z_inv; // 计算Z的平方倒数
  J(0, 0) = fx * z_inv; // 计算雅可比矩阵元素
  J(0, 1) = 0.0; // 计算雅可比矩阵元素
  J(0, 2) = -fx * x * z_inv_2; // 计算雅可比矩阵元素
  J(1, 0) = 0.0; // 计算雅可比矩阵元素
  J(1, 1) = fy * z_inv; // 计算雅可比矩阵元素
  J(1, 2) = -fy * y * z_inv_2; // 计算雅可比矩阵元素
}

// 获取图像补丁
void VIOManager::getImagePatch(cv::Mat img, V2D pc, float *patch_tmp, int level)
{
  const float u_ref = pc[0]; // 获取参考点的X坐标
  const float v_ref = pc[1]; // 获取参考点的Y坐标
  const int scale = (1 << level); // 计算缩放因子
  const int u_ref_i = floorf(pc[0] / scale) * scale; // 计算参考点的整数X坐标
  const int v_ref_i = floorf(pc[1] / scale) * scale; // 计算参考点的整数Y坐标
  const float subpix_u_ref = (u_ref - u_ref_i) / scale; // 计算X子像素偏移
  const float subpix_v_ref = (v_ref - v_ref_i) / scale; // 计算Y子像素偏移
  const float w_ref_tl = (1.0 - subpix_u_ref) * (1.0 - subpix_v_ref); // 左上角权重
  const float w_ref_tr = subpix_u_ref * (1.0 - subpix_v_ref); // 右上角权重
  const float w_ref_bl = (1.0 - subpix_u_ref) * subpix_v_ref; // 左下角权重
  const float w_ref_br = subpix_u_ref * subpix_v_ref; // 右下角权重
  for (int x = 0; x < patch_size; x++)
  {
    uint8_t *img_ptr = (uint8_t *)img.data + (v_ref_i - patch_size_half * scale + x * scale) * width + (u_ref_i - patch_size_half * scale); // 计算图像指针
    for (int y = 0; y < patch_size; y++, img_ptr += scale) // 遍历补丁
    {
      patch_tmp[patch_size_total * level + x * patch_size + y] =
          w_ref_tl * img_ptr[0] + w_ref_tr * img_ptr[scale] + w_ref_bl * img_ptr[scale * width] + w_ref_br * img_ptr[scale * width + scale]; // 计算补丁像素值
    }
  }
}

// 将点插入到体素地图中
void VIOManager::insertPointIntoVoxelMap(VisualPoint *pt_new)
{
  V3D pt_w(pt_new->pos_[0], pt_new->pos_[1], pt_new->pos_[2]); // 获取点的世界坐标
  double voxel_size = 0.5; // 定义体素大小
  float loc_xyz[3]; // 存储体素位置
  for (int j = 0; j < 3; j++)
  {
    loc_xyz[j] = pt_w[j] / voxel_size; // 计算体素位置
    if (loc_xyz[j] < 0) { loc_xyz[j] -= 1.0; } // 处理负值
  }
  VOXEL_LOCATION position((int64_t)loc_xyz[0], (int64_t)loc_xyz[1], (int64_t)loc_xyz[2]); // 创建体素位置
  auto iter = feat_map.find(position); // 查找体素位置
  if (iter != feat_map.end()) // 如果体素位置已存在
  {
    iter->second->voxel_points.push_back(pt_new); // 将新点添加到体素点集合
    iter->second->count++; // 更新体素点计数
  }
  else // 如果体素位置不存在
  {
    VOXEL_POINTS *ot = new VOXEL_POINTS(0); // 创建新的体素点集合
    ot->voxel_points.push_back(pt_new); // 添加新点
    feat_map[position] = ot; // 将体素位置与点集合关联
  }
}

// 获取变换矩阵的仿射单应矩阵
void VIOManager::getWarpMatrixAffineHomography(const vk::AbstractCamera &cam, const V2D &px_ref, const V3D &xyz_ref, const V3D &normal_ref,
                                                  const SE3 &T_cur_ref, const int level_ref, Matrix2d &A_cur_ref)
{
  // 创建单应矩阵
  const V3D t = T_cur_ref.inverse().translation(); // 获取平移部分
  const Eigen::Matrix3d H_cur_ref =
      T_cur_ref.rotation_matrix() * (normal_ref.dot(xyz_ref) * Eigen::Matrix3d::Identity() - t * normal_ref.transpose()); // 计算单应矩阵
  // 计算仿射变换矩阵A_ref_cur
  const int kHalfPatchSize = 4; // 定义补丁大小
  V3D f_du_ref(cam.cam2world(px_ref + Eigen::Vector2d(kHalfPatchSize, 0) * (1 << level_ref))); // 计算右侧补丁的世界坐标
  V3D f_dv_ref(cam.cam2world(px_ref + Eigen::Vector2d(0, kHalfPatchSize) * (1 << level_ref))); // 计算下侧补丁的世界坐标
  const V3D f_cur(H_cur_ref * xyz_ref); // 计算当前点的世界坐标
  const V3D f_du_cur = H_cur_ref * f_du_ref; // 计算右侧补丁的当前坐标
  const V3D f_dv_cur = H_cur_ref * f_dv_ref; // 计算下侧补丁的当前坐标
  V2D px_cur(cam.world2cam(f_cur)); // 转换当前点到相机坐标
  V2D px_du_cur(cam.world2cam(f_du_cur)); // 转换右侧补丁到相机坐标
  V2D px_dv_cur(cam.world2cam(f_dv_cur)); // 转换下侧补丁到相机坐标
  A_cur_ref.col(0) = (px_du_cur - px_cur) / kHalfPatchSize; // 计算仿射矩阵的第一列
  A_cur_ref.col(1) = (px_dv_cur - px_cur) / kHalfPatchSize; // 计算仿射矩阵的第二列
}

// 获取变换矩阵的仿射矩阵
void VIOManager::getWarpMatrixAffine(const vk::AbstractCamera &cam, const Vector2d &px_ref, const Vector3d &f_ref, const double depth_ref,
                                        const SE3 &T_cur_ref, const int level_ref, const int pyramid_level, const int halfpatch_size,
                                        Matrix2d &A_cur_ref)
{
  // 计算仿射变换矩阵A_ref_cur
  const Vector3d xyz_ref(f_ref * depth_ref); // 计算参考点的世界坐标
  Vector3d xyz_du_ref(cam.cam2world(px_ref + Vector2d(halfpatch_size, 0) * (1 << level_ref) * (1 << pyramid_level))); // 计算右侧补丁的世界坐标
  Vector3d xyz_dv_ref(cam.cam2world(px_ref + Vector2d(0, halfpatch_size) * (1 << level_ref) * (1 << pyramid_level))); // 计算下侧补丁的世界坐标
  xyz_du_ref *= xyz_ref[2] / xyz_du_ref[2]; // 归一化右侧补丁的深度
  xyz_dv_ref *= xyz_ref[2] / xyz_dv_ref[2]; // 归一化下侧补丁的深度
  const Vector2d px_cur(cam.world2cam(T_cur_ref * (xyz_ref))); // 转换当前点到相机坐标
  const Vector2d px_du(cam.world2cam(T_cur_ref * (xyz_du_ref))); // 转换右侧补丁到相机坐标
  const Vector2d px_dv(cam.world2cam(T_cur_ref * (xyz_dv_ref))); // 转换下侧补丁到相机坐标
  A_cur_ref.col(0) = (px_du - px_cur) / halfpatch_size; // 计算仿射矩阵的第一列
  A_cur_ref.col(1) = (px_dv - px_cur) / halfpatch_size; // 计算仿射矩阵的第二列
}

// 进行仿射变换
void VIOManager::warpAffine(const Matrix2d &A_cur_ref, const cv::Mat &img_ref, const Vector2d &px_ref, const int level_ref, const int search_level,
                               const int pyramid_level, const int halfpatch_size, float *patch)
{
  const int patch_size = halfpatch_size * 2; // 计算补丁的总大小
  const Matrix2f A_ref_cur = A_cur_ref.inverse().cast<float>(); // 计算逆仿射矩阵
  if (isnan(A_ref_cur(0, 0))) // 检查是否为NaN
  {
    printf("Affine warp is NaN, probably camera has no translation\n"); // 打印警告
    return; // 返回
  }

  float *patch_ptr = patch; // 指向补丁的指针
  for (int y = 0; y < patch_size; ++y) // 遍历补丁
  {
    for (int x = 0; x < patch_size; ++x) // 遍历补丁
    {
      Vector2f px_patch(x - halfpatch_size, y - halfpatch_size); // 计算补丁中的点
      px_patch *= (1 << search_level); // 根据搜索级别缩放
      px_patch *= (1 << pyramid_level); // 根据金字塔级别缩放
      const Vector2f px(A_ref_cur * px_patch + px_ref.cast<float>()); // 计算变换后的像素坐标
      if (px[0] < 0 || px[1] < 0 || px[0] >= img_ref.cols - 1 || px[1] >= img_ref.rows - 1) // 检查是否越界
        patch_ptr[patch_size_total * pyramid_level + y * patch_size + x] = 0; // 越界则设置为0
      else
        patch_ptr[patch_size_total * pyramid_level + y * patch_size + x] = (float)vk::interpolateMat_8u(img_ref, px[0], px[1]); // 进行双线性插值
    }
  }
}

// 获取最佳搜索级别
int VIOManager::getBestSearchLevel(const Matrix2d &A_cur_ref, const int max_level)
{
  // 计算其他图像中的补丁级别
  int search_level = 0; // 当前搜索级别
  double D = A_cur_ref.determinant(); // 计算行列式
  while (D > 3.0 && search_level < max_level) // 检查行列式和最大级别
  {
    search_level += 1; // 增加搜索级别
    D *= 0.25; // 更新行列式
  }
  return search_level; // 返回搜索级别
}

// 计算归一化互相关
double VIOManager::calculateNCC(float *ref_patch, float *cur_patch, int patch_size)
{
  double sum_ref = std::accumulate(ref_patch, ref_patch + patch_size, 0.0); // 计算参考补丁的总和
  double mean_ref = sum_ref / patch_size; // 计算参考补丁的均值

  double sum_cur = std::accumulate(cur_patch, cur_patch + patch_size, 0.0); // 计算当前补丁的总和
  double mean_curr = sum_cur / patch_size; // 计算当前补丁的均值

  double numerator = 0, demoniator1 = 0, demoniator2 = 0; // 初始化分子和分母
  for (int i = 0; i < patch_size; i++)
  {
    double n = (ref_patch[i] - mean_ref) * (cur_patch[i] - mean_curr); // 计算分子
    numerator += n; // 累加分子
    demoniator1 += (ref_patch[i] - mean_ref) * (ref_patch[i] - mean_ref); // 计算参考补丁的分母
    demoniator2 += (cur_patch[i] - mean_curr) * (cur_patch[i] - mean_curr); // 计算当前补丁的分母
  }
  return numerator / sqrt(demoniator1 * demoniator2 + 1e-10); // 返回归一化互相关
}
// 从视觉稀疏地图中检索点
void VIOManager::retrieveFromVisualSparseMap(cv::Mat img, vector<pointWithVar> &pg, const unordered_map<VOXEL_LOCATION, VoxelOctoTree *> &plane_map)
{
  if (feat_map.size() <= 0) return; // 如果特征地图为空，则返回
  double ts0 = omp_get_wtime(); // 记录开始时间

  // pg_down->reserve(feat_map.size());
  // downSizeFilter.setInputCloud(pg);
  // downSizeFilter.filter(*pg_down); // 对输入点云进行下采样（注释掉了）

  // resetRvizDisplay(); // 重置可视化显示（注释掉了）
  visual_submap->reset(); // 重置视觉子图

  // 控制是否包括之前帧的视觉子图
  sub_feat_map.clear(); // 清空子特征地图

  float voxel_size = 0.5; // 定义体素大小

  if (!normal_en) warp_map.clear(); // 如果不启用法线，则清空变换地图

  cv::Mat depth_img = cv::Mat::zeros(height, width, CV_32FC1); // 创建深度图像
  float *it = (float *)depth_img.data; // 获取深度图像的数据指针

  int loc_xyz[3]; // 存储体素坐标

  // printf("A0. initial depthmap: %.6lf \n", omp_get_wtime() - ts0);
  // double ts1 = omp_get_wtime();

  // printf("pg size: %zu \n", pg.size());

  for (int i = 0; i < pg.size(); i++) // 遍历输入点
  {
    V3D pt_w = pg[i].point_w; // 获取点的世界坐标

    for (int j = 0; j < 3; j++)
    {
      loc_xyz[j] = floor(pt_w[j] / voxel_size); // 计算体素坐标
      if (loc_xyz[j] < 0) { loc_xyz[j] -= 1.0; } // 处理负值
    }
    VOXEL_LOCATION position(loc_xyz[0], loc_xyz[1], loc_xyz[2]); // 创建体素位置

    auto iter = sub_feat_map.find(position); // 查找子特征地图中的位置
    if (iter == sub_feat_map.end()) { sub_feat_map[position] = 0; } // 如果位置不存在，则插入
    else { iter->second = 0; } // 如果位置存在，则重置

    V3D pt_c(new_frame_->w2f(pt_w)); // 将世界坐标转换为相机坐标

    if (pt_c[2] > 0) // 如果Z坐标大于0
    {
      V2D px; // 存储图像坐标
      px = new_frame_->cam_->world2cam(pt_c); // 将相机坐标转换为图像坐标

      if (new_frame_->cam_->isInFrame(px.cast<int>(), border)) // 检查点是否在相机视野内
      {
        float depth = pt_c[2]; // 获取深度值
        int col = int(px[0]); // 获取列索引
        int row = int(px[1]); // 获取行索引
        it[width * row + col] = depth; // 将深度值写入深度图像
      }
    }
    // t_depth += omp_get_wtime()-t2;
  }

  // printf("A1: %.6lf \n", omp_get_wtime() - ts1);
  vector<VOXEL_LOCATION> DeleteKeyList; // 存储需要删除的键值

  for (auto &iter : sub_feat_map) // 遍历子特征地图
  {
    VOXEL_LOCATION position = iter.first; // 获取体素位置

    auto corre_voxel = feat_map.find(position); // 查找特征地图中的对应体素
    if (corre_voxel != feat_map.end()) // 如果存在对应体素
    {
      bool voxel_in_fov = false; // 初始化视野标志
      std::vector<VisualPoint *> &voxel_points = corre_voxel->second->voxel_points; // 获取体素点集合
      int voxel_num = voxel_points.size(); // 获取体素点数量

      for (int i = 0; i < voxel_num; i++) // 遍历体素点
      {
        VisualPoint *pt = voxel_points[i]; // 获取体素点
        if (pt == nullptr) continue; // 如果点为空，则继续
        if (pt->obs_.size() == 0) continue; // 如果点没有观察数据，则继续

        V3D norm_vec(new_frame_->T_f_w_.rotation_matrix() * pt->normal_); // 计算法向量
        V3D dir(new_frame_->T_f_w_ * pt->pos_); // 计算点的方向
        if (dir[2] < 0) continue; // 如果Z坐标小于0，则继续

        V2D pc(new_frame_->w2c(pt->pos_)); // 将点的世界坐标转换为相机坐标
        if (new_frame_->cam_->isInFrame(pc.cast<int>(), border)) // 检查点是否在相机视野内
        {
          voxel_in_fov = true; // 设置视野标志
          int index = static_cast<int>(pc[1] / grid_size) * grid_n_width + static_cast<int>(pc[0] / grid_size); // 计算网格索引
          grid_num[index] = TYPE_MAP; // 更新网格状态
          Vector3d obs_vec(new_frame_->pos() - pt->pos_); // 计算观察向量
          float cur_dist = obs_vec.norm(); // 计算距离
          if (cur_dist <= map_dist[index]) // 如果当前距离小于等于记录的距离
          {
            map_dist[index] = cur_dist; // 更新距离
            retrieve_voxel_points[index] = pt; // 更新对应网格的体素点
          }
        }
      }
      if (!voxel_in_fov) { DeleteKeyList.push_back(position); } // 如果没有在视野内，则添加到删除列表
    }
  }

  // 光线投射模块
  if (raycast_en)
  {
    for (int i = 0; i < length; i++) // 遍历所有网格
    {
      if (grid_num[i] == TYPE_MAP || border_flag[i] == 1) continue; // 如果网格已标记为地图或是边界，则继续

      for (const auto &it : rays_with_sample_points[i]) // 遍历样本点
      {
        V3D sample_point_w = new_frame_->f2w(it); // 将样本点转换为世界坐标

        for (int j = 0; j < 3; j++)
        {
          loc_xyz[j] = floor(sample_point_w[j] / voxel_size); // 计算样本点的体素坐标
          if (loc_xyz[j] < 0) { loc_xyz[j] -= 1.0; } // 处理负值
        }

        VOXEL_LOCATION sample_pos(loc_xyz[0], loc_xyz[1], loc_xyz[2]); // 创建样本点位置

        auto corre_sub_feat_map = sub_feat_map.find(sample_pos); // 查找子特征地图中的样本点
        if (corre_sub_feat_map != sub_feat_map.end()) break; // 如果找到，则跳出循环

        auto corre_feat_map = feat_map.find(sample_pos); // 查找特征地图中的样本点
        if (corre_feat_map != feat_map.end()) // 如果存在
        {
          bool voxel_in_fov = false; // 初始化视野标志

          std::vector<VisualPoint *> &voxel_points = corre_feat_map->second->voxel_points; // 获取体素点集合
          int voxel_num = voxel_points.size(); // 获取体素点数量
          if (voxel_num == 0) continue; // 如果体素点数量为0则继续

          for (int j = 0; j < voxel_num; j++) // 遍历体素点
          {
            VisualPoint *pt = voxel_points[j]; // 获取体素点

            if (pt == nullptr) continue; // 如果点为空，则继续
            if (pt->obs_.size() == 0) continue; // 如果点没有观察数据，则继续

            V3D norm_vec(new_frame_->T_f_w_.rotation_matrix() * pt->normal_); // 计算法向量
            V3D dir(new_frame_->T_f_w_ * pt->pos_); // 计算点的方向
            if (dir[2] < 0) continue; // 如果Z坐标小于0，则继续
            dir.normalize(); // 归一化方向向量

            V2D pc(new_frame_->w2c(pt->pos_)); // 将点的世界坐标转换为相机坐标

            if (new_frame_->cam_->isInFrame(pc.cast<int>(), border)) // 检查点是否在相机视野内
            {
              voxel_in_fov = true; // 设置视野标志
              int index = static_cast<int>(pc[1] / grid_size) * grid_n_width + static_cast<int>(pc[0] / grid_size); // 计算网格索引
              grid_num[index] = TYPE_MAP; // 更新网格状态
              Vector3d obs_vec(new_frame_->pos() - pt->pos_); // 计算观察向量

              float cur_dist = obs_vec.norm(); // 计算距离

              if (cur_dist <= map_dist[index]) // 如果当前距离小于等于记录的距离
              {
                map_dist[index] = cur_dist; // 更新距离
                retrieve_voxel_points[index] = pt; // 更新对应网格的体素点
              }
            }
          }
          if (voxel_in_fov) sub_feat_map[sample_pos] = 0; // 如果体素在视野内，则更新子特征地图
          break; // 跳出循环
        }
        else // 如果特征地图中没有找到样本点
        {
          VOXEL_LOCATION sample_pos(loc_xyz[0], loc_xyz[1], loc_xyz[2]); // 创建样本点位置
          auto iter = plane_map.find(sample_pos); // 查找平面地图中的样本点
          if (iter != plane_map.end()) // 如果找到
          {
            VoxelOctoTree *current_octo;
            current_octo = iter->second->find_correspond(sample_point_w); // 找到对应的八叉树
            if (current_octo->plane_ptr_->is_plane_) // 如果是平面
            {
              pointWithVar plane_center; // 创建平面中心点
              VoxelPlane &plane = *current_octo->plane_ptr_; // 获取平面
              plane_center.point_w = plane.center_; // 设置平面中心点的世界坐标
              plane_center.normal = plane.normal_; // 设置平面法线
              visual_submap->add_from_voxel_map.push_back(plane_center); // 将平面中心点添加到视觉子图中
              break; // 跳出循环
            }
          }
        }
      }
      // if(add_sample) sample_points.push_back(sample_points_temp);
    }
  }

  for (auto &key : DeleteKeyList) // 遍历需要删除的键值
  {
    sub_feat_map.erase(key); // 从子特征地图中删除
  }

  for (int i = 0; i < length; i++) // 遍历所有网格
  {
    if (grid_num[i] == TYPE_MAP) // 如果网格标记为地图
    {
      VisualPoint *pt = retrieve_voxel_points[i]; // 获取对应的体素点

      V2D pc(new_frame_->w2c(pt->pos_)); // 将点的世界坐标转换为相机坐标

      V3D pt_cam(new_frame_->w2f(pt->pos_)); // 将点的世界坐标转换为相机坐标
      bool depth_continous = false; // 初始化深度连续性标志
      for (int u = -patch_size_half; u <= patch_size_half; u++) // 遍历补丁的范围
      {
        for (int v = -patch_size_half; v <= patch_size_half; v++)
        {
          if (u == 0 && v == 0) continue; // 跳过中心点

          float depth = it[width * (v + int(pc[1])) + u + int(pc[0])]; // 获取深度值

          if (depth == 0.) continue; // 如果深度值为0，则继续

          double delta_dist = abs(pt_cam[2] - depth); // 计算深度差

          if (delta_dist > 0.5) // 如果深度差大于0.5
          {
            depth_continous = true; // 设置深度连续性标志
            break; // 跳出循环
          }
        }
        if (depth_continous) break; // 如果已发现深度不连续，跳出循环
      }
      if (depth_continous) continue; // 如果深度不连续，则继续

      Feature *ref_ftr; // 定义参考特征
      std::vector<float> patch_wrap(warp_len); // 创建补丁包裹数组

      int search_level; // 搜索层级
      Matrix2d A_cur_ref_zero; // 当前参考的仿射矩阵

      if (!pt->is_normal_initialized_) continue; // 如果法线未初始化，则继续

      if (normal_en) // 如果启用法线
      {
        float phtometric_errors_min = std::numeric_limits<float>::max(); // 初始化光度误差最小值

        if (pt->obs_.size() == 1) // 如果观察数据只有一个
        {
          ref_ftr = *pt->obs_.begin(); // 设置参考特征
          pt->ref_patch = ref_ftr; // 设置参考补丁
          pt->has_ref_patch_ = true; // 标记为有参考补丁
        }
        else if (!pt->has_ref_patch_) // 如果没有参考补丁
        {
          for (auto it = pt->obs_.begin(), ite = pt->obs_.end(); it != ite; ++it) // 遍历观察数据
          {
            Feature *ref_patch_temp = *it; // 获取临时参考补丁
            float *patch_temp = ref_patch_temp->patch_; // 获取参考补丁数据
            float phtometric_errors = 0.0; // 初始化光度误差
            int count = 0; // 初始化计数器
            for (auto itm = pt->obs_.begin(), itme = pt->obs_.end(); itm != itme; ++itm) // 遍历观察数据
            {
              if ((*itm)->id_ == ref_patch_temp->id_) continue; // 跳过相同ID的补丁
              float *patch_cache = (*itm)->patch_; // 获取缓存补丁数据

              for (int ind = 0; ind < patch_size_total; ind++) // 计算光度误差
              {
                phtometric_errors += (patch_temp[ind] - patch_cache[ind]) * (patch_temp[ind] - patch_cache[ind]);
              }
              count++; // 增加计数
            }
            phtometric_errors = phtometric_errors / count; // 计算平均光度误差
            if (phtometric_errors < phtometric_errors_min) // 如果当前误差小于最小误差
            {
              phtometric_errors_min = phtometric_errors; // 更新最小误差
              ref_ftr = ref_patch_temp; // 更新参考特征
            }
          }
          pt->ref_patch = ref_ftr; // 设置参考补丁
          pt->has_ref_patch_ = true; // 标记为有参考补丁
        }
        else { ref_ftr = pt->ref_patch; } // 如果已有参考补丁，则直接使用
      }
      else // 如果不启用法线
      {
        if (!pt->getCloseViewObs(new_frame_->pos(), ref_ftr, pc)) continue; // 获取最近观察的特征
      }

      if (normal_en) // 如果启用法线
      {
        V3D norm_vec = (ref_ftr->T_f_w_.rotation_matrix() * pt->normal_).normalized(); // 计算法向量并归一化
        
        V3D pf(ref_ftr->T_f_w_ * pt->pos_); // 获取参考特征的世界坐标
        SE3 T_cur_ref = new_frame_->T_f_w_ * ref_ftr->T_f_w_.inverse(); // 计算当前和参考之间的变换

        getWarpMatrixAffineHomography(*cam, ref_ftr->px_, pf, norm_vec, T_cur_ref, 0, A_cur_ref_zero); // 获取仿射变换矩阵
        search_level = getBestSearchLevel(A_cur_ref_zero, 2); // 获取最佳搜索层级
      }
      else // 如果不启用法线
      {
        auto iter_warp = warp_map.find(ref_ftr->id_); // 查找变换地图
        if (iter_warp != warp_map.end()) // 如果找到
        {
          search_level = iter_warp->second->search_level; // 获取搜索层级
          A_cur_ref_zero = iter_warp->second->A_cur_ref; // 获取当前参考的变换矩阵
        }
        else // 如果未找到
        {
          getWarpMatrixAffine(*cam, ref_ftr->px_, ref_ftr->f_, (ref_ftr->pos() - pt->pos_).norm(), new_frame_->T_f_w_ * ref_ftr->T_f_w_.inverse(),
                              ref_ftr->level_, 0, patch_size_half, A_cur_ref_zero); // 获取仿射变换矩阵

          search_level = getBestSearchLevel(A_cur_ref_zero, 2); // 获取最佳搜索层级

          Warp *ot = new Warp(search_level, A_cur_ref_zero); // 创建新的变换对象
          warp_map[ref_ftr->id_] = ot; // 将变换对象添加到变换地图
        }
      }
      // t_4 += omp_get_wtime() - t_1;

      for (int pyramid_level = 0; pyramid_level <= patch_pyrimid_level - 1; pyramid_level++) // 遍历金字塔层级
      {
        warpAffine(A_cur_ref_zero, ref_ftr->img_, ref_ftr->px_, ref_ftr->level_, search_level, pyramid_level, patch_size_half, patch_wrap.data()); // 执行仿射变换
      }

      getImagePatch(img, pc, patch_buffer.data(), 0); // 从图像中获取补丁

      float error = 0.0; // 初始化误差
      for (int ind = 0; ind < patch_size_total; ind++) // 计算误差
      {
        error += (ref_ftr->inv_expo_time_ * patch_wrap[ind] - state->inv_expo_time * patch_buffer[ind]) *
                 (ref_ftr->inv_expo_time_ * patch_wrap[ind] - state->inv_expo_time * patch_buffer[ind]);
      }

      if (ncc_en) // 如果启用NCC
      {
        double ncc = calculateNCC(patch_wrap.data(), patch_buffer.data(), patch_size_total); // 计算归一化互相关
        if (ncc < ncc_thre) // 如果NCC小于阈值
        {
          continue; // 跳过当前点
        }
      }

      if (error > outlier_threshold * patch_size_total) continue; // 如果误差超过阈值，则继续

      visual_submap->voxel_points.push_back(pt); // 将体素点添加到视觉子图中
      visual_submap->propa_errors.push_back(error); // 添加误差到视觉子图
      visual_submap->search_levels.push_back(search_level); // 添加搜索层级到视觉子图
      visual_submap->errors.push_back(error); // 添加误差到视觉子图
      visual_submap->warp_patch.push_back(patch_wrap); // 添加补丁包裹到视觉子图
      visual_submap->inv_expo_list.push_back(ref_ftr->inv_expo_time_); // 添加反曝光时间到视觉子图
    }
  }
  total_points = visual_submap->voxel_points.size(); // 统计视觉子图中的点的数量

  printf("[ VIO ] Retrieve %d points from visual sparse map\n", total_points); // 输出检索的点的数量
}
//计算雅可比矩阵并更新EKF
void VIOManager::computeJacobianAndUpdateEKF(cv::Mat img)
{
  if (total_points == 0) return; // 如果总点数为0，则返回
  
  compute_jacobian_time = update_ekf_time = 0.0; // 初始化计算雅可比时间和更新EKF时间

  for (int level = patch_pyrimid_level - 1; level >= 0; level--) // 从金字塔层级的最高层遍历到最低层
  {
    if (inverse_composition_en) // 如果启用逆组合方法
    {
      has_ref_patch_cache = false; // 清空参考补丁缓存
      updateStateInverse(img, level); // 更新状态（逆组合）
    }
    else
      updateState(img, level); // 更新状态（正常组合）
  }
  state->cov -= G * state->cov; // 更新状态协方差
  updateFrameState(*state); // 更新帧状态
}
//获取最近观察的特征
void VIOManager::generateVisualMapPoints(cv::Mat img, vector<pointWithVar> &pg)
{
  if (pg.size() <= 10) return; // 如果点云大小小于等于10，则返回

  // double t0 = omp_get_wtime(); // 记录时间（可选）
  for (int i = 0; i < pg.size(); i++) // 遍历输入点云
  {
    if (pg[i].normal == V3D(0, 0, 0)) continue; // 如果法向量为零，则跳过

    V3D pt = pg[i].point_w; // 获取点的世界坐标
    V2D pc(new_frame_->w2c(pt)); // 将世界坐标转换为相机坐标

    if (new_frame_->cam_->isInFrame(pc.cast<int>(), border)) // 如果点在相机视野内
    {
      int index = static_cast<int>(pc[1] / grid_size) * grid_n_width + static_cast<int>(pc[0] / grid_size); // 计算网格索引

      if (grid_num[index] != TYPE_MAP) // 如果该网格不属于地图类型
      {
        float cur_value = vk::shiTomasiScore(img, pc[0], pc[1]); // 计算Shi-Tomasi特征值
        // if (cur_value < 5) continue; // 如果特征值小于5，则跳过
        if (cur_value > scan_value[index]) // 如果当前值大于扫描值
        {
          scan_value[index] = cur_value; // 更新扫描值
          append_voxel_points[index] = pg[i]; // 添加体素点
          grid_num[index] = TYPE_POINTCLOUD; // 标记为点云类型
        }
      }
    }
  }

  for (int j = 0; j < visual_submap->add_from_voxel_map.size(); j++) // 遍历视觉子图中的点
  {
    V3D pt = visual_submap->add_from_voxel_map[j].point_w; // 获取点的世界坐标
    V2D pc(new_frame_->w2c(pt)); // 转换为相机坐标

    if (new_frame_->cam_->isInFrame(pc.cast<int>(), border)) // 如果点在相机视野内
    {
      int index = static_cast<int>(pc[1] / grid_size) * grid_n_width + static_cast<int>(pc[0] / grid_size); // 计算网格索引

      if (grid_num[index] != TYPE_MAP) // 如果该网格不属于地图类型
      {
        float cur_value = vk::shiTomasiScore(img, pc[0], pc[1]); // 计算Shi-Tomasi特征值
        if (cur_value > scan_value[index]) // 如果当前值大于扫描值
        {
          scan_value[index] = cur_value; // 更新扫描值
          append_voxel_points[index] = visual_submap->add_from_voxel_map[j]; // 添加体素点
          grid_num[index] = TYPE_POINTCLOUD; // 标记为点云类型
        }
      }
    }
  }

  // double t_b1 = omp_get_wtime() - t0; // 记录时间（可选）
  // t0 = omp_get_wtime(); // 记录时间（可选）

  int add = 0; // 初始化添加计数
  for (int i = 0; i < length; i++) // 遍历所有网格
  {
    if (grid_num[i] == TYPE_POINTCLOUD) // 如果网格标记为点云类型
    {
      pointWithVar pt_var = append_voxel_points[i]; // 获取体素点
      V3D pt = pt_var.point_w; // 获取点的世界坐标

      V3D norm_vec(new_frame_->T_f_w_.rotation_matrix() * pt_var.normal); // 计算法向量
      V3D dir(new_frame_->T_f_w_ * pt); // 计算方向向量
      dir.normalize(); // 归一化方向向量
      double cos_theta = dir.dot(norm_vec); // 计算方向与法向量的夹角余弦
      // if(std::fabs(cos_theta)<0.34) continue; // 70度
      V2D pc(new_frame_->w2c(pt)); // 转换为相机坐标

      float *patch = new float[patch_size_total]; // 创建补丁数组
      getImagePatch(img, pc, patch, 0); // 从图像中获取补丁

      VisualPoint *pt_new = new VisualPoint(pt); // 创建新的视觉点

      Vector3d f = cam->cam2world(pc); // 将相机坐标转换为世界坐标
      Feature *ftr_new = new Feature(pt_new, patch, pc, f, new_frame_->T_f_w_, 0); // 创建新的特征
      ftr_new->img_ = img; // 设置特征图像
      ftr_new->id_ = new_frame_->id_; // 设置特征ID
      ftr_new->inv_expo_time_ = state->inv_expo_time; // 设置反曝光时间

      pt_new->addFrameRef(ftr_new); // 将特征添加到视觉点中
      pt_new->covariance_ = pt_var.var; // 设置协方差
      pt_new->is_normal_initialized_ = true; // 标记法线已初始化

      if (cos_theta < 0) { pt_new->normal_ = -pt_var.normal; } // 如果夹角余弦为负，则取法线反向
      else { pt_new->normal_ = pt_var.normal; } // 否则直接使用法线
      
      pt_new->previous_normal_ = pt_new->normal_; // 更新之前的法线

      insertPointIntoVoxelMap(pt_new); // 将新点插入到体素地图
      add += 1; // 增加添加计数
      // map_cur_frame.push_back(pt_new); // 可选：将当前帧的点添加到列表
    }
  }

  // double t_b2 = omp_get_wtime() - t0; // 记录时间（可选）

  printf("[ VIO ] Append %d new visual map points\n", add); // 输出添加的新视觉点数量
  // printf("pg.size: %d \n", pg.size()); // 可选：输出点云大小
  // printf("B1. : %.6lf \n", t_b1); // 可选：输出时间
  // printf("B2. : %.6lf \n", t_b2); // 可选：输出时间
}
//更新视觉子图中的点
void VIOManager::updateVisualMapPoints(cv::Mat img)
{
  if (total_points == 0) return; // 如果总点数为0，则返回

  int update_num = 0; // 初始化更新计数
  SE3 pose_cur = new_frame_->T_f_w_; // 获取当前帧的位姿
  for (int i = 0; i < total_points; i++) // 遍历视觉子图中的所有点
  {
    VisualPoint *pt = visual_submap->voxel_points[i]; // 获取视觉点
    if (pt == nullptr) continue; // 如果视觉点为空，则跳过
    if (pt->is_converged_) // 如果点已收敛
    { 
      pt->deleteNonRefPatchFeatures(); // 删除非参考补丁特征
      continue; // 跳过
    }

    V2D pc(new_frame_->w2c(pt->pos_)); // 将点的世界坐标转换为相机坐标
    bool add_flag = false; // 初始化添加标志
    
    float *patch_temp = new float[patch_size_total]; // 创建补丁数组
    getImagePatch(img, pc, patch_temp, 0); // 从图像中获取补丁
    // TODO: condition: distance and view_angle
    // Step 1: time
    Feature *last_feature = pt->obs_.back(); // 获取最后一个特征
    // if(new_frame_->id_ >= last_feature->id_ + 10) add_flag = true; // 10

    // Step 2: delta_pose
    SE3 pose_ref = last_feature->T_f_w_; // 获取参考帧的位姿
    SE3 delta_pose = pose_ref * pose_cur.inverse(); // 计算当前帧与参考帧之间的位姿变化
    double delta_p = delta_pose.translation().norm(); // 计算平移距离
    double delta_theta = (delta_pose.rotation_matrix().trace() > 3.0 - 1e-6) ? 0.0 : std::acos(0.5 * (delta_pose.rotation_matrix().trace() - 1)); // 计算旋转角度
    if (delta_p > 0.5 || delta_theta > 0.3) add_flag = true; // 如果平移距离大于0.5或旋转角度大于0.3，则设置添加标志

    // Step 3: pixel distance
    Vector2d last_px = last_feature->px_; // 获取最后特征的像素坐标
    double pixel_dist = (pc - last_px).norm(); // 计算当前像素与最后像素的距离
    if (pixel_dist > 40) add_flag = true; // 如果像素距离大于40，则设置添加标志

    // Maintain the size of 3D point observation features.
    if (pt->obs_.size() >= 30) // 如果观察特征数量达到30
    {
      Feature *ref_ftr; // 定义参考特征
      pt->findMinScoreFeature(new_frame_->pos(), ref_ftr); // 查找最小评分特征
      pt->deleteFeatureRef(ref_ftr); // 删除该特征
      // cout<<"pt->obs_.size() exceed 20 !!!!!!"<<endl; // 可选：输出信息
    }
    if (add_flag) // 如果需要添加
    {
      update_num += 1; // 增加更新计数
      update_flag[i] = 1; // 设置更新标志
      Vector3d f = cam->cam2world(pc); // 将相机坐标转换为世界坐标
      Feature *ftr_new = new Feature(pt, patch_temp, pc, f, new_frame_->T_f_w_, visual_submap->search_levels[i]); // 创建新的特征
      ftr_new->img_ = img; // 设置特征图像
      ftr_new->id_ = new_frame_->id_; // 设置特征ID
      ftr_new->inv_expo_time_ = state->inv_expo_time; // 设置反曝光时间
      pt->addFrameRef(ftr_new); // 将特征添加到视觉点中
    }
  }
  printf("[ VIO ] Update %d points in visual submap\n", update_num); // 输出更新的点数量
}
//更新参考补丁
void VIOManager::updateReferencePatch(const unordered_map<VOXEL_LOCATION, VoxelOctoTree *> &plane_map)
{
  if (total_points == 0) return; // 如果总点数为0，则返回

  for (int i = 0; i < visual_submap->voxel_points.size(); i++) // 遍历视觉子图中的所有点
  {
    VisualPoint *pt = visual_submap->voxel_points[i]; // 获取视觉点

    if (!pt->is_normal_initialized_) continue; // 如果法线未初始化，则跳过
    if (pt->is_converged_) continue; // 如果点已收敛，则跳过
    if (pt->obs_.size() <= 5) continue; // 如果观察特征数量小于等于5，则跳过
    if (update_flag[i] == 0) continue; // 如果更新标志为0，则跳过

    const V3D &p_w = pt->pos_; // 获取点的世界坐标
    float loc_xyz[3]; // 定义位置数组
    for (int j = 0; j < 3; j++)
    {
      loc_xyz[j] = p_w[j] / 0.5; // 计算体素位置
      if (loc_xyz[j] < 0) { loc_xyz[j] -= 1.0; } // 如果位置小于0，则减1
    }
    VOXEL_LOCATION position((int64_t)loc_xyz[0], (int64_t)loc_xyz[1], (int64_t)loc_xyz[2]); // 创建体素位置
    auto iter = plane_map.find(position); // 查找平面地图中的位置
    if (iter != plane_map.end()) // 如果找到
    {
      VoxelOctoTree *current_octo; // 定义当前八叉树
      current_octo = iter->second->find_correspond(p_w); // 找到对应的八叉树
      if (current_octo->plane_ptr_->is_plane_) // 如果是平面
      {
        VoxelPlane &plane = *current_octo->plane_ptr_; // 获取平面
        float dis_to_plane = plane.normal_(0) * p_w(0) + plane.normal_(1) * p_w(1) + plane.normal_(2) * p_w(2) + plane.d_; // 计算到平面的距离
        float dis_to_plane_abs = fabs(dis_to_plane); // 计算距离的绝对值
        float dis_to_center = (plane.center_(0) - p_w(0)) * (plane.center_(0) - p_w(0)) +
                              (plane.center_(1) - p_w(1)) * (plane.center_(1) - p_w(1)) + (plane.center_(2) - p_w(2)) * (plane.center_(2) - p_w(2)); // 计算到平面中心的距离
        float range_dis = sqrt(dis_to_center - dis_to_plane * dis_to_plane); // 计算范围距离
        if (range_dis <= 3 * plane.radius_) // 如果范围距离小于等于平面半径的三倍
        {
          Eigen::Matrix<double, 1, 6> J_nq; // 定义雅可比矩阵
          J_nq.block<1, 3>(0, 0) = p_w - plane.center_; // 设置雅可比矩阵的前3个元素为点到平面中心的向量
          J_nq.block<1, 3>(0, 3) = -plane.normal_; // 设置雅可比矩阵的后3个元素为平面法向量的负值
          double sigma_l = J_nq * plane.plane_var_ * J_nq.transpose(); // 计算方差
          sigma_l += plane.normal_.transpose() * pt->covariance_ * plane.normal_; // 更新方差

          if (dis_to_plane_abs < 3 * sqrt(sigma_l)) // 如果到平面的距离小于方差的三倍
          {
            // V3D norm_vec(new_frame_->T_f_w_.rotation_matrix() * plane.normal_);
            // V3D pf(new_frame_->T_f_w_ * pt->pos_);
            // V3D pf_ref(pt->ref_patch->T_f_w_ * pt->pos_);
            // V3D norm_vec_ref(pt->ref_patch->T_f_w_.rotation_matrix() *
            // plane.normal); double cos_ref = pf_ref.dot(norm_vec_ref);
            
            if (pt->previous_normal_.dot(plane.normal_) < 0) { pt->normal_ = -plane.normal_; } // 如果之前的法线与当前法线夹角大于90度，则取反
            else { pt->normal_ = plane.normal_; } // 否则直接使用当前法线

            double normal_update = (pt->normal_ - pt->previous_normal_).norm(); // 计算法线更新的幅度

            pt->previous_normal_ = pt->normal_; // 更新之前的法线

            if (normal_update < 0.0001 && pt->obs_.size() > 10) // 如果法线更新幅度小于0.0001且观察特征数量大于10
            {
              pt->is_converged_ = true; // 标记为收敛
              // visual_converged_point.push_back(pt); // 可选：将收敛的点添加到列表
            }
          }
        }
      }
    }

    float score_max = -1000.; // 初始化最大评分
    for (auto it = pt->obs_.begin(), ite = pt->obs_.end(); it != ite; ++it) // 遍历观察特征
    {
      Feature *ref_patch_temp = *it; // 获取参考补丁
      float *patch_temp = ref_patch_temp->patch_; // 获取补丁数据
      float NCC_up = 0.0; // 初始化NCC分子
      float NCC_down1 = 0.0; // 初始化NCC分母1
      float NCC_down2 = 0.0; // 初始化NCC分母2
      float NCC = 0.0; // 初始化NCC
      float score = 0.0; // 初始化评分
      int count = 0; // 初始化计数器

      V3D pf = ref_patch_temp->T_f_w_ * pt->pos_; // 获取点的世界坐标
      V3D norm_vec = ref_patch_temp->T_f_w_.rotation_matrix() * pt->normal_; // 获取法向量
      pf.normalize(); // 归一化世界坐标
      double cos_angle = pf.dot(norm_vec); // 计算法向量与点坐标的夹角余弦
      // if(fabs(cos_angle) < 0.86) continue; // 20度

      float ref_mean; // 定义参考补丁均值
      if (abs(ref_patch_temp->mean_) < 1e-6) // 如果均值小于1e-6
      {
        float ref_sum = std::accumulate(patch_temp, patch_temp + patch_size_total, 0.0); // 计算补丁的和
        ref_mean = ref_sum / patch_size_total; // 计算均值
        ref_patch_temp->mean_ = ref_mean; // 更新均值
      }

      for (auto itm = pt->obs_.begin(), itme = pt->obs_.end(); itm != itme; ++itm) // 遍历观察特征
      {
        if ((*itm)->id_ == ref_patch_temp->id_) continue; // 跳过相同ID的补丁
        float *patch_cache = (*itm)->patch_; // 获取缓存补丁数据

        float other_mean; // 定义其他补丁均值
        if (abs((*itm)->mean_) < 1e-6) // 如果均值小于1e-6
        {
          float other_sum = std::accumulate(patch_cache, patch_cache + patch_size_total, 0.0); // 计算补丁的和
          other_mean = other_sum / patch_size_total; // 计算均值
          (*itm)->mean_ = other_mean; // 更新均值
        }

        for (int ind = 0; ind < patch_size_total; ind++) // 计算NCC
        {
          NCC_up += (patch_temp[ind] - ref_mean) * (patch_cache[ind] - other_mean); // 更新NCC分子
          NCC_down1 += (patch_temp[ind] - ref_mean) * (patch_temp[ind] - ref_mean); // 更新NCC分母1
          NCC_down2 += (patch_cache[ind] - other_mean) * (patch_cache[ind] - other_mean); // 更新NCC分母2
        }
        NCC += fabs(NCC_up / sqrt(NCC_down1 * NCC_down2)); // 更新NCC
        count++; // 增加计数
      }

      NCC = NCC / count; // 计算平均NCC

      score = NCC + cos_angle; // 计算评分

      ref_patch_temp->score_ = score; // 更新参考补丁的评分

      if (score > score_max) // 如果当前评分大于最大评分
      {
        score_max = score; // 更新最大评分
        pt->ref_patch = ref_patch_temp; // 更新参考补丁
        pt->has_ref_patch_ = true; // 标记为有参考补丁
      }
    }
  }
}
//投影参考补丁到当前帧
void VIOManager::projectPatchFromRefToCur(const unordered_map<VOXEL_LOCATION, VoxelOctoTree *> &plane_map)
{
  if (total_points == 0) return; // 如果总点数为0，则返回
  // if(new_frame_->id_ != 2) return; //124

  int patch_size = 25; // 定义补丁大小
  string dir = string(ROOT_DIR) + "Log/ref_cur_combine/"; // 定义输出目录

  cv::Mat result = cv::Mat::zeros(height, width, CV_8UC1); // 初始化结果图像
  cv::Mat result_normal = cv::Mat::zeros(height, width, CV_8UC1); // 初始化法线图像
  cv::Mat result_dense = cv::Mat::zeros(height, width, CV_8UC1); // 初始化稠密图像

  cv::Mat img_photometric_error = new_frame_->img_.clone(); // 克隆当前帧图像

  uchar *it = (uchar *)result.data; // 获取结果图像数据指针
  uchar *it_normal = (uchar *)result_normal.data; // 获取法线图像数据指针
  uchar *it_dense = (uchar *)result_dense.data; // 获取稠密图像数据指针

  struct pixel_member // 像素成员结构
  {
    Vector2f pixel_pos; // 像素位置
    uint8_t pixel_value; // 像素值
  };

  int num = 0; // 初始化计数器
  for (int i = 0; i < visual_submap->voxel_points.size(); i++) // 遍历视觉子图中的所有点
  {
    VisualPoint *pt = visual_submap->voxel_points[i]; // 获取视觉点

    if (pt->is_normal_initialized_) // 如果法线已初始化
    {
      Feature *ref_ftr; // 定义参考特征
      ref_ftr = pt->ref_patch; // 获取参考补丁
      // Feature* ref_ftr;
      V2D pc(new_frame_->w2c(pt->pos_)); // 将点的世界坐标转换为相机坐标
      V2D pc_prior(new_frame_->w2c_prior(pt->pos_)); // 获取先前相机坐标

      V3D norm_vec(ref_ftr->T_f_w_.rotation_matrix() * pt->normal_); // 计算法向量
      V3D pf(ref_ftr->T_f_w_ * pt->pos_); // 获取点的世界坐标

      if (pf.dot(norm_vec) < 0) norm_vec = -norm_vec; // 如果法向量与点坐标的点积小于0，则反转法向量

      // norm_vec << norm_vec(1), norm_vec(0), norm_vec(2); // 可选：调整法向量顺序
      cv::Mat img_cur = new_frame_->img_; // 获取当前帧图像
      cv::Mat img_ref = ref_ftr->img_; // 获取参考帧图像

      SE3 T_cur_ref = new_frame_->T_f_w_ * ref_ftr->T_f_w_.inverse(); // 计算当前帧与参考帧之间的变换
      Matrix2d A_cur_ref; // 定义变换矩阵
      getWarpMatrixAffineHomography(*cam, ref_ftr->px_, pf, norm_vec, T_cur_ref, 0, A_cur_ref); // 获取仿射变换矩阵

      // const Matrix2f A_ref_cur = A_cur_ref.inverse().cast<float>(); // 可选：获取逆变换矩阵
      int search_level = getBestSearchLevel(A_cur_ref.inverse(), 2); // 获取最佳搜索层级

      double D = A_cur_ref.determinant(); // 计算变换矩阵的行列式
      if (D > 3) continue; // 如果行列式大于3，则跳过

      num++; // 增加计数

      cv::Mat ref_cur_combine_temp; // 定义合成图像
      int radius = 20; // 定义半径
      cv::hconcat(img_cur, img_ref, ref_cur_combine_temp); // 将当前帧和参考帧图像水平拼接
      cv::cvtColor(ref_cur_combine_temp, ref_cur_combine_temp, CV_GRAY2BGR); // 转换为彩色图像

      getImagePatch(img_cur, pc, patch_buffer.data(), 0); // 从当前图像中获取补丁

      float error_est = 0.0; // 初始化估计误差
      float error_gt = 0.0; // 初始化真实误差

      for (int ind = 0; ind < patch_size_total; ind++) // 计算补丁误差
      {
        error_est += (ref_ftr->inv_expo_time_ * visual_submap->warp_patch[i][ind] - state->inv_expo_time * patch_buffer[ind]) *
                     (ref_ftr->inv_expo_time_ * visual_submap->warp_patch[i][ind] - state->inv_expo_time * patch_buffer[ind]);
      }
      std::string ref_est = "ref_est " + std::to_string(1.0 / ref_ftr->inv_expo_time_); // 参考估计字符串
      std::string cur_est = "cur_est " + std::to_string(1.0 / state->inv_expo_time); // 当前估计字符串
      std::string cur_propa = "cur_gt " + std::to_string(error_gt); // 当前真实值字符串
      std::string cur_optimize = "cur_est " + std::to_string(error_est); // 当前优化估计字符串

      cv::putText(ref_cur_combine_temp, ref_est, cv::Point2f(ref_ftr->px_[0] + img_cur.cols - 40, ref_ftr->px_[1] + 40), cv::FONT_HERSHEY_COMPLEX, 0.4,
                  cv::Scalar(0, 255, 0), 1, 8, 0); // 在图像上绘制参考估计文本

      cv::putText(ref_cur_combine_temp, cur_est, cv::Point2f(pc[0] - 40, pc[1] + 40), cv::FONT_HERSHEY_COMPLEX, 0.4, cv::Scalar(0, 255, 0), 1, 8, 0); // 绘制当前估计文本
      cv::putText(ref_cur_combine_temp, cur_propa, cv::Point2f(pc[0] - 40, pc[1] + 60), cv::FONT_HERSHEY_COMPLEX, 0.4, cv::Scalar(0, 0, 255), 1, 8,
                  0); // 绘制当前真实值文本
      cv::putText(ref_cur_combine_temp, cur_optimize, cv::Point2f(pc[0] - 40, pc[1] + 80), cv::FONT_HERSHEY_COMPLEX, 0.4, cv::Scalar(0, 255, 0), 1, 8,
                  0); // 绘制当前优化估计文本

      cv::rectangle(ref_cur_combine_temp, cv::Point2f(ref_ftr->px_[0] + img_cur.cols - radius, ref_ftr->px_[1] - radius), // 绘制参考补丁矩形
                    cv::Point2f(ref_ftr->px_[0] + img_cur.cols + radius, ref_ftr->px_[1] + radius), cv::Scalar(0, 0, 255), 1);
      cv::rectangle(ref_cur_combine_temp, cv::Point2f(pc[0] - radius, pc[1] - radius), cv::Point2f(pc[0] + radius, pc[1] + radius), // 绘制当前补丁矩形
                    cv::Scalar(0, 255, 0), 1);
      cv::rectangle(ref_cur_combine_temp, cv::Point2f(pc_prior[0] - radius, pc_prior[1] - radius), // 绘制先前补丁矩形
                    cv::Point2f(pc_prior[0] + radius, pc_prior[1] + radius), cv::Scalar(255, 255, 255), 1);
      cv::circle(ref_cur_combine_temp, cv::Point2f(ref_ftr->px_[0] + img_cur.cols, ref_ftr->px_[1]), 1, cv::Scalar(0, 0, 255), -1, 8); // 绘制参考补丁圆
      cv::circle(ref_cur_combine_temp, cv::Point2f(pc[0], pc[1]), 1, cv::Scalar(0, 255, 0), -1, 8); // 绘制当前补丁圆
      cv::circle(ref_cur_combine_temp, cv::Point2f(pc_prior[0], pc_prior[1]), 1, cv::Scalar(255, 255, 255), -1, 8); // 绘制先前补丁圆
      cv::imwrite(dir + std::to_string(new_frame_->id_) + "_" + std::to_string(ref_ftr->id_) + "_" + std::to_string(num) + ".png", // 保存合成图像
                  ref_cur_combine_temp);

      std::vector<std::vector<pixel_member>> pixel_warp_matrix; // 定义像素变换矩阵

      for (int y = 0; y < patch_size; ++y) // 遍历补丁的高度
      {
        vector<pixel_member> pixel_warp_vec; // 定义像素变换向量
        for (int x = 0; x < patch_size; ++x) // 遍历补丁的宽度
        {
          Vector2f px_patch(x - patch_size / 2, y - patch_size / 2); // 计算补丁像素坐标
          px_patch *= (1 << search_level); // 根据搜索层级调整坐标
          const Vector2f px_ref(px_patch + ref_ftr->px_.cast<float>()); // 计算参考补丁的像素坐标
          uint8_t pixel_value = (uint8_t)vk::interpolateMat_8u(img_ref, px_ref[0], px_ref[1]); // 获取参考补丁的像素值

          const Vector2f px(A_cur_ref.cast<float>() * px_patch + pc.cast<float>()); // 计算当前补丁的像素坐标
          if (px[0] < 0 || px[1] < 0 || px[0] >= img_cur.cols - 1 || px[1] >= img_cur.rows - 1) // 如果像素坐标超出范围，则跳过
            continue;
          else
          {
            pixel_member pixel_warp; // 创建像素成员
            pixel_warp.pixel_pos << px[0], px[1]; // 设置像素位置
            pixel_warp.pixel_value = pixel_value; // 设置像素值
            pixel_warp_vec.push_back(pixel_warp); // 添加到向量中
          }
        }
        pixel_warp_matrix.push_back(pixel_warp_vec); // 添加到变换矩阵中
      }

      float x_min = 1000; // 初始化最小x坐标
      float y_min = 1000; // 初始化最小y坐标
      float x_max = 0; // 初始化最大x坐标
      float y_max = 0; // 初始化最大y坐标

      for (int i = 0; i < pixel_warp_matrix.size(); i++) // 遍历像素变换矩阵
      {
        vector<pixel_member> pixel_warp_row = pixel_warp_matrix[i]; // 获取当前行
        for (int j = 0; j < pixel_warp_row.size(); j++) // 遍历当前行的像素
        {
          float x_temp = pixel_warp_row[j].pixel_pos[0]; // 获取x坐标
          float y_temp = pixel_warp_row[j].pixel_pos[1]; // 获取y坐标
          if (x_temp < x_min) x_min = x_temp; // 更新最小x坐标
          if (y_temp < y_min) y_min = y_temp; // 更新最小y坐标
          if (x_temp > x_max) x_max = x_temp; // 更新最大x坐标
          if (y_temp > y_max) y_max = y_temp; // 更新最大y坐标
        }
      }
      int x_min_i = floor(x_min); // 计算最小x坐标的整数值
      int y_min_i = floor(y_min); // 计算最小y坐标的整数值
      int x_max_i = ceil(x_max); // 计算最大x坐标的整数值
      int y_max_i = ceil(y_max); // 计算最大y坐标的整数值
      Matrix2f A_cur_ref_Inv = A_cur_ref.inverse().cast<float>(); // 获取逆变换矩阵
      for (int i = x_min_i; i < x_max_i; i++) // 遍历最小到最大x坐标
      {
        for (int j = y_min_i; j < y_max_i; j++) // 遍历最小到最大y坐标
        {
          Eigen::Vector2f pc_temp(i, j); // 创建临时像素坐标
          Vector2f px_patch = A_cur_ref_Inv * (pc_temp - pc.cast<float>()); // 计算补丁像素坐标
          if (px_patch[0] > (-patch_size / 2 * (1 << search_level)) && px_patch[0] < (patch_size / 2 * (1 << search_level)) && // 如果像素坐标在范围内
              px_patch[1] > (-patch_size / 2 * (1 << search_level)) && px_patch[1] < (patch_size / 2 * (1 << search_level)))
          {
            const Vector2f px_ref(px_patch + ref_ftr->px_.cast<float>()); // 计算参考补丁像素坐标
            uint8_t pixel_value = (uint8_t)vk::interpolateMat_8u(img_ref, px_ref[0], px_ref[1]); // 获取参考补丁像素值
            it_normal[width * j + i] = pixel_value; // 更新法线图像数据
          }
        }
      }
    }
  }
  for (int i = 0; i < visual_submap->voxel_points.size(); i++) // 遍历视觉子图中的所有点
  {
    VisualPoint *pt = visual_submap->voxel_points[i]; // 获取视觉点

    if (!pt->is_normal_initialized_) continue; // 如果法线未初始化，则跳过

    Feature *ref_ftr; // 定义参考特征
    V2D pc(new_frame_->w2c(pt->pos_)); // 将点的世界坐标转换为相机坐标
    ref_ftr = pt->ref_patch; // 获取参考补丁

    Matrix2d A_cur_ref; // 定义变换矩阵
    getWarpMatrixAffine(*cam, ref_ftr->px_, ref_ftr->f_, (ref_ftr->pos() - pt->pos_).norm(), new_frame_->T_f_w_ * ref_ftr->T_f_w_.inverse(), 0, 0,
                        patch_size_half, A_cur_ref); // 获取仿射变换矩阵
    int search_level = getBestSearchLevel(A_cur_ref.inverse(), 2); // 获取最佳搜索层级
    double D = A_cur_ref.determinant(); // 计算变换矩阵的行列式
    if (D > 3) continue; // 如果行列式大于3，则跳过

    cv::Mat img_cur = new_frame_->img_; // 获取当前帧图像
    cv::Mat img_ref = ref_ftr->img_; // 获取参考帧图像
    for (int y = 0; y < patch_size; ++y) // 遍历补丁的高度
    {
      for (int x = 0; x < patch_size; ++x) // 遍历补丁的宽度
      {
        Vector2f px_patch(x - patch_size / 2, y - patch_size / 2); // 计算补丁像素坐标
        px_patch *= (1 << search_level); // 根据搜索层级调整坐标
        const Vector2f px_ref(px_patch + ref_ftr->px_.cast<float>()); // 计算参考补丁像素坐标
        uint8_t pixel_value = (uint8_t)vk::interpolateMat_8u(img_ref, px_ref[0], px_ref[1]); // 获取参考补丁像素值

        const Vector2f px(A_cur_ref.cast<float>() * px_patch + pc.cast<float>()); // 计算当前补丁像素坐标
        if (px[0] < 0 || px[1] < 0 || px[0] >= img_cur.cols - 1 || px[1] >= img_cur.rows - 1) // 如果像素坐标超出范围，则跳过
          continue;
        else
        {
          int col = int(px[0]); // 获取当前像素的列
          int row = int(px[1]); // 获取当前像素的行
          it[width * row + col] = pixel_value; // 更新结果图像数据
        }
      }
    }
  }
  cv::Mat ref_cur_combine; // 定义合成图像
  cv::Mat ref_cur_combine_normal; // 定义法线合成图像
  cv::Mat ref_cur_combine_error; // 定义误差合成图像

  cv::hconcat(result, new_frame_->img_, ref_cur_combine); // 将结果图像与当前帧图像水平拼接
  cv::hconcat(result_normal, new_frame_->img_, ref_cur_combine_normal); // 将法线图像与当前帧图像水平拼接

  cv::cvtColor(ref_cur_combine, ref_cur_combine, CV_GRAY2BGR); // 转换为彩色合成图像
  cv::cvtColor(ref_cur_combine_normal, ref_cur_combine_normal, CV_GRAY2BGR); // 转换为彩色法线合成图像
  cv::absdiff(img_photometric_error, result_normal, img_photometric_error); // 计算光度误差
  cv::hconcat(img_photometric_error, new_frame_->img_, ref_cur_combine_error); // 将光度误差与当前帧图像水平拼接

  cv::imwrite(dir + std::to_string(new_frame_->id_) + "_0_" + ".png", ref_cur_combine); // 保存合成图像
  cv::imwrite(dir + std::to_string(new_frame_->id_) + +"_0_" +
                  "photometric"
                  ".png",
              ref_cur_combine_error); // 保存光度误差图像
  cv::imwrite(dir + std::to_string(new_frame_->id_) + "_0_" + "normal" + ".png", ref_cur_combine_normal); // 保存法线合成图像
}

//预计算参考补丁
void VIOManager::precomputeReferencePatches(int level)
{
  double t1 = omp_get_wtime(); // 记录开始时间
  if (total_points == 0) return; // 如果总点数为0，则返回
  
  MD(1, 2) Jimg; // 定义图像雅可比矩阵
  MD(2, 3) Jdpi; // 定义深度雅可比矩阵
  MD(1, 3) Jdphi, Jdp, JdR, Jdt; // 定义其他雅可比矩阵

  const int H_DIM = total_points * patch_size_total; // 计算雅可比矩阵的维度

  H_sub_inv.resize(H_DIM, 6); // 初始化H_sub_inv矩阵
  H_sub_inv.setZero(); // 将H_sub_inv矩阵置零
  M3D p_w_hat; // 定义3D点的雅可比矩阵

  for (int i = 0; i < total_points; i++) // 遍历所有视觉点
  {
    const int scale = (1 << level); // 计算当前层级的缩放因子

    VisualPoint *pt = visual_submap->voxel_points[i]; // 获取当前视觉点
    cv::Mat img = pt->ref_patch->img_; // 获取当前视觉点的参考补丁图像

    if (pt == nullptr) continue; // 如果视觉点为空，则跳过

    double depth((pt->pos_ - pt->ref_patch->pos()).norm()); // 计算深度
    V3D pf = pt->ref_patch->f_ * depth; // 计算点的世界坐标
    V2D pc = pt->ref_patch->px_; // 获取补丁的像素坐标
    M3D R_ref_w = pt->ref_patch->T_f_w_.rotation_matrix(); // 获取参考补丁的旋转矩阵

    computeProjectionJacobian(pf, Jdpi); // 计算投影雅可比矩阵
    p_w_hat << SKEW_SYM_MATRX(pt->pos_); // 计算点的雅可比矩阵

    const float u_ref = pc[0]; // 获取参考补丁的u坐标
    const float v_ref = pc[1]; // 获取参考补丁的v坐标
    const int u_ref_i = floorf(pc[0] / scale) * scale; // 计算u坐标的缩放索引
    const int v_ref_i = floorf(pc[1] / scale) * scale; // 计算v坐标的缩放索引
    const float subpix_u_ref = (u_ref - u_ref_i) / scale; // 计算u坐标的子像素
    const float subpix_v_ref = (v_ref - v_ref_i) / scale; // 计算v坐标的子像素
    const float w_ref_tl = (1.0 - subpix_u_ref) * (1.0 - subpix_v_ref); // 计算左上角权重
    const float w_ref_tr = subpix_u_ref * (1.0 - subpix_v_ref); // 计算右上角权重
    const float w_ref_bl = (1.0 - subpix_u_ref) * subpix_v_ref; // 计算左下角权重
    const float w_ref_br = subpix_u_ref * subpix_v_ref; // 计算右下角权重

    for (int x = 0; x < patch_size; x++) // 遍历补丁的宽度
    {
      uint8_t *img_ptr = (uint8_t *)img.data + (v_ref_i + x * scale - patch_size_half * scale) * width + u_ref_i - patch_size_half * scale; // 获取补丁图像指针
      for (int y = 0; y < patch_size; ++y, img_ptr += scale) // 遍历补丁的高度
      {
        float du = // 计算u方向的梯度
            0.5f *
            ((w_ref_tl * img_ptr[scale] + w_ref_tr * img_ptr[scale * 2] + w_ref_bl * img_ptr[scale * width + scale] +
              w_ref_br * img_ptr[scale * width + scale * 2]) -
             (w_ref_tl * img_ptr[-scale] + w_ref_tr * img_ptr[0] + w_ref_bl * img_ptr[scale * width - scale] + w_ref_br * img_ptr[scale * width]));
        float dv = // 计算v方向的梯度
            0.5f *
            ((w_ref_tl * img_ptr[scale * width] + w_ref_tr * img_ptr[scale + scale * width] + w_ref_bl * img_ptr[width * scale * 2] +
              w_ref_br * img_ptr[width * scale * 2 + scale]) -
             (w_ref_tl * img_ptr[-scale * width] + w_ref_tr * img_ptr[-scale * width + scale] + w_ref_bl * img_ptr[0] + w_ref_br * img_ptr[scale]));

        Jimg << du, dv; // 更新图像雅可比矩阵
        Jimg = Jimg * (1.0 / scale); // 归一化雅可比矩阵

        JdR = Jimg * Jdpi * R_ref_w * p_w_hat; // 计算旋转雅可比矩阵
        Jdt = -Jimg * Jdpi * R_ref_w; // 计算平移雅可比矩阵

        H_sub_inv.block<1, 6>(i * patch_size_total + x * patch_size + y, 0) << JdR, Jdt; // 更新H_sub_inv矩阵
      }
    }
  }
  has_ref_patch_cache = true; // 标记参考补丁缓存已计算
}

//更新逆状态
void VIOManager::updateStateInverse(cv::Mat img, int level)
{
  if (total_points == 0) return; // 如果总点数为0，则返回
  StatesGroup old_state = (*state); // 保存旧状态
  V2D pc; // 定义相机坐标
  MD(1, 2) Jimg; // 定义图像雅可比矩阵
  MD(2, 3) Jdpi; // 定义深度雅可比矩阵
  MD(1, 3) Jdphi, Jdp, JdR, Jdt; // 定义其他雅可比矩阵
  VectorXd z; // 定义测量向量
  MatrixXd H_sub; // 定义雅可比矩阵
  bool EKF_end = false; // 定义EKF结束标志
  float last_error = std::numeric_limits<float>::max(); // 初始化上一次误差为最大值
  compute_jacobian_time = update_ekf_time = 0.0; // 初始化计算时间
  M3D P_wi_hat; // 定义3D点的雅可比矩阵
  bool z_init = true; // 定义z初始化标志
  const int H_DIM = total_points * patch_size_total; // 计算雅可比矩阵的维度

  z.resize(H_DIM); // 初始化测量向量
  z.setZero(); // 将测量向量置零

  H_sub.resize(H_DIM, 6); // 初始化雅可比矩阵
  H_sub.setZero(); // 将雅可比矩阵置零

  for (int iteration = 0; iteration < max_iterations; iteration++) // 迭代更新状态
  {
    double t1 = omp_get_wtime(); // 记录开始时间
    double count_outlier = 0; // 初始化异常值计数
    if (has_ref_patch_cache == false) precomputeReferencePatches(level); // 如果没有缓存，预计算参考补丁
    int n_meas = 0; // 初始化测量数量
    float error = 0.0; // 初始化误差
    M3D Rwi(state->rot_end); // 获取当前状态的旋转矩阵
    V3D Pwi(state->pos_end); // 获取当前状态的位移
    P_wi_hat << SKEW_SYM_MATRX(Pwi); // 计算点的雅可比矩阵
    Rcw = Rci * Rwi.transpose(); // 计算当前帧的旋转矩阵
    Pcw = -Rci * Rwi.transpose() * Pwi + Pci; // 计算当前帧的位移

    M3D p_hat; // 定义3D点的雅可比矩阵

    for (int i = 0; i < total_points; i++) // 遍历所有视觉点
    {
      float patch_error = 0.0; // 初始化补丁误差

      const int scale = (1 << level); // 计算缩放因子

      VisualPoint *pt = visual_submap->voxel_points[i]; // 获取当前视觉点

      if (pt == nullptr) continue; // 如果视觉点为空，则跳过

      V3D pf = Rcw * pt->pos_ + Pcw; // 计算点的世界坐标
      pc = cam->world2cam(pf); // 将世界坐标转换为相机坐标

      const float u_ref = pc[0]; // 获取u坐标
      const float v_ref = pc[1]; // 获取v坐标
      const int u_ref_i = floorf(pc[0] / scale) * scale; // 计算u坐标的缩放索引
      const int v_ref_i = floorf(pc[1] / scale) * scale; // 计算v坐标的缩放索引
      const float subpix_u_ref = (u_ref - u_ref_i) / scale; // 计算u坐标的子像素
      const float subpix_v_ref = (v_ref - v_ref_i) / scale; // 计算v坐标的子像素
      const float w_ref_tl = (1.0 - subpix_u_ref) * (1.0 - subpix_v_ref); // 计算左上角权重
      const float w_ref_tr = subpix_u_ref * (1.0 - subpix_v_ref); // 计算右上角权重
      const float w_ref_bl = (1.0 - subpix_u_ref) * subpix_v_ref; // 计算左下角权重
      const float w_ref_br = subpix_u_ref * subpix_v_ref; // 计算右下角权重

      vector<float> P = visual_submap->warp_patch[i]; // 获取当前视觉点的补丁
      for (int x = 0; x < patch_size; x++) // 遍历补丁的宽度
      {
        uint8_t *img_ptr = (uint8_t *)img.data + (v_ref_i + x * scale - patch_size_half * scale) * width + u_ref_i - patch_size_half * scale; // 获取补丁图像指针
        for (int y = 0; y < patch_size; ++y, img_ptr += scale) // 遍历补丁的高度
        {
          double res = w_ref_tl * img_ptr[0] + w_ref_tr * img_ptr[scale] + w_ref_bl * img_ptr[scale * width] +
                       w_ref_br * img_ptr[scale * width + scale] - P[patch_size_total * level + x * patch_size + y]; // 计算残差
          z(i * patch_size_total + x * patch_size + y) = res; // 更新测量向量
          patch_error += res * res; // 更新补丁误差
          MD(1, 3) J_dR = H_sub_inv.block<1, 3>(i * patch_size_total + x * patch_size + y, 0); // 获取旋转雅可比矩阵
          MD(1, 3) J_dt = H_sub_inv.block<1, 3>(i * patch_size_total + x * patch_size + y, 3); // 获取平移雅可比矩阵
          JdR = J_dR * Rwi + J_dt * P_wi_hat * Rwi; // 计算旋转雅可比矩阵
          Jdt = J_dt * Rwi; // 计算平移雅可比矩阵
          H_sub.block<1, 6>(i * patch_size_total + x * patch_size + y, 0) << JdR, Jdt; // 更新雅可比矩阵
          n_meas++; // 增加测量数量
        }
      }
      visual_submap->errors[i] = patch_error; // 更新视觉子图中的误差
      error += patch_error; // 更新总误差
    }

    error = error / n_meas; // 计算平均误差

    compute_jacobian_time += omp_get_wtime() - t1; // 计算雅可比矩阵的时间

    double t3 = omp_get_wtime(); // 记录计算时间

    if (error <= last_error) // 如果当前误差小于等于上一次误差
    {
      old_state = (*state); // 保存当前状态
      last_error = error; // 更新上一次误差

      auto &&H_sub_T = H_sub.transpose(); // 获取雅可比矩阵的转置
      H_T_H.setZero(); // 将H_T_H矩阵置零
      G.setZero(); // 将G矩阵置零
      H_T_H.block<6, 6>(0, 0) = H_sub_T * H_sub; // 计算H_T_H矩阵
      MD(DIM_STATE, DIM_STATE) &&K_1 = (H_T_H + (state->cov / img_point_cov).inverse()).inverse(); // 计算卡尔曼增益
      auto &&HTz = H_sub_T * z; // 计算H_Tz
      auto vec = (*state_propagat) - (*state); // 计算状态传播的差值
      G.block<DIM_STATE, 6>(0, 0) = K_1.block<DIM_STATE, 6>(0, 0) * H_T_H.block<6, 6>(0, 0); // 更新G矩阵
      auto solution = -K_1.block<DIM_STATE, 6>(0, 0) * HTz + vec - G.block<DIM_STATE, 6>(0, 0) * vec.block<6, 1>(0, 0); // 计算状态更新
      (*state) += solution; // 更新状态
      auto &&rot_add = solution.block<3, 1>(0, 0); // 获取旋转增量
      auto &&t_add = solution.block<3, 1>(3, 0); // 获取位移增量

      if ((rot_add.norm() * 57.3f < 0.001f) && (t_add.norm() * 100.0f < 0.001f)) { EKF_end = true; } // 如果旋转和位移增量都小于阈值，则结束EKF
    }
    else
    {
      (*state) = old_state; // 如果当前误差大于上一次误差，则恢复旧状态
      EKF_end = true; // 设置EKF结束标志
    }

    update_ekf_time += omp_get_wtime() - t3; // 计算EKF更新的时间

    if (iteration == max_iterations || EKF_end) break; // 如果达到最大迭代次数或结束标志，则跳出循环
  }
}
//更新状态
void VIOManager::updateState(cv::Mat img, int level)
{
  if (total_points == 0) return; // 如果总点数为0，则返回
  StatesGroup old_state = (*state); // 保存旧状态

  VectorXd z; // 定义测量向量
  MatrixXd H_sub; // 定义雅可比矩阵
  bool EKF_end = false; // 定义EKF结束标志
  float last_error = std::numeric_limits<float>::max(); // 初始化上一次误差为最大值

  const int H_DIM = total_points * patch_size_total; // 计算雅可比矩阵的维度
  z.resize(H_DIM); // 初始化测量向量
  z.setZero(); // 将测量向量置零
  H_sub.resize(H_DIM, 7); // 初始化雅可比矩阵
  H_sub.setZero(); // 将雅可比矩阵置零

  for (int iteration = 0; iteration < max_iterations; iteration++) // 迭代更新状态
  {
    double t1 = omp_get_wtime(); // 记录开始时间

    M3D Rwi(state->rot_end); // 获取当前状态的旋转矩阵
    V3D Pwi(state->pos_end); // 获取当前状态的位移
    Rcw = Rci * Rwi.transpose(); // 计算当前帧的旋转矩阵
    Pcw = -Rci * Rwi.transpose() * Pwi + Pci; // 计算当前帧的位移
    Jdp_dt = Rci * Rwi.transpose(); // 计算平移雅可比矩阵
    
    float error = 0.0; // 初始化误差
    int n_meas = 0; // 初始化测量数量
    // int max_threads = omp_get_max_threads();
    // int desired_threads = std::min(max_threads, total_points);
    // omp_set_num_threads(desired_threads);
  
    #ifdef MP_EN
      omp_set_num_threads(MP_PROC_NUM); // 设置并行处理的线程数
      #pragma omp parallel for reduction(+:error, n_meas) // 并行处理
    #endif
    for (int i = 0; i < total_points; i++) // 遍历所有视觉点
    {
      // printf("thread is %d, i=%d, i address is %p\n", omp_get_thread_num(), i, &i);
      MD(1, 2) Jimg; // 定义图像雅可比矩阵
      MD(2, 3) Jdpi; // 定义深度雅可比矩阵
      MD(1, 3) Jdphi, Jdp, JdR, Jdt; // 定义其他雅可比矩阵

      float patch_error = 0.0; // 初始化补丁误差
      int search_level = visual_submap->search_levels[i]; // 获取搜索层级
      int pyramid_level = level + search_level; // 计算金字塔层级
      int scale = (1 << pyramid_level); // 计算缩放因子
      float inv_scale = 1.0f / scale; // 计算缩放因子的倒数

      VisualPoint *pt = visual_submap->voxel_points[i]; // 获取当前视觉点

      if (pt == nullptr) continue; // 如果视觉点为空，则跳过

      V3D pf = Rcw * pt->pos_ + Pcw; // 计算点的世界坐标
      V2D pc = cam->world2cam(pf); // 将世界坐标转换为相机坐标

      computeProjectionJacobian(pf, Jdpi); // 计算投影雅可比矩阵
      M3D p_hat; // 定义3D点的雅可比矩阵
      p_hat << SKEW_SYM_MATRX(pf); // 计算点的雅可比矩阵

      float u_ref = pc[0]; // 获取u坐标
      float v_ref = pc[1]; // 获取v坐标
      int u_ref_i = floorf(pc[0] / scale) * scale; // 计算u坐标的缩放索引
      int v_ref_i = floorf(pc[1] / scale) * scale; // 计算v坐标的缩放索引
      float subpix_u_ref = (u_ref - u_ref_i) / scale; // 计算u坐标的子像素
      float subpix_v_ref = (v_ref - v_ref_i) / scale; // 计算v坐标的子像素
      float w_ref_tl = (1.0 - subpix_u_ref) * (1.0 - subpix_v_ref); // 计算左上角权重
      float w_ref_tr = subpix_u_ref * (1.0 - subpix_v_ref); // 计算右上角权重
      float w_ref_bl = (1.0 - subpix_u_ref) * subpix_v_ref; // 计算左下角权重
      float w_ref_br = subpix_u_ref * subpix_v_ref; // 计算右下角权重

      vector<float> P = visual_submap->warp_patch[i]; // 获取当前视觉点的补丁
      double inv_ref_expo = visual_submap->inv_expo_list[i]; // 获取参考曝光时间
      // ROS_ERROR("inv_ref_expo: %.3lf, state->inv_expo_time: %.3lf\n", inv_ref_expo, state->inv_expo_time);

      for (int x = 0; x < patch_size; x++) // 遍历补丁的宽度
      {
        uint8_t *img_ptr = (uint8_t *)img.data + (v_ref_i + x * scale - patch_size_half * scale) * width + u_ref_i - patch_size_half * scale; // 获取补丁图像指针
        for (int y = 0; y < patch_size; ++y, img_ptr += scale) // 遍历补丁的高度
        {
          float du = // 计算u方向的梯度
              0.5f *
              ((w_ref_tl * img_ptr[scale] + w_ref_tr * img_ptr[scale * 2] + w_ref_bl * img_ptr[scale * width + scale] +
                w_ref_br * img_ptr[scale * width + scale * 2]) -
               (w_ref_tl * img_ptr[-scale] + w_ref_tr * img_ptr[0] + w_ref_bl * img_ptr[scale * width - scale] + w_ref_br * img_ptr[scale * width]));
          float dv = // 计算v方向的梯度
              0.5f *
              ((w_ref_tl * img_ptr[scale * width] + w_ref_tr * img_ptr[scale + scale * width] + w_ref_bl * img_ptr[width * scale * 2] +
                w_ref_br * img_ptr[width * scale * 2 + scale]) -
               (w_ref_tl * img_ptr[-scale * width] + w_ref_tr * img_ptr[-scale * width + scale] + w_ref_bl * img_ptr[0] + w_ref_br * img_ptr[scale]));

          Jimg << du, dv; // 更新图像雅可比矩阵
          Jimg = Jimg * state->inv_expo_time; // 乘以反曝光时间
          Jimg = Jimg * inv_scale; // 归一化雅可比矩阵
          Jdphi = Jimg * Jdpi * p_hat; // 计算旋转雅可比矩阵
          Jdp = -Jimg * Jdpi; // 计算平移雅可比矩阵
          JdR = Jdphi * Jdphi_dR + Jdp * Jdp_dR; // 计算旋转雅可比矩阵
          Jdt = Jdp * Jdp_dt; // 计算平移雅可比矩阵

          double cur_value = // 计算当前补丁的值
              w_ref_tl * img_ptr[0] + w_ref_tr * img_ptr[scale] + w_ref_bl * img_ptr[scale * width] + w_ref_br * img_ptr[scale * width + scale];
          double res = state->inv_expo_time * cur_value - inv_ref_expo * P[patch_size_total * level + x * patch_size + y]; // 计算残差

          z(i * patch_size_total + x * patch_size + y) = res; // 更新测量向量

          patch_error += res * res; // 更新补丁误差
          n_meas += 1; // 增加测量数量
          
          if (exposure_estimate_en) { H_sub.block<1, 7>(i * patch_size_total + x * patch_size + y, 0) << JdR, Jdt, cur_value; } // 如果启用曝光估计，更新H_sub矩阵
          else { H_sub.block<1, 6>(i * patch_size_total + x * patch_size + y, 0) << JdR, Jdt; } // 否则只更新旋转和平移雅可比矩阵
        }
      }
      visual_submap->errors[i] = patch_error; // 更新视觉子图中的误差
      error += patch_error; // 更新总误差
    }

    error = error / n_meas; // 计算平均误差
    
    compute_jacobian_time += omp_get_wtime() - t1; // 计算雅可比矩阵的时间

    // printf("\nPYRAMID LEVEL %i\n---------------\n", level);
    // std::cout << "It. " << iteration
    //           << "\t last_error = " << last_error
    //           << "\t new_error = " << error
    //           << std::endl;

    double t3 = omp_get_wtime(); // 记录计算时间

    if (error <= last_error) // 如果当前误差小于等于上一次误差
    {
      old_state = (*state); // 保存当前状态
      last_error = error; // 更新上一次误差

      auto &&H_sub_T = H_sub.transpose(); // 获取雅可比矩阵的转置
      H_T_H.setZero(); // 将H_T_H矩阵置零
      G.setZero(); // 将G矩阵置零
      H_T_H.block<7, 7>(0, 0) = H_sub_T * H_sub; // 计算H_T_H矩阵
      MD(DIM_STATE, DIM_STATE) &&K_1 = (H_T_H + (state->cov / img_point_cov).inverse()).inverse(); // 计算卡尔曼增益
      auto &&HTz = H_sub_T * z; // 计算H_Tz
      // K = K_1.block<DIM_STATE,6>(0,0) * H_sub_T;
      auto vec = (*state_propagat) - (*state); // 计算状态传播的差值
      G.block<DIM_STATE, 7>(0, 0) = K_1.block<DIM_STATE, 7>(0, 0) * H_T_H.block<7, 7>(0, 0); // 更新G矩阵
      MD(DIM_STATE, 1)
      solution = -K_1.block<DIM_STATE, 7>(0, 0) * HTz + vec - G.block<DIM_STATE, 7>(0, 0) * vec.block<7, 1>(0, 0); // 计算状态更新

      (*state) += solution; // 更新状态
      auto &&rot_add = solution.block<3, 1>(0, 0); // 获取旋转增量
      auto &&t_add = solution.block<3, 1>(3, 0); // 获取位移增量

      auto &&expo_add = solution.block<1, 1>(6, 0); // 获取曝光增量
      // if ((rot_add.norm() * 57.3f < 0.001f) && (t_add.norm() * 100.0f < 0.001f) && (expo_add.norm() < 0.001f)) EKF_end = true;
      if ((rot_add.norm() * 57.3f < 0.001f) && (t_add.norm() * 100.0f < 0.001f))  EKF_end = true; // 如果旋转和位移增量都小于阈值，则结束EKF
    }
    else
    {
      (*state) = old_state; // 如果当前误差大于上一次误差，则恢复旧状态
      EKF_end = true; // 设置EKF结束标志
    }

    update_ekf_time += omp_get_wtime() - t3; // 计算EKF更新的时间

    if (iteration == max_iterations || EKF_end) break; // 如果达到最大迭代次数或结束标志，则跳出循环
  }
  // if (state->inv_expo_time < 0.0)  {ROS_ERROR("reset expo time!!!!!!!!!!\n"); state->inv_expo_time = 0.0;}
}

void VIOManager::updateFrameState(StatesGroup state)
{
  M3D Rwi(state.rot_end); // 获取当前状态的旋转矩阵
  V3D Pwi(state.pos_end); // 获取当前状态的位移
  Rcw = Rci * Rwi.transpose(); // 计算当前帧的旋转矩阵
  Pcw = -Rci * Rwi.transpose() * Pwi + Pci; // 计算当前帧的位移
  new_frame_->T_f_w_ = SE3(Rcw, Pcw); // 更新当前帧的变换矩阵
}

void VIOManager::plotTrackedPoints()
{
  int total_points = visual_submap->voxel_points.size(); // 获取视觉子图中的点数量
  if (total_points == 0) return; // 如果点数量为0，则返回
  // int inlier_count = 0;
  // for (int i = 0; i < img_cp.rows / grid_size; i++)
  // {
  //   cv::line(img_cp, cv::Poaint2f(0, grid_size * i), cv::Point2f(img_cp.cols, grid_size * i), cv::Scalar(255, 255, 255), 1, CV_AA);
  // }
  // for (int i = 0; i < img_cp.cols / grid_size; i++)
  // {
  //   cv::line(img_cp, cv::Point2f(grid_size * i, 0), cv::Point2f(grid_size * i, img_cp.rows), cv::Scalar(255, 255, 255), 1, CV_AA);
  // }
  // for (int i = 0; i < img_cp.rows / grid_size; i++)
  // {
  //   cv::line(img_cp, cv::Point2f(0, grid_size * i), cv::Point2f(img_cp.cols, grid_size * i), cv::Scalar(255, 255, 255), 1, CV_AA);
  // }
  // for (int i = 0; i < img_cp.cols / grid_size; i++)
  // {
  //   cv::line(img_cp, cv::Point2f(grid_size * i, 0), cv::Point2f(grid_size * i, img_cp.rows), cv::Scalar(255, 255, 255), 1, CV_AA);
  // }
  for (int i = 0; i < total_points; i++) // 遍历所有视觉点
  {
    VisualPoint *pt = visual_submap->voxel_points[i]; // 获取当前视觉点
    V2D pc(new_frame_->w2c(pt->pos_)); // 将世界坐标转换为相机坐标

    if (visual_submap->errors[i] <= visual_submap->propa_errors[i]) // 如果误差小于等于传播误差
    {
      // inlier_count++;
      cv::circle(img_cp, cv::Point2f(pc[0], pc[1]), 7, cv::Scalar(0, 255, 0), -1, 8); // 绘制绿色圆圈表示跟踪成功
    }
    else
    {
      cv::circle(img_cp, cv::Point2f(pc[0], pc[1]), 7, cv::Scalar(255, 0, 0), -1, 8); // 绘制红色圆圈表示跟踪失败
    }
  }
  // std::string text = std::to_string(inlier_count) + " " + std::to_string(total_points);
  // cv::Point2f origin;
  // origin.x = img_cp.cols - 110;
  // origin.y = 20;
  // cv::putText(img_cp, text, origin, cv::FONT_HERSHEY_COMPLEX, 0.7, cv::Scalar(0, 255, 0), 2, 8, 0);
}
//获取插值像素
V3F VIOManager::getInterpolatedPixel(cv::Mat img, V2D pc)
{
  const float u_ref = pc[0]; // 获取u坐标
  const float v_ref = pc[1]; // 获取v坐标
  const int u_ref_i = floorf(pc[0]); // 计算u坐标的整数部分
  const int v_ref_i = floorf(pc[1]); // 计算v坐标的整数部分
  const float subpix_u_ref = (u_ref - u_ref_i); // 计算u坐标的子像素
  const float subpix_v_ref = (v_ref - v_ref_i); // 计算v坐标的子像素
  const float w_ref_tl = (1.0 - subpix_u_ref) * (1.0 - subpix_v_ref); // 计算左上角权重
  const float w_ref_tr = subpix_u_ref * (1.0 - subpix_v_ref); // 计算右上角权重
  const float w_ref_bl = (1.0 - subpix_u_ref) * subpix_v_ref; // 计算左下角权重
  const float w_ref_br = subpix_u_ref * subpix_v_ref; // 计算右下角权重
  uint8_t *img_ptr = (uint8_t *)img.data + ((v_ref_i)*width + (u_ref_i)) * 3; // 获取图像指针
  float B = w_ref_tl * img_ptr[0] + w_ref_tr * img_ptr[0 + 3] + w_ref_bl * img_ptr[width * 3] + w_ref_br * img_ptr[width * 3 + 0 + 3]; // 计算B通道值
  float G = w_ref_tl * img_ptr[1] + w_ref_tr * img_ptr[1 + 3] + w_ref_bl * img_ptr[1 + width * 3] + w_ref_br * img_ptr[width * 3 + 1 + 3]; // 计算G通道值
  float R = w_ref_tl * img_ptr[2] + w_ref_tr * img_ptr[2 + 3] + w_ref_bl * img_ptr[2 + width * 3] + w_ref_br * img_ptr[width * 3 + 2 + 3]; // 计算R通道值
  V3F pixel(B, G, R); // 返回像素值
  return pixel;
}
//截止到当前帧的特征点
void VIOManager::dumpDataForColmap()
{
  static int cnt = 1; // 初始化计数器
  std::ostringstream ss; // 定义字符串流
  ss << std::setw(5) << std::setfill('0') << cnt; // 格式化计数器
  std::string cnt_str = ss.str(); // 获取格式化后的字符串
  std::string image_path = std::string(ROOT_DIR) + "Log/Colmap/images/" + cnt_str + ".png"; // 定义图像路径
  
  cv::Mat img_rgb_undistort; // 定义去畸变图像
  pinhole_cam->undistortImage(img_rgb, img_rgb_undistort); // 进行去畸变处理
  cv::imwrite(image_path, img_rgb_undistort); // 保存去畸变图像
  
  Eigen::Quaterniond q(new_frame_->T_f_w_.rotation_matrix()); // 获取当前帧的四元数表示
  Eigen::Vector3d t = new_frame_->T_f_w_.translation(); // 获取当前帧的位移
  fout_colmap << cnt << " " // 写入数据到Colmap文件
            << std::fixed << std::setprecision(6)  // 保证浮点数精度为6位
            << q.w() << " " << q.x() << " " << q.y() << " " << q.z() << " "
            << t.x() << " " << t.y() << " " << t.z() << " "
            << 1 << " "  // CAMERA_ID (假设相机ID为1)
            << cnt_str << ".png" << std::endl;
  fout_colmap << "0.0 0.0 -1" << std::endl; // 写入相机的外参
  cnt++; // 增加计数器
}
//更新视觉地图点
void VIOManager::processFrame(cv::Mat &img, vector<pointWithVar> &pg, const unordered_map<VOXEL_LOCATION, VoxelOctoTree *> &feat_map, double img_time)
{
  if (width != img.cols || height != img.rows) // 如果图像尺寸不匹配
  {
    if (img.empty()) printf("[ VIO ] Empty Image!\n"); // 如果图像为空，则输出提示
    cv::resize(img, img, cv::Size(img.cols * image_resize_factor, img.rows * image_resize_factor), 0, 0, CV_INTER_LINEAR); // 缩放图像
  }
  img_rgb = img.clone(); // 复制图像
  img_cp = img.clone(); // 复制图像用于显示
  // img_test = img.clone(); // 可选：用于测试

  if (img.channels() == 3) cv::cvtColor(img, img, CV_BGR2GRAY); // 如果图像为彩色，则转换为灰度图

  new_frame_.reset(new Frame(cam, img)); // 创建新的帧对象
  updateFrameState(*state); // 更新帧状态
  
  resetGrid(); // 重置网格

  double t1 = omp_get_wtime(); // 记录开始时间

  retrieveFromVisualSparseMap(img, pg, feat_map); // 从视觉稀疏地图中检索特征

  double t2 = omp_get_wtime(); // 记录时间

  computeJacobianAndUpdateEKF(img); // 计算雅可比矩阵并更新扩展卡尔曼滤波器（EKF）

  double t3 = omp_get_wtime(); // 记录计算时间

  generateVisualMapPoints(img, pg); // 生成视觉地图点

  double t4 = omp_get_wtime(); // 记录时间
  
  plotTrackedPoints(); // 绘制被跟踪的点

  if (plot_flag) projectPatchFromRefToCur(feat_map); // 如果启用绘图标志，则将参考补丁投影到当前帧

  double t5 = omp_get_wtime(); // 记录时间

  updateVisualMapPoints(img); // 更新视觉地图点

  double t6 = omp_get_wtime(); // 记录时间

  updateReferencePatch(feat_map); // 更新参考补丁

  double t7 = omp_get_wtime(); // 记录时间
  
  if(colmap_output_en) dumpDataForColmap(); // 如果启用Colmap输出，则转储数据

  frame_count++; // 增加帧计数
  ave_total = ave_total * (frame_count - 1) / frame_count + (t7 - t1 - (t5 - t4)) / frame_count; // 更新平均时间

  // printf("[ VIO ] feat_map.size(): %zu\n", feat_map.size());
  // printf("\033[1;32m[ VIO time ]: current frame: retrieveFromVisualSparseMap time: %.6lf secs.\033[0m\n", t2 - t1);
  // printf("\033[1;32m[ VIO time ]: current frame: computeJacobianAndUpdateEKF time: %.6lf secs, comp H: %.6lf secs, ekf: %.6lf secs.\033[0m\n", t3 - t2, computeH, ekf_time);
  // printf("\033[1;32m[ VIO time ]: current frame: generateVisualMapPoints time: %.6lf secs.\033[0m\n", t4 - t3);
  // printf("\033[1;32m[ VIO time ]: current frame: updateVisualMapPoints time: %.6lf secs.\033[0m\n", t6 - t5);
  // printf("\033[1;32m[ VIO time ]: current frame: updateReferencePatch time: %.6lf secs.\033[0m\n", t7 - t6);
  // printf("\033[1;32m[ VIO time ]: current total time: %.6lf, average total time: %.6lf secs.\033[0m\n", t7 - t1 - (t5 - t4), ave_total);

  // ave_build_residual_time = ave_build_residual_time * (frame_count - 1) / frame_count + (t2 - t1) / frame_count;
  // ave_ekf_time = ave_ekf_time * (frame_count - 1) / frame_count + (t3 - t2) / frame_count;
 
  // cout << BLUE << "ave_build_residual_time: " << ave_build_residual_time << RESET << endl;
  // cout << BLUE << "ave_ekf_time: " << ave_ekf_time << RESET << endl;
  
  printf("\033[1;34m+-------------------------------------------------------------+\033[0m\n");
  printf("\033[1;34m|                         VIO Time                            |\033[0m\n");
  printf("\033[1;34m+-------------------------------------------------------------+\033[0m\n");
  printf("\033[1;34m| %-29s | %-27zu |\033[0m\n", "Sparse Map Size", feat_map.size()); // 输出稀疏地图的大小
  printf("\033[1;34m+-------------------------------------------------------------+\033[0m\n");
  printf("\033[1;34m| %-29s | %-27s |\033[0m\n", "Algorithm Stage", "Time (secs)"); // 输出算法阶段和时间
  printf("\033[1;34m+-------------------------------------------------------------+\033[0m\n");
  printf("\033[1;32m| %-29s | %-27lf |\033[0m\n", "retrieveFromVisualSparseMap", t2 - t1); // 输出从视觉稀疏地图检索的时间
  printf("\033[1;32m| %-29s | %-27lf |\033[0m\n", "computeJacobianAndUpdateEKF", t3 - t2); // 输出计算雅可比矩阵和更新EKF的时间
  printf("\033[1;32m| %-27s   | %-27lf |\033[0m\n", "-> computeJacobian", compute_jacobian_time); // 输出计算雅可比矩阵的时间
  printf("\033[1;32m| %-27s   | %-27lf |\033[0m\n", "-> updateEKF", update_ekf_time); // 输出更新EKF的时间
  printf("\033[1;32m| %-29s | %-27lf |\033[0m\n", "generateVisualMapPoints", t4 - t3); // 输出生成视觉地图点的时间
  printf("\033[1;32m| %-29s | %-27lf |\033[0m\n", "updateVisualMapPoints", t6 - t5); // 输出更新视觉地图点的时间
  printf("\033[1;32m| %-29s | %-27lf |\033[0m\n", "updateReferencePatch", t7 - t6); // 输出更新参考补丁的时间
  printf("\033[1;34m+-------------------------------------------------------------+\033[0m\n");
  printf("\033[1;32m| %-29s | %-27lf |\033[0m\n", "Current Total Time", t7 - t1 - (t5 - t4)); // 输出当前总时间
  printf("\033[1;32m| %-29s | %-27lf |\033[0m\n", "Average Total Time", ave_total); // 输出平均总时间
  printf("\033[1;34m+-------------------------------------------------------------+\033[0m\n");

  // std::string text = std::to_string(int(1 / (t7 - t1 - (t5 - t4)))) + " HZ";
  // cv::Point2f origin;
  // origin.x = 20;
  // origin.y = 20;
  // cv::putText(img_cp, text, origin, cv::FONT_HERSHEY_COMPLEX, 0.6, cv::Scalar(255, 255, 255), 1, 8, 0);
  // cv::imwrite("/home/chunran/Desktop/raycasting/" + std::to_string(new_frame_->id_) + ".png", img_cp);
}