/* 
This file is part of FAST-LIVO2: Fast, Direct LiDAR-Inertial-Visual Odometry.

Developer: Chunran Zheng <zhengcr@connect.hku.hk>

For commercial use, please contact me at <zhengcr@connect.hku.hk> or
Prof. Fu Zhang at <fuzhang@hku.hk>.

This file is subject to the terms and conditions outlined in the 'LICENSE' file,
which is included as part of this source code package.
*/

#include <boost/bind.hpp>
#include "feature.h" // 导入特征的定义
#include "frame.h"   // 导入帧的定义
#include "visual_point.h" // 导入视觉点的定义
#include <stdexcept> // 导入异常处理
#include <vikit/math_utils.h> // 导入数学工具
#include <vikit/performance_monitor.h> // 导入性能监控工具
#include <vikit/vision.h> // 导入视觉相关工具

// 静态成员变量初始化
int Frame::frame_counter_ = 0;

// 帧的构造函数
Frame::Frame(vk::AbstractCamera *cam, const cv::Mat &img)
    : id_(frame_counter_++),  // 为帧分配唯一的ID
      cam_(cam) // 关联相机对象
{
  initFrame(img); // 初始化帧
}

// 帧的析构函数
Frame::~Frame()
{
  // 删除所有特征以释放内存
  std::for_each(fts_.begin(), fts_.end(), [&](Feature *i) { delete i; });
}

// 初始化帧的函数
void Frame::initFrame(const cv::Mat &img)
{
  // 检查图像是否为空
  if (img.empty()) { throw std::runtime_error("Frame: provided image is empty"); }

  // 检查图像大小是否与相机模型一致
  if (img.cols != cam_->width() || img.rows != cam_->height())
  {
    throw std::runtime_error("Frame: provided image has not the same size as the camera model");
  }

  // 检查图像类型是否为单通道灰度图像
  if (img.type() != CV_8UC1) { throw std::runtime_error("Frame: provided image is not grayscale"); }

  img_ = img; // 将图像赋值给帧
}

/// 帧类的实用函数
namespace frame_utils
{

  // 创建图像金字塔
  void createImgPyramid(const cv::Mat &img_level_0, int n_levels, ImgPyr &pyr)
  {
    pyr.resize(n_levels); // 调整金字塔大小
    pyr[0] = img_level_0; // 将第0层设为输入图像
    for (int i = 1; i < n_levels; ++i)
    {
      // 创建更小的图像，尺寸减半
      pyr[i] = cv::Mat(pyr[i - 1].rows / 2, pyr[i - 1].cols / 2, CV_8U);
      vk::halfSample(pyr[i - 1], pyr[i]); // 进行下采样
    }
  }

} // namespace frame_utils
