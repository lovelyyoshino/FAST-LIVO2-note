/* 
This file is part of FAST-LIVO2: Fast, Direct LiDAR-Inertial-Visual Odometry.

Developer: Chunran Zheng <zhengcr@connect.hku.hk>

For commercial use, please contact me at <zhengcr@connect.hku.hk> or
Prof. Fu Zhang at <fuzhang@hku.hk>.

This file is subject to the terms and conditions outlined in the 'LICENSE' file,
which is included as part of this source code package.
*/

#include "visual_point.h"
#include "feature.h"
#include <stdexcept>
#include <vikit/math_utils.h>

VisualPoint::VisualPoint(const Vector3d &pos)
    : pos_(pos), previous_normal_(Vector3d::Zero()), normal_(Vector3d::Zero()),
      is_converged_(false), is_normal_initialized_(false), has_ref_patch_(false)
{
    // 构造函数，初始化视觉点的位置和其他成员变量
}

VisualPoint::~VisualPoint() 
{
    // 析构函数，释放观察特征的内存
    for (auto it = obs_.begin(), ite = obs_.end(); it != ite; ++it)
    {
        delete(*it); // 删除每个观察特征
    }
    obs_.clear(); // 清空观察特征列表
    ref_patch = nullptr; // 将参考补丁指针置为nullptr
}
//添加观察特征
void VisualPoint::addFrameRef(Feature *ftr)
{
    // 将特征添加到观察列表的前面
    obs_.push_front(ftr);
}
//删除观察特征
void VisualPoint::deleteFeatureRef(Feature *ftr)
{
    // 删除观察列表中的特征
    if (ref_patch == ftr) // 如果待删除特征是参考补丁
    {
        ref_patch = nullptr; // 置空参考补丁指针
        has_ref_patch_ = false; // 更新参考补丁状态
    }
    for (auto it = obs_.begin(), ite = obs_.end(); it != ite; ++it)
    {
        if ((*it) == ftr) // 找到待删除特征
        {
            delete((*it)); // 删除特征
            obs_.erase(it); // 从观察列表中移除
            return;
        }
    }
}
//获取与当前帧位置相近的观察特征
bool VisualPoint::getCloseViewObs(const Vector3d &framepos, Feature *&ftr, const Vector2d &cur_px) const
{
    // 获取与当前帧位置相近的观察特征
    if (obs_.size() <= 0) return false; // 如果没有观察特征，返回false

    Vector3d obs_dir(framepos - pos_); // 计算观察方向
    obs_dir.normalize(); // 归一化观察方向
    auto min_it = obs_.begin(); // 初始化最小角度的迭代器
    double min_cos_angle = 0; // 初始化最小余弦角度

    for (auto it = obs_.begin(), ite = obs_.end(); it != ite; ++it)
    {
        Vector3d dir((*it)->T_f_w_.inverse().translation() - pos_); // 计算特征方向
        dir.normalize(); // 归一化特征方向
        double cos_angle = obs_dir.dot(dir); // 计算观察方向与特征方向的余弦角度
        if (cos_angle > min_cos_angle) // 更新最小余弦角度
        {
            min_cos_angle = cos_angle;
            min_it = it;
        }
    }
    ftr = *min_it; // 返回最接近的特征

    // 检查观察角度是否超过60°
    if (min_cos_angle < 0.5) // 如果余弦角度小于0.5（即大于60°）
    {
        return false; // 返回false
    }

    return true; // 返回true，表示找到有效的观察特征
}
//查找具有最小评分的特征
void VisualPoint::findMinScoreFeature(const Vector3d &framepos, Feature *&ftr) const
{
    // 查找具有最小评分的特征
    auto min_it = obs_.begin(); // 初始化最小评分的迭代器
    float min_score = std::numeric_limits<float>::max(); // 初始化最小评分

    for (auto it = obs_.begin(), ite = obs_.end(); it != ite; ++it)
    {
        if ((*it)->score_ < min_score) // 如果当前特征的评分小于最小评分
        {
            min_score = (*it)->score_; // 更新最小评分
            min_it = it; // 更新最小评分特征的迭代器
        }
    }
    ftr = *min_it; // 返回最小评分的特征
}
//删除非参考补丁的特征
void VisualPoint::deleteNonRefPatchFeatures()
{
    // 删除非参考补丁的特征
    for (auto it = obs_.begin(); it != obs_.end();)
    {
        if (*it != ref_patch) // 如果当前特征不是参考补丁
        {
            delete *it; // 删除特征
            it = obs_.erase(it); // 从观察列表中移除并更新迭代器
        }
        else
        {
            ++it; // 移动到下一个特征
        }
    }
}