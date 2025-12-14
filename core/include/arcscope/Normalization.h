#pragma once
#include <vector>
#include <cstddef>

namespace arcscope {

// ---------- robust statistics ----------
double median(std::vector<double> v);                 // copy by value
double mad(const std::vector<double>& v, double med); // median(|x-med|)
double quantile(std::vector<double> values, double q); // quantile calculation

// robust z-score: (x-med)/(MAD+eps), eps 内部处理
std::vector<double> robust_z(const std::vector<double>& x);

// robust sigmoid normalize: sigmoid( robust_z(x) )
std::vector<double> robust_sigmoid01(const std::vector<double>& x);

// ---------- utilities ----------
double dot(const std::vector<double>& a, const std::vector<double>& b);

// moving average smoothing (centered, edge-clamped)
std::vector<double> smooth_series(const std::vector<double>& x, std::size_t window);

// linear interpolation to a target time grid (for continuous features)
std::vector<double> resample_linear(const std::vector<double>& src_t,
                                    const std::vector<double>& src_x,
                                    const std::vector<double>& dst_t);

// nearest neighbor resample (for masks / discrete)
std::vector<double> resample_nearest(const std::vector<double>& src_t,
                                     const std::vector<double>& src_x,
                                     const std::vector<double>& dst_t);

} // namespace arcscope
