#include "arcscope/Normalization.h"
#include <algorithm>
#include <cmath>
#include <stdexcept>

namespace arcscope {

static constexpr double kMadEps = 1e-6; // arcscope.md 3.2: Z = (X - median) / (MAD + 1e-6)

double median(std::vector<double> v) {
    if (v.empty()) return 0.0;
    const std::size_t n = v.size();
    std::nth_element(v.begin(), v.begin() + n/2, v.end());
    double med = v[n/2];
    if (n % 2 == 0) {
        auto it = std::max_element(v.begin(), v.begin() + n/2);
        med = 0.5 * (med + *it);
    }
    return med;
}

double mad(const std::vector<double>& v, double med) {
    if (v.empty()) return 0.0;
    std::vector<double> dev(v.size());
    for (std::size_t i = 0; i < v.size(); ++i) dev[i] = std::abs(v[i] - med);
    return median(std::move(dev));
}

double quantile(std::vector<double> v, double q) {
    if (v.empty()) return 0.0;
    q = std::clamp(q, 0.0, 1.0);

    const double pos = q * (v.size() - 1);
    const std::size_t k = static_cast<std::size_t>(std::floor(pos));
    const std::size_t k2 = std::min(k + 1, v.size() - 1);
    const double frac = pos - static_cast<double>(k);

    std::nth_element(v.begin(), v.begin() + k, v.end());
    const double a = v[k];
    if (k2 == k) return a;

    std::nth_element(v.begin(), v.begin() + k2, v.end());
    const double b = v[k2];
    return a + frac * (b - a);
}

std::vector<double> robust_z(const std::vector<double>& x) {
    if (x.empty()) return {};
    const double med = median(std::vector<double>(x.begin(), x.end()));
    const double madv = mad(x, med);
    const double den = (std::isfinite(madv) ? madv : 0.0) + kMadEps;
    std::vector<double> z(x.size());
    for (std::size_t i = 0; i < x.size(); ++i) z[i] = (x[i] - med) / den;
    return z;
}

static inline double sigmoid01(double z) {
    if (z >= 40.0) return 1.0;
    if (z <= -40.0) return 0.0;
    return 1.0 / (1.0 + std::exp(-z));
}

std::vector<double> robust_sigmoid01(const std::vector<double>& x) {
    auto z = robust_z(x);
    for (auto& v : z) v = sigmoid01(v);
    return z;
}

double dot(const std::vector<double>& a, const std::vector<double>& b) {
    const std::size_t n = std::min(a.size(), b.size());
    double s = 0.0;
    for (std::size_t i = 0; i < n; ++i) s += a[i] * b[i];
    return s;
}

std::vector<double> smooth_series(const std::vector<double>& x, std::size_t window) {
    if (x.empty()) return {};
    if (window <= 1) return x;
    const std::size_t n = x.size();
    const std::size_t half = window / 2;
    std::vector<double> y(n, 0.0);
    for (std::size_t i = 0; i < n; ++i) {
        const std::size_t a = (i > half) ? (i - half) : 0;
        const std::size_t b = std::min(n - 1, i + half);
        double acc = 0.0;
        for (std::size_t j = a; j <= b; ++j) acc += x[j];
        y[i] = acc / static_cast<double>(b - a + 1);
    }
    return y;
}

// Sort by time, then deduplicate and merge (t, x) pairs with identical timestamps.
// For duplicate times, merge values by mean (supports >2 duplicates).
static void sort_and_dedup_time_series(std::vector<double>& t, std::vector<double>& x) {
    if (t.size() <= 1) return;
    if (t.size() != x.size()) throw std::runtime_error("t/x size mismatch");

    std::vector<std::pair<double, double>> pairs;
    pairs.reserve(t.size());
    for (std::size_t i = 0; i < t.size(); ++i) pairs.emplace_back(t[i], x[i]);

    std::sort(pairs.begin(), pairs.end(), [](const auto& a, const auto& b) {
        const bool a_finite = std::isfinite(a.first);
        const bool b_finite = std::isfinite(b.first);
        if (a_finite != b_finite) return a_finite; // finite times first
        return a.first < b.first;
    });

    std::vector<double> t_out;
    std::vector<double> x_out;
    t_out.reserve(pairs.size());
    x_out.reserve(pairs.size());

    t_out.push_back(pairs[0].first);
    x_out.push_back(pairs[0].second);
    int count = 1;

    for (std::size_t i = 1; i < pairs.size(); ++i) {
        if (std::abs(pairs[i].first - t_out.back()) < 1e-9) {
            x_out.back() = (x_out.back() * count + pairs[i].second) / double(count + 1);
            ++count;
        } else {
            t_out.push_back(pairs[i].first);
            x_out.push_back(pairs[i].second);
            count = 1;
        }
    }

    t = std::move(t_out);
    x = std::move(x_out);
}

std::vector<double> resample_linear(const std::vector<double>& src_t,
                                    const std::vector<double>& src_x,
                                    const std::vector<double>& dst_t) {
    if (dst_t.empty()) return {};
    if (src_t.empty() || src_x.empty() || src_t.size() != src_x.size()) {
        return std::vector<double>(dst_t.size(), 0.0);
    }

    // Sort+deduplicate timestamps to prevent interpolation issues (resampling assumes monotonic t)
    std::vector<double> t_clean = src_t;
    std::vector<double> x_clean = src_x;
    sort_and_dedup_time_series(t_clean, x_clean);

    const std::size_t n = t_clean.size();
    std::vector<double> out(dst_t.size(), 0.0);

    std::size_t j = 0;
    for (std::size_t i = 0; i < dst_t.size(); ++i) {
        const double t = dst_t[i];

        while (j + 1 < n && t_clean[j + 1] < t) ++j;

        if (t <= t_clean.front()) { out[i] = x_clean.front(); continue; }
        if (t >= t_clean.back())  { out[i] = x_clean.back();  continue; }

        const double t0 = t_clean[j];
        const double t1 = t_clean[j + 1];
        const double x0 = x_clean[j];
        const double x1 = x_clean[j + 1];
        const double w  = (t1 > t0) ? ((t - t0) / (t1 - t0)) : 0.0;
        out[i] = x0 + w * (x1 - x0);
    }
    return out;
}

std::vector<double> resample_nearest(const std::vector<double>& src_t,
                                     const std::vector<double>& src_x,
                                     const std::vector<double>& dst_t) {
    if (dst_t.empty()) return {};
    if (src_t.empty() || src_x.empty() || src_t.size() != src_x.size()) {
        return std::vector<double>(dst_t.size(), 0.0);
    }

    // Sort+deduplicate timestamps to prevent nearest-neighbor issues (resampling assumes monotonic t)
    std::vector<double> t_clean = src_t;
    std::vector<double> x_clean = src_x;
    sort_and_dedup_time_series(t_clean, x_clean);

    std::vector<double> out(dst_t.size(), 0.0);
    std::size_t j = 0;
    for (std::size_t i = 0; i < dst_t.size(); ++i) {
        const double t = dst_t[i];
        while (j + 1 < t_clean.size() && std::abs(t_clean[j + 1] - t) <= std::abs(t_clean[j] - t)) ++j;
        out[i] = x_clean[j];
    }
    return out;
}

} // namespace arcscope
