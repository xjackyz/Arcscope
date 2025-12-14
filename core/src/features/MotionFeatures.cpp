#include "arcscope/features/MotionFeatures.h"
#include "arcscope/Normalization.h"
#include <algorithm>
#include <cmath>

namespace arcscope {
namespace features {

namespace {

std::vector<double> sanitize_cut_times(const std::vector<double>& cutTimesSec, double durationSec) {
    std::vector<double> out;
    out.reserve(cutTimesSec.size());
    for (double t : cutTimesSec) {
        if (!std::isfinite(t)) continue;
        if (t <= 0.0) continue;
        if (t >= durationSec) continue;
        out.push_back(t);
    }
    std::sort(out.begin(), out.end());
    out.erase(std::unique(out.begin(), out.end(), [](double a, double b) { return std::abs(a - b) < 1e-6; }), out.end());
    return out;
}

// Fallback: derive approximate cut timestamps from per-second cut counts.
// Used only when Bridge does not provide cutTimesSec (e.g., external callers).
std::vector<double> cut_times_from_counts(const std::vector<FeatureSample>& samples) {
    std::vector<double> out;
    const std::size_t N = samples.size();
    out.reserve(N);
    for (std::size_t sec = 0; sec < N; ++sec) {
        const int cutCount = std::max(0, (int)std::llround(std::max(0.0, samples[sec].cutDensity)));
        if (cutCount <= 0) continue;
        for (int j = 0; j < cutCount; ++j) {
            const double t = (double)sec + (double)(j + 1) / (double)(cutCount + 1);
            out.push_back(t);
        }
    }
    std::sort(out.begin(), out.end());
    out.erase(std::unique(out.begin(), out.end(), [](double a, double b) { return std::abs(a - b) < 1e-6; }), out.end());
    return out;
}

} // namespace

std::vector<ShotSegment> MotionAnalyzer::detectShots(const std::vector<FeatureSample>& samples,
                                                     const std::vector<double>& cutTimesSec) const {
    std::vector<ShotSegment> shots;
    if (samples.empty()) {
        return shots;
    }

    const std::size_t N = samples.size();
    const double duration = (double)N;

    // Boundaries = [0] + cutTimes + [duration]
    std::vector<double> boundaries;
    boundaries.reserve(cutTimesSec.size() + 2);
    boundaries.push_back(0.0);

    std::vector<double> cuts = sanitize_cut_times(cutTimesSec, duration);
    if (cuts.empty()) {
        // Compatibility path: if caller doesn't provide exact cut timestamps, fall back to counts.
        cuts = sanitize_cut_times(cut_times_from_counts(samples), duration);
    }
    for (double t : cuts) boundaries.push_back(t);
    if (boundaries.back() < duration - 1e-9) boundaries.push_back(duration);

    // Build shot segments
    for (std::size_t i = 0; i + 1 < boundaries.size(); ++i) {
        ShotSegment shot;
        shot.startTime = boundaries[i];
        shot.endTime = boundaries[i + 1];
        shot.duration = shot.endTime - shot.startTime;

        // Compute avg motion by overlap with 1-second bins [k, k+1).
        double acc = 0.0;
        const double len = std::max(shot.duration, 1e-9);
        const std::size_t k0 = std::min<std::size_t>(N, (std::size_t)std::floor(shot.startTime));
        const std::size_t k1 = std::min<std::size_t>(N, (std::size_t)std::ceil(shot.endTime));
        for (std::size_t k = k0; k < k1; ++k) {
            const double a = std::max(shot.startTime, (double)k);
            const double b = std::min(shot.endTime, (double)k + 1.0);
            const double w = std::max(0.0, b - a);
            acc += w * samples[k].motionAmplitude;
        }
        shot.avgMotion = acc / len;

        shots.push_back(shot);
    }

    return shots;
}

double MotionAnalyzer::computeLocalASL(const std::vector<ShotSegment>& shots,
                                       double centerTime,
                                       double windowSize) const {
    double halfWindow = windowSize / 2.0;
    double sumDuration = 0.0;
    int count = 0;

    for (const auto& shot : shots) {
        double shotCenter = (shot.startTime + shot.endTime) / 2.0;
        if (shotCenter >= centerTime - halfWindow && shotCenter <= centerTime + halfWindow) {
            sumDuration += shot.duration;
            count++;
        }
    }

    return (count > 0) ? (sumDuration / count) : 0.0;
}

double MotionAnalyzer::computeCutDensity(const std::vector<ShotSegment>& shots,
                                         double centerTime,
                                         double windowSize) const {
    double halfWindow = windowSize / 2.0;
    int cutCount = 0;

    for (const auto& shot : shots) {
        if (shot.startTime >= centerTime - halfWindow && shot.startTime <= centerTime + halfWindow) {
            cutCount++;
        }
    }

    return static_cast<double>(cutCount) / windowSize;
}

}  // namespace features
}  // namespace arcscope
