#pragma once

#include "../FilmTypes.h"
#include <vector>

namespace arcscope {
namespace features {

/**
 * ColorAnalyzer: Extract perceptual color features from per-second samples
 *
 * Input: Per-second CAM16-UCS observations provided by Bridge (VideoFeatureExtractor)
 * Output: ColorFeatures with CAM16-derived descriptors and state transitions
 *
 * Implements arcscope.md sections 3.4 and 4.3
 */
class ColorAnalyzer {
public:
    ColorAnalyzer() = default;

    /**
     * Analyze basic color samples to extract detailed features
     * @param samples Basic per-second color measurements
     * @return Vector of ColorFeatures with expanded metrics
     */
    std::vector<ColorFeatures> analyze(const std::vector<FeatureSample>& samples) const;

private:
    // Cluster color states using robust K-means with adaptive K, empty cluster handling, and temporal smoothing
    std::vector<int> clusterColorStates(const std::vector<FeatureSample>& samples, int K) const;
};

}  // namespace features
}  // namespace arcscope
