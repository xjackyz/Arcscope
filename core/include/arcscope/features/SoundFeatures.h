#pragma once

#include "../FilmTypes.h"
#include <vector>

namespace arcscope {
namespace features {

/**
 * SoundFeatures module: Extract audio features from per-second samples
 *
 * Input: FeatureSample populated by Bridge (AudioAnalysisEngine at full sample rate)
 * Output: AudioFeatures copied to Core analysis context (no estimation in Core)
 *
 * Implements arcscope.md sections 3.4.2 and 4.2
 */
class SoundAnalyzer {
public:
    SoundAnalyzer() = default;

    /**
     * Analyze a sequence of basic audio samples to extract features
     * @param samples Basic per-second audio measurements from Bridge layer
     * @return Vector of AudioFeatures (with honest defaults for unavailable features)
     */
    std::vector<AudioFeatures> analyze(const std::vector<FeatureSample>& samples) const;
};

}  // namespace features
}  // namespace arcscope
