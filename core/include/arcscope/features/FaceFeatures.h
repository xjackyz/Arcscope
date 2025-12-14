#pragma once

#include "../FilmTypes.h"
#include <vector>

namespace arcscope {
namespace features {

/**
 * Face analysis output structure
 */
struct FaceAnalysisOut {
    std::vector<FaceFeatures> perSecond;        // Per-second features
    std::vector<float> facePresenceMask;        // 0/1 mask per second
    int dominantTrackId{-1};                    // -1 if none
    bool faceEnabled{false};                    // dominantTrackId >= 0
};

/**
 * FaceAnalyzer: Identify dominant face trajectory and compute per-second features
 *
 * Implements arcscope.md section 4.5: "叙事视角:只看主角"
 *
 * Key principles:
 * 1. Multi-object tracking generates candidate trajectories
 * 2. Adaptive quantile-based filtering: T_i >= Q_T(0.75) && C_i >= Q_C(0.75)
 * 3. Score_i = T_i * S_i * C_i selects dominant track
 * 4. If no track passes threshold → hasDominantFace = false, FacePresence全0
 */
class FaceAnalyzer {
public:
    FaceAnalyzer() = default;

    /**
     * Analyze face presence samples to identify dominant trajectory
     * @param samples Basic per-second face measurements from Bridge layer
     * @param tracks Pre-computed face tracks from MOT algorithm
     * @return FaceAnalysisOut with per-second features and face enabled status
     */
    FaceAnalysisOut analyze(const std::vector<FeatureSample>& samples,
                            const std::vector<FaceTrack>& tracks) const;

    /**
     * Identify dominant face track using adaptive quantile thresholds
     * Strict implementation of arcscope.md 4.5.2:
     * - T_q = Q_T(0.75), C_q = Q_C(0.75)
     * - Score_i = T_i * S_i * C_i for tracks where T_i >= T_q && C_i >= C_q
     * @param tracks All candidate face trajectories
     * @param filmDurationSec Total film duration in seconds
     * @return Index of dominant track, or -1 if none qualifies
     */
    int identifyDominantTrack(const std::vector<FaceTrack>& tracks,
                             double filmDurationSec) const;

private:
    // Compute close-up weight based on face area relative to median
    double computeCloseUpWeight(double faceArea, double medianArea) const;

    // Calculate track score: T_i * S_i * C_i
    double calculateTrackScore(const FaceTrack& track) const;
};

}  // namespace features
}  // namespace arcscope
