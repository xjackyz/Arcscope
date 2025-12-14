#pragma once

#include "../FilmTypes.h"
#include <vector>

namespace arcscope {
namespace features {

/**
 * MotionAnalyzer: Shot detection and motion analysis
 *
 * Input:
 *   - Cut timestamps (hard cuts) detected at analysis FPS in Bridge
 *   - Per-second motion amplitude (optical flow magnitude) for shot-level avg motion
 * Output: ShotSegments (start/end/duration/avgMotion)
 *
 * Implements arcscope.md section 4.1: Pace curve foundations
 */
class MotionAnalyzer {
public:
    MotionAnalyzer() = default;

    /**
     * Build shot segments from cut timestamps.
     * @param samples Per-second feature samples (for shot avgMotion aggregation)
     * @param cutTimesSec Cut timestamps in seconds (monotonic, unique, within [0,duration])
     * @return Vector of ShotSegment
     */
    std::vector<ShotSegment> detectShots(const std::vector<FeatureSample>& samples,
                                         const std::vector<double>& cutTimesSec) const;

    /**
     * Compute local average shot length in a window
     * @param shots Detected shot boundaries
     * @param centerTime Center of analysis window
     * @param windowSize Window duration in seconds
     * @return Average shot length within window
     */
    double computeLocalASL(const std::vector<ShotSegment>& shots,
                           double centerTime,
                           double windowSize) const;

    /**
     * Compute cut density (cuts per second) in a window
     * @param shots Detected shot boundaries
     * @param centerTime Center of analysis window
     * @param windowSize Window duration in seconds
     * @return Cuts per second
     */
    double computeCutDensity(const std::vector<ShotSegment>& shots,
                             double centerTime,
                             double windowSize) const;
};

}  // namespace features
}  // namespace arcscope
