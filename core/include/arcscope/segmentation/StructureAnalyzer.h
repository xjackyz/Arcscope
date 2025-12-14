#pragma once

#include "../FilmTypes.h"
#include <vector>

namespace arcscope {
namespace segmentation {

/**
 * StructureAnalyzer: Multi-scale segmentation (Scene/Sequence/Act)
 *
 * Implements arcscope.md section 5: 多尺度结构分段
 *
 * Scene: 节奏-情绪耦合的变点检测
 * Sequence: 联合聚合相邻Scene
 * Act: 基于总时长的顶层分段
 */
class StructureAnalyzer {
public:
    StructureAnalyzer() = default;

    /**
     * Detect scene boundary times on the given curve time axis.
     *
     * The returned vector is strictly increasing and includes both start and end boundaries.
     * If curves are sampled at 1Hz center times (t=k+0.5), boundaries are seconds boundaries (t=k).
     *
     * This is exposed to support arcscope.md 5.6 emotion-progress reparameterization:
     * run boundary detection on reparameterized curves (u-axis), then map boundaries back to seconds.
     */
    std::vector<double> detectSceneBoundaries(const std::vector<CurveSeries>& curves,
                                              double windowSize) const;

    /**
     * Build scene segments from externally provided boundaries (in seconds).
     * Boundaries must be sorted, unique, and within the curve time domain.
     */
    std::vector<SceneSegment> scenesFromBoundaries(const std::vector<CurveSeries>& curves,
                                                   const std::vector<double>& boundaries) const;

    /**
     * Detect Scene boundaries using multi-modal change point detection
     * @param curves Pace/Sound/Color/Info/Arousal curves
     * @param windowSize Adaptive window size (from lambda)
     * @return Vector of SceneSegment
     */
    std::vector<SceneSegment> detectScenes(const std::vector<CurveSeries>& curves,
                                            double windowSize) const;

    /**
     * Group Scenes into Sequences using greedy agglomeration
     * @param scenes Scene segments
     * @param targetCount Target number of sequences (default: auto)
     * @return Vector of SequenceSegment
     */
    std::vector<SequenceSegment> groupSequences(const std::vector<SceneSegment>& scenes,
                                                 int targetCount = -1) const;

    /**
     * Group Sequences into Acts based on total duration
     * @param sequences Sequence segments
     * @param totalDuration Total film duration
     * @return Vector of ActSegment
     */
    std::vector<ActSegment> groupActs(const std::vector<SequenceSegment>& sequences,
                                       double totalDuration) const;

    /**
     * Classify Scene type based on fingerprint (fallback, uses quantiles computed internally)
     * @param scene Scene segment
     * @return SceneType enum
     */
    SceneType classifySceneType(const SceneSegment& scene) const;

private:
    // --- Quantile pack for scene labeling (arcscope.md 5.4) ---
    struct SceneQuantiles {
        double qA30 = 0.0;        // Q_A(0.3)
        double qA75 = 0.0;        // Q_A(0.75)
        double qP30 = 0.0;        // Q_P(0.3)

        double qAprime70 = 0.0;   // Q_A'(0.7)
        double qAprime30 = 0.0;   // Q_A'(0.3)
        double qAbsAprime80 = 0.0;// Q_|A'|(0.8)

        double qDPA30 = 0.0;      // Q_DPA(0.3)
        double qDPA80 = 0.0;      // Q_DPA(0.8)
    };

    // Compute one change-point score S(i) using doc 5.2 ingredients
    double computeChangePointScore(const std::vector<CurveSeries>& curves,
                                  size_t index,
                                  double windowSize) const;

    // Compute scene fingerprint (avg features)
    SceneSegment computeSceneFingerprint(const std::vector<CurveSeries>& curves,
                                          double startTime,
                                          double endTime) const;

    // Distance between two scene fingerprints
    double sceneDistance(const SceneSegment& a, const SceneSegment& b) const;

    // Merge two scene segments
    SceneSegment mergeScenes(const SceneSegment& a, const SceneSegment& b) const;

    // Quantile-driven classification (doc 5.4)
    SceneQuantiles computeSceneQuantiles(const std::vector<SceneSegment>& scenes) const;
    SceneType classifySceneType(const SceneSegment& scene, const SceneQuantiles& Q) const;

private:
    // Cache last scenes so groupActs() can compute fingerprints without changing its signature.
    mutable std::vector<SceneSegment> lastScenes_;
};

}  // namespace segmentation
}  // namespace arcscope
