#pragma once

#include "arcscope/CurveEngine.h"
#include "arcscope/Diagnostics.h"
#include "arcscope/MetadataEngine.h"
#include "arcscope/SQLiteStore.h"
#include "arcscope/features/SoundFeatures.h"
#include "arcscope/features/ColorFeatures.h"
#include "arcscope/features/FaceFeatures.h"
#include "arcscope/features/MotionFeatures.h"
#include "arcscope/features/InfoFeatures.h"
#include "arcscope/segmentation/StructureAnalyzer.h"

#include <memory>
#include <string>

namespace arcscope {

/**
 * FilmEngine: Complete film analysis pipeline
 *
 * Orchestrates: Features → Curves → Diagnostics → Segmentation → Storage
 */
class FilmEngine {
  public:
    explicit FilmEngine(std::string databasePath);

    /**
     * Analyze film and store results
     * @param request Film analysis request with basic FeatureSample array
     * @return Complete analysis result
     */
    FilmAnalysisResult analyze_and_store(const FilmAnalysisRequest& request);

  private:
    // Build AnalysisContext from basic FeatureSample array and FaceTrack data
    AnalysisContext buildAnalysisContext(const std::vector<FeatureSample>& samples,
                                         const std::vector<FaceTrack>& tracks,
                                         const std::vector<double>& cutTimesSec) const;

    // Compute scalar metrics from curves
    FilmMetrics compute_metrics(const std::vector<CurveSeries>& curves) const;

    // Feature analyzers
    features::SoundAnalyzer soundAnalyzer_;
    features::ColorAnalyzer colorAnalyzer_;
    features::FaceAnalyzer faceAnalyzer_;
    features::MotionAnalyzer motionAnalyzer_;
    features::InfoAnalyzer infoAnalyzer_;

    // Segmentation
    segmentation::StructureAnalyzer structureAnalyzer_;

    // Curve and diagnostics engines
    CurveEngine curveEngine_;
    MetadataEngine metadataEngine_;
    DiagnosticsEngine diagnostics_;

    // Storage
    std::unique_ptr<SQLiteStore> store_;
};

}  // namespace arcscope
