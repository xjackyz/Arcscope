#pragma once

#include "arcscope/FilmTypes.h"

#include <vector>

namespace arcscope {

/**
 * CurveEngine: Build core curves + metadata overlays from detailed features
 *
 * Implements arcscope.md Chapter 4: 曲线系统设计
 *
 * Input: AnalysisContext with AudioFeatures, ColorFeatures, FaceFeatures, etc.
 * Output: core curves + FacePresence mask (arcscope.md 4.5.3 / 8.x persistence).
 *
 * All curves use PCA in robust Z-space to automatically determine feature weights.
 * All PCA PC1 directions are sign-aligned to reference features.
 */
class CurveEngine {
  public:
    /**
     * Build curves from analysis context
     * @param context Internal analysis context with detailed features
     * @return Vector of 6 CurveSeries
     */
    std::vector<CurveSeries> build_curves(const AnalysisContext& context) const;

  private:
    // Build individual curves
    CurveSeries buildPaceCurve(const AnalysisContext& ctx) const;
    CurveSeries buildSoundCurve(const AnalysisContext& ctx) const;
    CurveSeries buildColorCurve(const AnalysisContext& ctx) const; // ColorMood(t)
    CurveSeries buildInfoCurve(const AnalysisContext& ctx) const;
    CurveSeries buildFaceAffectCurve(const AnalysisContext& ctx) const;
    CurveSeries buildFacePresenceCurve(const AnalysisContext& ctx) const;
    // Arousal curve depends on other curves (to avoid redundant computation)
    CurveSeries buildArousalCurve(const AnalysisContext& ctx,
                                  const CurveSeries& pace,
                                  const CurveSeries& sound,
                                  const CurveSeries& color,
                                  const CurveSeries& info,
                                  const CurveSeries& faceAffect) const;
};

}  // namespace arcscope
