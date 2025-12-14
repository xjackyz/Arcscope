#pragma once

#include "arcscope/FilmTypes.h"

#include <vector>

namespace arcscope {

/**
 * MetadataEngine: Build metadata overlays (arcscope.md 5.6 / 10.7)
 *
 * Overlays are display-only tracks that NEVER enter PCA/Arousal/Diagnostics.
 * They are stored and loaded separately from core curves.
 */
class MetadataEngine {
  public:
    std::vector<MetadataOverlaySeries> build_overlays(const AnalysisContext& ctx,
                                                      const std::vector<CurveSeries>& curves) const;
};

}  // namespace arcscope
