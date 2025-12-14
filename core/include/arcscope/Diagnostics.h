#pragma once

#include "arcscope/FilmTypes.h"

#include <vector>

namespace arcscope {

class DiagnosticsEngine {
  public:
    std::vector<IssueSegment> detect(const std::vector<CurveSeries>& curves,
                                     const std::vector<MetadataOverlaySeries>& overlays) const;
};

}  // namespace arcscope
