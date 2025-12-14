#pragma once

#include "../FilmTypes.h"
#include <vector>

namespace arcscope {
namespace features {

/**
 * InfoAnalyzer: Information density and cognitive load analysis
 *
 * Input: Subtitle density, exposition probability, dialogue density
 * Output: Verbal/Visual/Event load metrics
 *
 * Implements arcscope.md section 4.4
 *
 * Note: Full NLP pipeline (tokenization, NER, concept tracking) is
 * implemented in Bridge layer using macOS NaturalLanguage.
 * This analyzer works with pre-extracted InfoFeatures.
 */
class InfoAnalyzer {
public:
    InfoAnalyzer() = default;

    /**
     * Analyze info density samples
     * Converts basic FeatureSamples to detailed InfoFeatures.
     * Bridge layer will populate additional NLP features separately.
     *
     * @param samples Per-second feature samples
     * @return InfoFeatures with basic metrics populated
     */
    std::vector<InfoFeatures> analyze(const std::vector<FeatureSample>& samples) const;
};

/**
 * Build InfoDensity curve from InfoFeatures using PCA
 *
 * Combines VerbalLoad, VisualLoad, EventLoad, and Exposition
 * into a unified info density measure via PCA or weighted average.
 *
 * Implements arcscope.md section 4.4.6
 *
 * @param infoFeatures Per-second info features from NLP pipeline
 * @return Normalized [0,1] InfoDensity curve
 */
std::vector<double> buildInfoDensityCurve(const std::vector<InfoFeatures>& infoFeatures);

}  // namespace features
}  // namespace arcscope
