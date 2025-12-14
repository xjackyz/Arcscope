#pragma once

#include "SubtitleCue.h"
#include <vector>

namespace arcscope {
namespace nlp {

/**
 * Per-second NLP features extracted from subtitle text
 * All values are 1Hz center-aligned (k+0.5 seconds)
 */
struct NLPSecondFeatures {
    double timeSeconds{0.0};              // Center time (k+0.5)
    int wordCount{0};                     // True word count from tokenizer
    double sentenceComplexity{0.0};       // [0,1] sentence complexity score
    double newConceptRatio{0.0};          // [0,1] ratio of new concepts
    double entityDensity{0.0};            // [0,1] named entity density (optional)
};

/**
 * Subtitle NLP Processor
 *
 * Level 0 (no NLP): Returns empty vector, caller uses fallback proxy
 * Level 1 (has NLP): Returns per-second NLP features
 *
 * Input: SubtitleCues with text_utf8 field populated
 * Output: Per-second NLPSecondFeatures aligned to 1Hz grid
 */
class SubtitleNLPProcessor {
public:
    SubtitleNLPProcessor();
    ~SubtitleNLPProcessor();

    /**
     * Process subtitle cues and generate per-second NLP features
     *
     * @param cues Subtitle cues with text_utf8 field populated
     * @param durationSeconds Total video duration
     * @return Per-second NLP features (length = ceil(durationSeconds))
     *         Returns empty vector if no valid text or NLP unavailable
     */
    std::vector<NLPSecondFeatures> process(
        const std::vector<arcscope::bridge::SubtitleCue>& cues,
        double durationSeconds
    );

private:
    class Impl;
    Impl* impl_{nullptr};
};

}  // namespace nlp
}  // namespace arcscope
