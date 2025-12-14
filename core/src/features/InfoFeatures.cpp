#include "arcscope/features/InfoFeatures.h"
#include "arcscope/Normalization.h"
#include <algorithm>
#include <cmath>

namespace arcscope {
namespace features {

namespace {

double hue_entropy01_from_hist12(const std::array<double, 12>& hist) {
    double sum = 0.0;
    for (double v : hist) sum += std::max(0.0, v);
    if (sum <= 1e-12) return 0.0;
    double h = 0.0;
    for (double v : hist) {
        const double p = std::max(0.0, v) / sum;
        if (p > 1e-12) h -= p * std::log(p);
    }
    const double hmax = std::log(12.0);
    return (hmax > 1e-12) ? std::clamp(h / hmax, 0.0, 1.0) : 0.0;
}

double chroma_l1_change01(const std::vector<double>& a, const std::vector<double>& b) {
    if (a.size() != 12 || b.size() != 12) return 0.0;
    double acc = 0.0;
    for (int i = 0; i < 12; ++i) acc += std::abs(a[(std::size_t)i] - b[(std::size_t)i]);
    // L1 distance between distributions in [0,2]; map to [0,1].
    return std::clamp(0.5 * acc, 0.0, 1.0);
}

}  // namespace

std::vector<InfoFeatures> InfoAnalyzer::analyze(const std::vector<FeatureSample>& samples) const {
    // Convert FeatureSamples to InfoFeatures
    // Implements Level 0 (no NLP) / Level 1 (has NLP) dual-mode system
    std::vector<InfoFeatures> result;
    result.reserve(samples.size());

    for (std::size_t i = 0; i < samples.size(); ++i) {
        const auto& sample = samples[i];
        InfoFeatures info;
        info.timeSeconds = sample.timeSeconds;

        // DUAL MODE: Check if NLP features are available
        // Level 1 (has NLP): nlpWordCount > 0 indicates valid NLP processing
        // Level 0 (no NLP): fall back to proxy features
        const bool hasNLP = (sample.nlpWordCount > 0);

        if (hasNLP) {
            // Level 1: Use NLP features directly
            info.wordCount = sample.nlpWordCount;
            info.sentenceComplexity = sample.nlpSentenceComplexity;
            info.newConceptRatio = sample.nlpNewConceptRatio;
            // entityDensity can be used for additional analysis if needed
        } else {
            // Level 0: Use proxy features
            info.wordCount = 0;
            info.sentenceComplexity = 0.5;  // Neutral default
            info.newConceptRatio = 0.0;
        }

        // FUSION: Compute dialogue density from raw observations (audio + subtitle)
        // Proxy: If no speech detection, use audio RMS + subtitle density
        info.speechDuty = 0.5 * sample.audioRms + 0.5 * sample.subtitleDensity;

        // FUSION: Compute concept density from raw observations (subtitle + visual)
        // Avoid using intermediate fusion results - compute directly from raw observations
        double visualComplexity = 1.0 - sample.grayness;  // Color complexity as proxy
        info.conceptDensity = 0.6 * sample.subtitleDensity + 0.2 * sample.audioRms + 0.2 * visualComplexity;

        info.expositionScore = sample.expositionProbability;  // Use exposition from sample

        // Verbal auxiliary (doc 4.4.2): dialogue clarity is an input feature.
        info.dialogueClarity = std::clamp(sample.audioDialogueClarity, 0.0, 1.0);

        // Visual subfeatures (doc 4.4.3):
        // - VisualEntropy: normalized hue-distribution entropy (0..1), computed from CAM16 hue histogram.
        info.visualEntropy = hue_entropy01_from_hist12(sample.cam16HueHist);
        // - ShotScaleNumeric / FaceChange depend on dominant face track; set in FilmEngine after FaceAnalyzer runs.
        info.shotScaleNumeric = 0.0;
        info.faceChange = 0.0;
        // - CameraMotionComplexity: motion magnitude + local change proxy (still observation-level).
        const double m = std::clamp(sample.motionAmplitude, 0.0, 1.0);
        info.cameraMotionComplexity = m;

        // Event subfeatures (doc 4.4.4):
        info.soundEvents = std::clamp(sample.audioTransient, 0.0, 1.0);
        info.cutDensity = std::max(0.0, sample.cutDensity);
        // MusicChange: chroma distribution change (no fixed thresholds).
        // Note: Bridge provides per-second chroma vectors.
        info.musicChange = (i > 0) ? chroma_l1_change01(samples[i - 1].audioChroma, sample.audioChroma) : 0.0;
        // ObjectMotionJump: magnitude change in optical flow (no fixed thresholds).
        info.objectMotionJump = (i > 0) ? std::clamp(std::abs(m - std::clamp(samples[i - 1].motionAmplitude, 0.0, 1.0)), 0.0, 1.0) : 0.0;

        result.push_back(info);
    }

    return result;
}

// ========== Info Load Calculation (used by CurveEngine) ==========

namespace {

// Extract individual load components from InfoFeatures
std::vector<double> extractVerbalLoad(const std::vector<InfoFeatures>& features) {
    // Verbal load: words + speech duty + concept density + new concepts + complexity + exposition
    std::vector<double> verbalLoad;
    verbalLoad.reserve(features.size());

    for (const auto& f : features) {
        // Combine linguistic features (unnormalized)
        // FIXED: Do not divide by arbitrary constants that will flatten contributions
        double load =
            (f.wordCount / 20.0) +              // Normalize word count (20 words/sec ≈ fast speech)
            f.speechDuty +                       // Already [0,1]
            f.conceptDensity +                   // Already [0,1] from analyze(), do NOT divide by 10
            f.newConceptRatio +                  // Already [0,1]
            f.sentenceComplexity +               // Already [0,1]
            f.expositionScore +                  // Already [0,1]
            f.dialogueClarity;                   // Already [0,1] (doc 4.4.2)

        verbalLoad.push_back(load);
    }

    return verbalLoad;
}

std::vector<double> extractVisualLoad(const std::vector<InfoFeatures>& features) {
    // Visual load: entropy + shot scale difficulty + camera motion + face change
    std::vector<double> visualLoad;
    visualLoad.reserve(features.size());

    for (const auto& f : features) {
        // FIXED: Close-ups (shotScaleNumeric=1) are easier to read, so they should
        // *reduce* visual load. Use (1 - shotScaleNumeric) to represent difficulty.
        double shotScaleDifficulty = 1.0 - f.shotScaleNumeric;  // 1=wide (harder), 0=closeup (easier)

        double load =
            f.visualEntropy +                   // Already normalized
            shotScaleDifficulty +               // Wider shots = more difficult
            f.cameraMotionComplexity +          // Already normalized
            f.faceChange;                        // Rate of appearance/disappearance

        visualLoad.push_back(load);
    }

    return visualLoad;
}

std::vector<double> extractEventLoad(const std::vector<InfoFeatures>& features) {
    // Event load: sound events + music changes + motion jumps + cut density
    std::vector<double> eventLoad;
    eventLoad.reserve(features.size());

    for (const auto& f : features) {
        double load =
            f.soundEvents +                     // Non-speech sound density
            f.musicChange +                     // Musical transitions
            f.cutDensity +                      // CutDensity (doc 4.4.4)
            f.objectMotionJump;                 // Motion discontinuities

        eventLoad.push_back(load);
    }

    return eventLoad;
}

std::vector<double> extractExpositionCurve(const std::vector<InfoFeatures>& features) {
    std::vector<double> expo;
    expo.reserve(features.size());

    for (const auto& f : features) {
        expo.push_back(f.expositionScore);
    }

    return expo;
}

} // anonymous namespace

// PCA helper: standardize columns using median + MAD
std::vector<std::vector<double>> standardizeColumns(const std::vector<std::vector<double>>& columns) {
    if (columns.empty()) {
        return {};
    }
    std::vector<std::vector<double>> standardized = columns;
    for (auto& column : standardized) {
        const double med = median(column);
        const double mad_value = mad(column, med);
        if (mad_value > 1e-6) {
            for (double& v : column) {
                v = (v - med) / mad_value;
            }
        }
    }
    return standardized;
}

// PCA helper: project onto first principal component
std::vector<double> projectFirstComponent(const std::vector<std::vector<double>>& columns) {
    if (columns.empty()) {
        return {};
    }
    const size_t featureCount = columns.size();
    const size_t sampleCount = columns.front().size();

    if (featureCount == 1) {
        return robust_sigmoid01(columns.front());
    }

    auto standardized = standardizeColumns(columns);

    // Build covariance matrix
    std::vector<double> cov(featureCount * featureCount, 0.0);
    for (size_t i = 0; i < featureCount; ++i) {
        for (size_t j = i; j < featureCount; ++j) {
            double acc = 0.0;
            for (size_t n = 0; n < sampleCount; ++n) {
                acc += standardized[i][n] * standardized[j][n];
            }
            const double value = acc / static_cast<double>(sampleCount);
            cov[i * featureCount + j] = value;
            cov[j * featureCount + i] = value;
        }
    }

    // Power iteration for dominant eigenvector
    std::vector<double> eigen(featureCount, 1.0 / std::sqrt(static_cast<double>(featureCount)));
    for (int iter = 0; iter < 32; ++iter) {
        std::vector<double> next(featureCount, 0.0);
        for (size_t i = 0; i < featureCount; ++i) {
            for (size_t j = 0; j < featureCount; ++j) {
                next[i] += cov[i * featureCount + j] * eigen[j];
            }
        }
        double norm = 0.0;
        for (double v : next) {
            norm += v * v;
        }
        norm = std::sqrt(std::max(norm, 1e-12));
        for (double& v : next) {
            v /= norm;
        }
        eigen = next;
    }

    // Project data onto first PC
    std::vector<double> projected(sampleCount, 0.0);
    for (size_t n = 0; n < sampleCount; ++n) {
        double value = 0.0;
        for (size_t f = 0; f < featureCount; ++f) {
            value += standardized[f][n] * eigen[f];
        }
        projected[n] = value;
    }

    return robust_sigmoid01(projected);
}

// Build complete InfoDensity curve using PCA
// This is called by CurveEngine::buildInfoCurve()
std::vector<double> buildInfoDensityCurve(const std::vector<InfoFeatures>& infoFeatures) {
    if (infoFeatures.empty()) {
        return {};
    }

    size_t N = infoFeatures.size();

    // Extract sub-components
    auto verbalRaw = extractVerbalLoad(infoFeatures);
    auto visualRaw = extractVisualLoad(infoFeatures);
    auto eventRaw = extractEventLoad(infoFeatures);
    auto expoRaw = extractExpositionCurve(infoFeatures);

    // Normalize each component using robust sigmoid
    auto zVerbal = robust_sigmoid01(verbalRaw);
    auto zVisual = robust_sigmoid01(visualRaw);
    auto zEvent = robust_sigmoid01(eventRaw);
    auto zExpo = robust_sigmoid01(expoRaw);

    // Build feature matrix: [zVerbal, zVisual, zEvent, zExpo]
    std::vector<std::vector<double>> matrix(4);
    matrix[0] = zVerbal;
    matrix[1] = zVisual;
    matrix[2] = zEvent;
    matrix[3] = zExpo;

    // PCA: extract PC1 as InfoDensity main axis (4.4.6)
    // Use adaptive PCA weights instead of hardcoded values
    return projectFirstComponent(matrix);
}

}  // namespace features
}  // namespace arcscope
