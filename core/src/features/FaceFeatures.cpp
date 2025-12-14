#include "arcscope/features/FaceFeatures.h"
#include "arcscope/Normalization.h"

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <limits>
#include <numeric>
#include <unordered_map>
#include <vector>

namespace arcscope {
namespace features {

namespace {

inline double clamp01(double x) {
    if (!std::isfinite(x)) return 0.0;
    return std::clamp(x, 0.0, 1.0);
}

inline double safe_div(double a, double b, double fallback = 0.0) {
    return (std::abs(b) > 1e-12) ? (a / b) : fallback;
}

// Robust-ish normalization to [0,1] using quantiles (no placeholders)
double qnorm01(double x, const std::vector<double>& xs, double q_lo, double q_hi) {
    if (xs.empty()) return clamp01(x);
    double lo = arcscope::quantile(xs, q_lo);
    double hi = arcscope::quantile(xs, q_hi);
    if (!(hi > lo + 1e-12)) return 0.5;
    return clamp01((x - lo) / (hi - lo));
}

// Sigmoid on log-ratio, stable close-up weighting
double sigmoid(double x) {
    // clamp x to avoid overflow
    x = std::clamp(x, -12.0, 12.0);
    return 1.0 / (1.0 + std::exp(-x));
}

// Build per-second index for one track.
// If multiple points fall into same second, keep the "best" point (confidence first, then area).
std::unordered_map<int, size_t> build_second_index(const FaceTrack& tr) {
    std::unordered_map<int, size_t> idx;
    idx.reserve(tr.timePoints.size());

    const bool hasConf = (!tr.confidences.empty() && tr.confidences.size() == tr.timePoints.size());
    const bool hasArea = (!tr.faceAreas.empty()   && tr.faceAreas.size()   == tr.timePoints.size());

    for (size_t i = 0; i < tr.timePoints.size(); ++i) {
        double t = tr.timePoints[i];
        if (!std::isfinite(t)) continue;
        int sec = (int)std::floor(t);
        if (sec < 0) continue;

        auto it = idx.find(sec);
        if (it == idx.end()) {
            idx[sec] = i;
            continue;
        }

        size_t cur = it->second;
        double conf_i = hasConf ? clamp01(tr.confidences[i]) : clamp01(tr.avgConfidence);
        double conf_c = hasConf ? clamp01(tr.confidences[cur]) : clamp01(tr.avgConfidence);

        double area_i = hasArea ? std::max(0.0, tr.faceAreas[i]) : std::max(0.0, tr.avgFaceArea);
        double area_c = hasArea ? std::max(0.0, tr.faceAreas[cur]) : std::max(0.0, tr.avgFaceArea);

        // Prefer higher confidence; tie-break by larger area
        if (conf_i > conf_c + 1e-6 || (std::abs(conf_i - conf_c) <= 1e-6 && area_i > area_c)) {
            idx[sec] = i;
        }
    }
    return idx;
}

// Compute a continuity score in [0,1] based on gaps between occupied seconds
double continuity_score(const FaceTrack& tr) {
    if (tr.timePoints.empty()) return 0.0;
    std::vector<int> secs;
    secs.reserve(tr.timePoints.size());
    for (double t : tr.timePoints) {
        if (!std::isfinite(t)) continue;
        secs.push_back((int)std::floor(t));
    }
    if (secs.empty()) return 0.0;
    std::sort(secs.begin(), secs.end());
    secs.erase(std::unique(secs.begin(), secs.end()), secs.end());
    if (secs.size() <= 2) return 1.0;

    int gaps = 0;
    for (size_t i = 1; i < secs.size(); ++i) {
        if (secs[i] - secs[i - 1] >= 2) gaps++;
    }
    // More gaps => lower continuity
    double g = (double)gaps / (double)(secs.size() - 1);
    return clamp01(1.0 - 1.5 * g);
}

} // namespace

FaceAnalysisOut
FaceAnalyzer::analyze(const std::vector<FeatureSample>& samples,
                      const std::vector<FaceTrack>& tracks) const {
    FaceAnalysisOut out;
    out.perSecond.reserve(samples.size());
    out.facePresenceMask.reserve(samples.size());

    // Calculate film duration from samples
    double filmDurationSec = 0.0;
    if (!samples.empty()) {
        // timeSeconds is center time (k+0.5), so film duration ≈ last center + 0.5
        filmDurationSec = samples.back().timeSeconds + 0.5;
    }

    // Identify dominant track
    out.dominantTrackId = identifyDominantTrack(tracks, filmDurationSec);
    out.faceEnabled = (out.dominantTrackId >= 0 && out.dominantTrackId < (int)tracks.size());

    // Build per-second lookup for dominant track
    std::unordered_map<int, size_t> domIdx;
    if (out.faceEnabled) {
        domIdx = build_second_index(tracks[out.dominantTrackId]);
    }

    // Calculate median face area for close-up weighting (only if face enabled)
    double medianArea = 1.0;
    if (out.faceEnabled) {
        std::vector<double> areas;
        areas.reserve(1024);

        const auto& tr = tracks[out.dominantTrackId];
        if (tr.faceAreas.size() == tr.timePoints.size() && !tr.faceAreas.empty()) {
            for (double a : tr.faceAreas) {
                if (std::isfinite(a) && a > 1e-9) areas.push_back(a);
            }
        } else if (std::isfinite(tr.avgFaceArea) && tr.avgFaceArea > 1e-9) {
            for (int k = 0; k < 8; ++k) areas.push_back(tr.avgFaceArea);
        }

        if (!areas.empty()) {
            medianArea = arcscope::median(areas);
            if (!(medianArea > 1e-9)) medianArea = 1.0;
        }
    }

    // Build per-second features
    for (const auto& s : samples) {
        FaceFeatures ff;
        ff.timeSeconds = s.timeSeconds;
        ff.dominantTrackId = out.dominantTrackId;

        const int sec = (int)std::floor(s.timeSeconds);

        if (out.faceEnabled) {
            // Dominant track exists
            const auto& tr = tracks[out.dominantTrackId];
            auto it = domIdx.find(sec);

            if (it != domIdx.end()) {
                // Dominant face present in this second
                const size_t i = it->second;

                ff.facePresence = 1.0;

                // Face area
                if (tr.faceAreas.size() == tr.timePoints.size() && i < tr.faceAreas.size()) {
                    ff.faceArea = std::max(0.0, tr.faceAreas[i]);
                } else {
                    ff.faceArea = std::max(0.0, tr.avgFaceArea);
                }

                // Valence
                if (tr.valences.size() == tr.timePoints.size() && i < tr.valences.size()) {
                    ff.faceValence = std::clamp(tr.valences[i], -1.0, 1.0);
                } else {
                    ff.faceValence = 0.0;
                }

                // Arousal
                if (tr.arousals.size() == tr.timePoints.size() && i < tr.arousals.size()) {
                    ff.faceArousal = clamp01(tr.arousals[i]);
                } else {
                    ff.faceArousal = 0.0;
                }

                // Expression intensity: E = |V| (no fixed weights)
                // This is observation-level; final combination happens in CurveEngine
                ff.expressionIntensity = std::abs(ff.faceValence);

                // Close-up weight: document formula w_close = clip(S / (median(S) + 1e-6), 0, 1)
                ff.closeUpWeight = clamp01(ff.faceArea / (medianArea + 1e-6));

                // Head pose + downcast proxy (observation-level, for metadata overlays)
                if (tr.yaws.size() == tr.timePoints.size() && i < tr.yaws.size()) ff.yaw = tr.yaws[i];
                if (tr.rolls.size() == tr.timePoints.size() && i < tr.rolls.size()) ff.roll = tr.rolls[i];
                if (tr.downcastRatios.size() == tr.timePoints.size() && i < tr.downcastRatios.size()) {
                    ff.downcastRatio = std::clamp(tr.downcastRatios[i], 0.0, 1.0);
                } else {
                    ff.downcastRatio = 0.5;
                }

            } else {
                // Dominant face absent this second
                ff.facePresence = 0.0;
                ff.faceArea = 0.0;
                ff.faceValence = 0.0;
                ff.faceArousal = 0.0;
                ff.expressionIntensity = 0.0;
                ff.closeUpWeight = 0.0;
                ff.yaw = 0.0;
                ff.roll = 0.0;
                ff.downcastRatio = 0.5;
            }
        } else {
            // No reliable dominant track: entire film disabled
            // All seconds output zeros as per design document
            ff.facePresence = 0.0;
            ff.faceArea = 0.0;
            ff.faceValence = 0.0;
            ff.faceArousal = 0.0;
            ff.expressionIntensity = 0.0;
            ff.closeUpWeight = 0.0;
            ff.yaw = 0.0;
            ff.roll = 0.0;
            ff.downcastRatio = 0.5;
        }

        out.perSecond.push_back(std::move(ff));
        out.facePresenceMask.push_back(ff.facePresence);
    }

    return out;
}

int FaceAnalyzer::identifyDominantTrack(const std::vector<FaceTrack>& tracks,
                                        double filmDurationSec) const {
    if (tracks.empty() || !(filmDurationSec > 0.0)) return -1;

    // Step 1: Compute T_i, S_i, C_i for each track
    std::vector<double> T, S, C;
    T.reserve(tracks.size());
    S.reserve(tracks.size());
    C.reserve(tracks.size());

    for (const auto& tr : tracks) {
        // T_i = appearanceDuration / filmDuration, clamped to [0,1]
        // Note: appearanceDuration is now count of unique seconds (from Bridge)
        double Ti = std::clamp(tr.appearanceDuration / filmDurationSec, 0.0, 1.0);

        // S_i = avgFaceArea (already a ratio: face_area / frame_area)
        double Si = std::max(0.0, tr.avgFaceArea);

        // C_i = avgConfidence, clamped to [0,1]
        double Ci = clamp01(tr.avgConfidence);

        T.push_back(Ti);
        S.push_back(Si);
        C.push_back(Ci);
    }

    // Step 2: Compute adaptive thresholds Q_T(0.75) and Q_C(0.75)
    double Tq = arcscope::quantile(T, 0.75);
    double Cq = arcscope::quantile(C, 0.75);

    // Step 3: Filter candidates and find max score
    int best = -1;
    double bestScore = -1.0;

    for (int i = 0; i < (int)tracks.size(); ++i) {
        double Ti = T[i];
        double Si = S[i];
        double Ci = C[i];

        // Only consider tracks that pass both thresholds
        if (!(Ti >= Tq && Ci >= Cq)) continue;

        // Score_i = T_i * S_i * C_i
        double score = Ti * Si * Ci;

        if (score > bestScore) {
            bestScore = score;
            best = i;
        }
    }

    // Step 4: Fallback if no candidates passed thresholds
    // Instead of returning -1 immediately, find the track with max composite score
    if (best == -1) {
        for (int i = 0; i < (int)tracks.size(); ++i) {
            double Ti = T[i];
            double Si = S[i];
            double Ci = C[i];
            double score = Ti * Si * Ci;

            if (score > bestScore) {
                bestScore = score;
                best = i;
            }
        }
    }

    // Step 5: Sanity gate - ensure minimum Ti and Ci thresholds
    // If the best track is still unreliable (too short appearance or too low confidence), disable face
    if (best >= 0) {
        double bestTi = T[best];
        double bestCi = C[best];

        // Minimum thresholds: at least 1% of film duration and 15% confidence
        if (bestTi < 0.01 || bestCi < 0.15) {
            return -1;
        }
    }

    // If best == -1 at this point, no track had any score → entire film disabled
    return best;
}

double FaceAnalyzer::computeCloseUpWeight(double faceArea, double medianArea) const {
    // Document formula: w_close = clip(S / (median(S) + 1e-6), 0, 1)
    if (!(faceArea > 0.0)) return 0.0;
    return clamp01(faceArea / (medianArea + 1e-6));
}

double FaceAnalyzer::calculateTrackScore(const FaceTrack& track) const {
    // Kept for API compatibility; still valid but not used as the only criterion anymore.
    return std::max(0.0, track.appearanceDuration) *
           std::max(0.0, track.avgFaceArea) *
           clamp01(track.avgConfidence);
}

}  // namespace features
}  // namespace arcscope
