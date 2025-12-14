#include "arcscope/segmentation/StructureAnalyzer.h"
#include "arcscope/Normalization.h"

#include <algorithm>
#include <cmath>
#include <limits>
#include <numeric>
#include <stdexcept>

namespace arcscope {
namespace segmentation {

namespace {

const CurveSeries* findCurve(const std::vector<CurveSeries>& curves, CurveKind kind) {
    for (const auto& c : curves) {
        if (c.kind == kind) return &c;
    }
    return nullptr;
}

std::vector<double> diff1(const std::vector<double>& v) {
    std::vector<double> d(v.size(), 0.0);
    for (size_t i = 1; i < v.size(); ++i) d[i] = v[i] - v[i - 1];
    return d;
}

std::vector<double> diff2(const std::vector<double>& v) {
    std::vector<double> d2(v.size(), 0.0);
    for (size_t i = 1; i + 1 < v.size(); ++i) d2[i] = v[i + 1] - 2.0 * v[i] + v[i - 1];
    return d2;
}

double mean(const std::vector<double>& v) {
    if (v.empty()) return 0.0;
    return std::accumulate(v.begin(), v.end(), 0.0) / static_cast<double>(v.size());
}

double quantile_copy(std::vector<double> v, double q) {
    if (v.empty()) return 0.0;
    q = std::clamp(q, 0.0, 1.0);
    const size_t n = v.size();
    const size_t k = static_cast<size_t>(std::llround(q * static_cast<double>(n - 1)));
    std::nth_element(v.begin(), v.begin() + k, v.end());
    return v[k];
}

double estimate_dt(const std::vector<double>& times) {
    if (times.size() < 2) return 1.0;
    std::vector<double> diffs;
    diffs.reserve(times.size() - 1);
    for (size_t i = 1; i < times.size(); ++i) diffs.push_back(times[i] - times[i - 1]);
    const double dt = arcscope::median(std::move(diffs));
    return (std::isfinite(dt) && dt > 1e-9) ? dt : 1.0;
}

// median+MAD standardize (robust z) per column
std::vector<std::vector<double>> robust_standardize_columns(const std::vector<std::vector<double>>& cols) {
    if (cols.empty()) return {};
    const size_t N = cols.front().size();
    for (const auto& c : cols) {
        if (c.size() != N) throw std::runtime_error("robust_standardize_columns: size mismatch");
    }

    std::vector<std::vector<double>> out = cols;
    for (auto& col : out) {
        const double m = arcscope::median(col);
        const double s = arcscope::mad(col, m);
        const double denom = (s > 1e-12 ? s : 1e-12);
        for (double& x : col) x = (x - m) / denom;
    }
    return out;
}

std::vector<double> pca_first_component_weights(const std::vector<std::vector<double>>& standardized_cols) {
    const size_t F = standardized_cols.size();
    const size_t N = standardized_cols.front().size();
    if (F == 0) return {};
    if (F == 1) return {1.0};

    std::vector<double> cov(F * F, 0.0);
    for (size_t i = 0; i < F; ++i) {
        for (size_t j = i; j < F; ++j) {
            double acc = 0.0;
            for (size_t k = 0; k < N; ++k) acc += standardized_cols[i][k] * standardized_cols[j][k];
            const double v = acc / static_cast<double>(std::max<size_t>(1, N - 1));
            cov[i * F + j] = v;
            cov[j * F + i] = v;
        }
    }

    std::vector<double> w(F, 1.0 / std::sqrt(static_cast<double>(F)));
    for (int it = 0; it < 48; ++it) {
        std::vector<double> nxt(F, 0.0);
        for (size_t i = 0; i < F; ++i) {
            for (size_t j = 0; j < F; ++j) nxt[i] += cov[i * F + j] * w[j];
        }
        double n2 = 0.0;
        for (double x : nxt) n2 += x * x;
        const double nrm = std::sqrt(std::max(n2, 1e-12));
        for (double& x : nxt) x /= nrm;
        w = std::move(nxt);
    }
    return w;
}

// Robust z helper: use raw_score if available (already Z-space), otherwise compute
std::vector<double> get_z_score(const CurveSeries* curve) {
    if (!curve) return {};

    // Prefer raw_score (already in Z-space from PCA)
    if (!curve->raw_score.empty()) {
        return curve->raw_score;
    }

    // Fallback: compute from values
    if (curve->values.empty()) return {};
    const double m = arcscope::median(curve->values);
    const double s = arcscope::mad(curve->values, m);
    const double denom = (s > 1e-12 ? s : 1e-12);
    std::vector<double> out(curve->values.size());
    for (size_t i = 0; i < curve->values.size(); ++i) out[i] = (curve->values[i] - m) / denom;
    return out;
}

}  // namespace


// --- 5.4 quantiles ---
StructureAnalyzer::SceneQuantiles
StructureAnalyzer::computeSceneQuantiles(const std::vector<SceneSegment>& scenes) const {
    SceneQuantiles Q;
    std::vector<double> A, P, Aprime, absAprime, DPA;
    A.reserve(scenes.size());
    P.reserve(scenes.size());
    Aprime.reserve(scenes.size());
    absAprime.reserve(scenes.size());
    DPA.reserve(scenes.size());

    for (const auto& s : scenes) {
        A.push_back(s.avgArousal);
        P.push_back(s.avgPace);
        Aprime.push_back(s.avgArousalSlope);
        absAprime.push_back(std::abs(s.avgArousalSlope));
        DPA.push_back(s.avgPaceSoundDivergence); // we store |D_PA| here (see fingerprint)
    }

    Q.qA30 = quantile_copy(A, 0.30);
    Q.qA75 = quantile_copy(A, 0.75);
    Q.qP30 = quantile_copy(P, 0.30);

    Q.qAprime70 = quantile_copy(Aprime, 0.70);
    Q.qAprime30 = quantile_copy(Aprime, 0.30);
    Q.qAbsAprime80 = quantile_copy(absAprime, 0.80);

    Q.qDPA30 = quantile_copy(DPA, 0.30);
    Q.qDPA80 = quantile_copy(DPA, 0.80);

    return Q;
}

SceneType StructureAnalyzer::classifySceneType(const SceneSegment& scene) const {
    // Fallback: if you call this without a scene-population context, we classify using a self-quantile:
    // treat this scene as "Unknown". Real classification happens in detectScenes() with quantiles.
    return SceneType::Unknown;
}

SceneType StructureAnalyzer::classifySceneType(const SceneSegment& s, const SceneQuantiles& Q) const {
    // Doc 5.4 rules (quantile-driven, no fixed 0.3/0.7/0.8 constants except quantile levels)

    // LowEnergy Basin: mean(A)<Q_A(0.3) and mean(P)<Q_P(0.3)
    if (s.avgArousal < Q.qA30 && s.avgPace < Q.qP30) {
        return SceneType::LowEnergyBasin;
    }

    // BuildUp: mean(A')>Q_A'(0.7) and D_PA low
    if (s.avgArousalSlope > Q.qAprime70 && s.avgPaceSoundDivergence < Q.qDPA30) {
        return SceneType::BuildUp;
    }

    // HighEnergy Peak: mean(A)>Q_A(0.75)
    if (s.avgArousal > Q.qA75) {
        return SceneType::HighEnergyPeak;
    }

    // Emotional Reversal: |mean(A')| > Q_|A'|(0.8)
    if (std::abs(s.avgArousalSlope) > Q.qAbsAprime80) {
        return SceneType::EmotionalReversal;
    }

    // Misaligned Tension: D_PA > Q_DPA(0.8)
    if (s.avgPaceSoundDivergence > Q.qDPA80) {
        return SceneType::MisalignedTension;
    }

    // Release: mean(A') < Q_A'(0.3)
    if (s.avgArousalSlope < Q.qAprime30) {
        return SceneType::Release;
    }

    return SceneType::Unknown;
}


// --- 5.2 change-point score (single index) ---
// Note: detectScenes() uses vectorized version for efficiency; this keeps header contract complete.
double StructureAnalyzer::computeChangePointScore(const std::vector<CurveSeries>& curves,
                                                 size_t index,
                                                 double /*windowSize*/) const {
    const auto* pace    = findCurve(curves, CurveKind::Pace);
    const auto* sound   = findCurve(curves, CurveKind::Sound);
    const auto* color   = findCurve(curves, CurveKind::Color);
    const auto* info    = findCurve(curves, CurveKind::Info);
    const auto* arousal = findCurve(curves, CurveKind::Arousal);

    if (!pace || !arousal) return 0.0;
    const size_t N = pace->values.size();
    if (index >= N) return 0.0;
    if (arousal->values.size() != N) return 0.0;
    if (sound && sound->values.size() != N) sound = nullptr;
    if (color && color->values.size() != N) color = nullptr;
    if (info  && info->values.size()  != N) info  = nullptr;

    // We compute the four features at this index (requires previous point for diffs)
    const auto zP = get_z_score(pace);
    const auto zA = get_z_score(arousal);
    const auto zS = get_z_score(sound);
    const auto zC = get_z_score(color);
    const auto zI = get_z_score(info);

    const auto A1 = diff1(arousal->values);
    const auto A2 = diff2(arousal->values);

    const size_t i = index;
    const double xdiff =
        (i == 0) ? 0.0 :
        std::sqrt(
            std::pow(zP[i] - zP[i - 1], 2) +
            std::pow(zS[i] - zS[i - 1], 2) +
            std::pow(zC[i] - zC[i - 1], 2) +
            std::pow(zI[i] - zI[i - 1], 2)
        );

    const double absA1  = std::abs(A1[i]);
    const double absA2  = std::abs(A2[i]);
    const double absDPA = std::abs(zP[i] - zA[i]);

    // For a single point we can't do PCA reliably; we return the raw sum in standardized space.
    // Full PCA weights are computed in detectScenes() (the correct place).
    const double s = xdiff + absA1 + absA2 + absDPA;
    return s;
}


// --- 5.2 scene detection (vectorized, PCA weights learned from film itself) ---
std::vector<double> StructureAnalyzer::detectSceneBoundaries(const std::vector<CurveSeries>& curves,
                                                             double windowSize) const {
    const auto* pace    = findCurve(curves, CurveKind::Pace);
    const auto* sound   = findCurve(curves, CurveKind::Sound);
    const auto* color   = findCurve(curves, CurveKind::Color);
    const auto* info    = findCurve(curves, CurveKind::Info);
    const auto* arousal = findCurve(curves, CurveKind::Arousal);

    if (!pace || !arousal || pace->values.empty() || arousal->values.empty()) return {};
    const size_t N = pace->values.size();
    if (pace->times.size() != N) return {};
    if (arousal->values.size() != N) return {};
    if (sound && sound->values.size() != N) sound = nullptr;
    if (color && color->values.size() != N) color = nullptr;
    if (info  && info->values.size()  != N) info  = nullptr;

    const double dt = estimate_dt(pace->times);
    const double axisStart = pace->times.front() - 0.5 * dt;
    const double axisEnd   = pace->times.back()  + 0.5 * dt;

    // emotion dynamics
    const auto A1 = diff1(arousal->values);
    const auto A2 = diff2(arousal->values);

    // z curves (use raw_score which is already in Z-space)
    const auto zP = get_z_score(pace);
    const auto zA = get_z_score(arousal);
    const auto zS = get_z_score(sound);
    const auto zC = get_z_score(color);
    const auto zI = get_z_score(info);

    // feature1: ||x_t - x_{t-1}||_2, x=[P,S,C,I] in robust-z space
    std::vector<double> xdiff(N, 0.0);
    for (size_t i = 1; i < N; ++i) {
        const double dp = zP[i] - zP[i - 1];
        const double ds = zS[i] - zS[i - 1];
        const double dc = zC[i] - zC[i - 1];
        const double di = zI[i] - zI[i - 1];
        xdiff[i] = std::sqrt(dp * dp + ds * ds + dc * dc + di * di);
    }

    // feature2..4: |A'|, |A''|, |D_PA|
    std::vector<double> absA1(N), absA2(N), absDPA(N);
    for (size_t i = 0; i < N; ++i) {
        absA1[i]  = std::abs(A1[i]);
        absA2[i]  = std::abs(A2[i]);
        absDPA[i] = std::abs(zP[i] - zA[i]);   // D_PA = Z_P - Z_A
    }

    // Robust standardize each feature column, then PCA weights
    auto cols = std::vector<std::vector<double>>{ xdiff, absA1, absA2, absDPA };
    cols = robust_standardize_columns(cols);
    const auto w = pca_first_component_weights(cols);

    // Score S(t) = Σ w_j * feature_j(t) (in standardized space)
    std::vector<double> S(N, 0.0);
    for (size_t i = 0; i < N; ++i) {
        S[i] = w[0] * cols[0][i] + w[1] * cols[1][i] + w[2] * cols[2][i] + w[3] * cols[3][i];
    }

    const double thr = quantile_copy(S, 0.85);
    const double minGap = std::max(dt, windowSize / 2.0);

    std::vector<double> boundaries;
    boundaries.reserve(64);
    boundaries.push_back(axisStart);

    for (size_t i = 1; i + 1 < N; ++i) {
        const bool isPeak = (S[i] > thr) && (S[i] > S[i - 1]) && (S[i] > S[i + 1]);
        if (!isPeak) continue;
        const double boundaryTime = pace->times[i] - 0.5 * dt; // convert center-aligned sample time to seconds boundary
        if (boundaryTime <= boundaries.back() + 1e-9) continue;
        if (boundaryTime - boundaries.back() >= minGap) boundaries.push_back(boundaryTime);
    }

    boundaries.push_back(axisEnd);
    return boundaries;
}

std::vector<SceneSegment> StructureAnalyzer::scenesFromBoundaries(const std::vector<CurveSeries>& curves,
                                                                 const std::vector<double>& boundaries) const {
    if (boundaries.size() < 2) return {};
    std::vector<SceneSegment> scenes;
    scenes.reserve(boundaries.size() - 1);
    for (size_t i = 0; i + 1 < boundaries.size(); ++i) {
        scenes.push_back(computeSceneFingerprint(curves, boundaries[i], boundaries[i + 1]));
    }

    lastScenes_ = scenes;
    const auto Q = computeSceneQuantiles(scenes);
    for (auto& s : scenes) s.type = classifySceneType(s, Q);
    return scenes;
}

std::vector<SceneSegment> StructureAnalyzer::detectScenes(const std::vector<CurveSeries>& curves,
                                                         double windowSize) const {
    const auto boundaries = detectSceneBoundaries(curves, windowSize);
    return scenesFromBoundaries(curves, boundaries);
}


// --- 5.3 sequences: greedy merge nearest neighbors by fingerprint distance ---
std::vector<SequenceSegment> StructureAnalyzer::groupSequences(const std::vector<SceneSegment>& scenes,
                                                              int targetCount) const {
    if (scenes.empty()) return {};

    // cache
    lastScenes_ = scenes;

    if (targetCount < 0) {
        targetCount = std::min(12, std::max(4, static_cast<int>(scenes.size()) / 3));
    }
    targetCount = std::max(1, targetCount);

    std::vector<SequenceSegment> seqs;
    seqs.reserve(scenes.size());
    for (size_t i = 0; i < scenes.size(); ++i) {
        SequenceSegment s;
        s.startTime = scenes[i].startTime;
        s.endTime   = scenes[i].endTime;
        s.sceneIndices = { static_cast<int>(i) };
        seqs.push_back(s);
    }

    auto seqFingerprint = [&](const SequenceSegment& seg) -> SceneSegment {
        SceneSegment agg;
        bool init = false;
        for (int idx : seg.sceneIndices) {
            if (!init) { agg = scenes[static_cast<size_t>(idx)]; init = true; }
            else       { agg = mergeScenes(agg, scenes[static_cast<size_t>(idx)]); }
        }
        return agg;
    };

    while (static_cast<int>(seqs.size()) > targetCount && seqs.size() > 1) {
        double best = std::numeric_limits<double>::infinity();
        size_t bestIdx = 0;

        for (size_t i = 0; i + 1 < seqs.size(); ++i) {
            const auto a = seqFingerprint(seqs[i]);
            const auto b = seqFingerprint(seqs[i + 1]);
            const double d = sceneDistance(a, b);
            if (d < best) { best = d; bestIdx = i; }
        }

        seqs[bestIdx].endTime = seqs[bestIdx + 1].endTime;
        seqs[bestIdx].sceneIndices.insert(
            seqs[bestIdx].sceneIndices.end(),
            seqs[bestIdx + 1].sceneIndices.begin(),
            seqs[bestIdx + 1].sceneIndices.end()
        );
        seqs.erase(seqs.begin() + bestIdx + 1);
    }

    return seqs;
}


// --- 5.3 acts: greedy merge nearest neighbors using cached scenes (no API change) ---
std::vector<ActSegment> StructureAnalyzer::groupActs(const std::vector<SequenceSegment>& sequences,
                                                     double totalDuration) const {
    if (sequences.empty()) return {};

    int actCount = 3;
    if (totalDuration < 60.0 * 60.0) actCount = 2;
    else if (totalDuration < 150.0 * 60.0) actCount = 3;
    else actCount = 4;
    actCount = std::max(1, actCount);

    std::vector<ActSegment> acts;
    acts.reserve(sequences.size());
    for (size_t i = 0; i < sequences.size(); ++i) {
        ActSegment a;
        a.startTime = sequences[i].startTime;
        a.endTime   = sequences[i].endTime;
        a.actNumber = static_cast<int>(i) + 1;
        a.sequenceIndices = { static_cast<int>(i) };
        acts.push_back(a);
    }

    // We need fingerprint distances. We compute Act fingerprint by merging all scenes covered by its sequences.
    // This requires lastScenes_ and SequenceSegment.sceneIndices to be valid.
    auto actFingerprint = [&](const ActSegment& act) -> SceneSegment {
        SceneSegment agg;
        bool init = false;

        for (int seqIdx : act.sequenceIndices) {
            const auto& seq = sequences[static_cast<size_t>(seqIdx)];
            for (int sceneIdx : seq.sceneIndices) {
                const auto& sc = lastScenes_[static_cast<size_t>(sceneIdx)];
                if (!init) { agg = sc; init = true; }
                else       { agg = mergeScenes(agg, sc); }
            }
        }
        return agg;
    };

    while (static_cast<int>(acts.size()) > actCount && acts.size() > 1) {
        double best = std::numeric_limits<double>::infinity();
        size_t bestIdx = 0;

        for (size_t i = 0; i + 1 < acts.size(); ++i) {
            const auto a = actFingerprint(acts[i]);
            const auto b = actFingerprint(acts[i + 1]);
            const double d = sceneDistance(a, b);
            if (d < best) { best = d; bestIdx = i; }
        }

        acts[bestIdx].endTime = acts[bestIdx + 1].endTime;
        acts[bestIdx].sequenceIndices.insert(
            acts[bestIdx].sequenceIndices.end(),
            acts[bestIdx + 1].sequenceIndices.begin(),
            acts[bestIdx + 1].sequenceIndices.end()
        );
        acts.erase(acts.begin() + bestIdx + 1);
    }

    for (size_t i = 0; i < acts.size(); ++i) acts[i].actNumber = static_cast<int>(i) + 1;
    return acts;
}


// --- fingerprint ---
SceneSegment StructureAnalyzer::computeSceneFingerprint(const std::vector<CurveSeries>& curves,
                                                       double startTime,
                                                       double endTime) const {
    SceneSegment scene;
    scene.startTime = startTime;
    scene.endTime   = endTime;

    const auto* pace    = findCurve(curves, CurveKind::Pace);
    const auto* sound   = findCurve(curves, CurveKind::Sound);
    const auto* arousal = findCurve(curves, CurveKind::Arousal);

    if (!pace || !arousal) return scene;

    const size_t N = pace->times.size();
    if (arousal->values.size() != N) return scene;

    // D_PA uses robust-z of pace and arousal
    const auto zP = get_z_score(pace);
    const auto zA = get_z_score(arousal);

    std::vector<double> Pv, Sv, Av, Aprimev, Dpav;
    Pv.reserve(256); Sv.reserve(256); Av.reserve(256); Aprimev.reserve(256); Dpav.reserve(256);

    for (size_t i = 0; i < N; ++i) {
        const double t = pace->times[i];
        if (t >= startTime && t < endTime) {
            Pv.push_back(pace->values[i]);
            if (sound && sound->values.size() == N) Sv.push_back(sound->values[i]);
            else Sv.push_back(0.0);

            Av.push_back(arousal->values[i]);
            if (i > 0) Aprimev.push_back(arousal->values[i] - arousal->values[i - 1]);

            // Store |D_PA| into avgPaceSoundDivergence (field name old, meaning updated)
            Dpav.push_back(zP[i] - zA[i]);
        }
    }

    scene.avgPace = mean(Pv);
    scene.avgSound = mean(Sv);
    scene.avgArousal = mean(Av);
    scene.avgArousalSlope = mean(Aprimev);

    scene.avgPaceSoundDivergence = std::abs(mean(Dpav)); // now means |mean(D_PA)|

    return scene;
}


// --- distance / merge ---
double StructureAnalyzer::sceneDistance(const SceneSegment& a, const SceneSegment& b) const {
    return std::sqrt(
        std::pow(a.avgPace - b.avgPace, 2) +
        std::pow(a.avgSound - b.avgSound, 2) +
        std::pow(a.avgArousal - b.avgArousal, 2) +
        std::pow(a.avgArousalSlope - b.avgArousalSlope, 2) +
        std::pow(a.avgPaceSoundDivergence - b.avgPaceSoundDivergence, 2));
}

SceneSegment StructureAnalyzer::mergeScenes(const SceneSegment& a, const SceneSegment& b) const {
    SceneSegment m;
    m.startTime = a.startTime;
    m.endTime   = b.endTime;

    const double da = std::max(1e-6, a.endTime - a.startTime);
    const double db = std::max(1e-6, b.endTime - b.startTime);
    const double wA = da / (da + db);
    const double wB = db / (da + db);

    m.avgPace = wA * a.avgPace + wB * b.avgPace;
    m.avgSound = wA * a.avgSound + wB * b.avgSound;
    m.avgArousal = wA * a.avgArousal + wB * b.avgArousal;
    m.avgArousalSlope = wA * a.avgArousalSlope + wB * b.avgArousalSlope;
    m.avgPaceSoundDivergence = wA * a.avgPaceSoundDivergence + wB * b.avgPaceSoundDivergence;

    m.type = SceneType::Unknown;
    return m;
}

}  // namespace segmentation
}  // namespace arcscope
