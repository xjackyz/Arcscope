#include "arcscope/CurveEngine.h"
#include "arcscope/Normalization.h"
#include "arcscope/CoreContract.h"

#include <algorithm>
#include <cmath>
#include <stdexcept>
#include <numeric>

namespace {

using Matrix = std::vector<std::vector<double>>;

std::vector<double> robust_z_masked(const std::vector<double>& x, const std::vector<double>& mask01) {
    if (x.empty() || mask01.empty()) return {};
    const std::size_t N = std::min(x.size(), mask01.size());
    std::vector<double> filtered;
    filtered.reserve(N);
    std::vector<std::size_t> idx;
    idx.reserve(N);
    for (std::size_t i = 0; i < N; ++i) {
        if (mask01[i] > 0.0) {
            filtered.push_back(x[i]);
            idx.push_back(i);
        }
    }
    std::vector<double> out(N, 0.0);
    if (filtered.empty()) return out;
    auto zf = arcscope::robust_z(filtered);
    for (std::size_t j = 0; j < idx.size(); ++j) out[idx[j]] = zf[j];
    return out;
}

double mean_abs_slope_1hz(const std::vector<double>& z) {
    if (z.size() < 2) return 0.0;
    double acc = 0.0;
    for (std::size_t k = 1; k < z.size(); ++k) acc += std::abs(z[k] - z[k - 1]);
    return acc / static_cast<double>(z.size() - 1);
}

// Adaptive window size computation using contract constants (doc 3.3.1)
// See: CoreContract.h for formula and rationale
double window_from_lambda(double lambda) {
    return arcscope::contract::adaptive_window_size(lambda);
}

struct ColorAffectBinding {
    bool valid{false};
    double a_w{0.0};
    double a_b{0.0};
    double a_c{0.0};
    double meanW{0.0};
    double meanB{0.0};
    double meanS{0.0};
    double meanV{0.0};
    double meanA{0.0};
};

ColorAffectBinding fit_color_affect_binding(const std::vector<double>& warmth,
                                            const std::vector<double>& brightness,
                                            const std::vector<double>& saturation,
                                            const std::vector<double>& facePresence,
                                            const std::vector<double>& faceValence,
                                            const std::vector<double>& faceArousal) {
    ColorAffectBinding out;
    const std::size_t N = std::min({warmth.size(), brightness.size(), saturation.size(),
                                    facePresence.size(), faceValence.size(), faceArousal.size()});
    if (N == 0) return out;

    // Collect masked observations (facePresence is 0/1 from FaceAnalyzer; no thresholds).
    std::vector<double> W, B, S, V, A;
    W.reserve(N); B.reserve(N); S.reserve(N); V.reserve(N); A.reserve(N);
    for (std::size_t i = 0; i < N; ++i) {
        if (facePresence[i] > 0.0) {
            W.push_back(warmth[i]);
            B.push_back(brightness[i]);
            S.push_back(saturation[i]);
            V.push_back(faceValence[i]);
            A.push_back(faceArousal[i]);
        }
    }
    if (W.empty()) return out;

    auto mean = [](const std::vector<double>& v) -> double {
        if (v.empty()) return 0.0;
        double acc = 0.0;
        for (double x : v) acc += x;
        return acc / static_cast<double>(v.size());
    };

    out.meanW = mean(W);
    out.meanB = mean(B);
    out.meanS = mean(S);
    out.meanV = mean(V);
    out.meanA = mean(A);

    // Valence model: V ~= meanV + a_w*(W-meanW) + a_b*(B-meanB)
    double Sww = 0.0, Sbb = 0.0, Swb = 0.0;
    double Svw = 0.0, Svb = 0.0;
    for (std::size_t i = 0; i < W.size(); ++i) {
        const double wc = W[i] - out.meanW;
        const double bc = B[i] - out.meanB;
        const double vc = V[i] - out.meanV;
        Sww += wc * wc;
        Sbb += bc * bc;
        Swb += wc * bc;
        Svw += vc * wc;
        Svb += vc * bc;
    }

    const double det = Sww * Sbb - Swb * Swb;
    if (std::abs(det) > 1e-12) {
        out.a_w = (Svw * Sbb - Svb * Swb) / det;
        out.a_b = (Svb * Sww - Svw * Swb) / det;
    } else {
        out.a_w = 0.0;
        out.a_b = 0.0;
    }

    // Arousal model: A ~= meanA + a_c*(S-meanS)
    double Sss = 0.0, Sas = 0.0;
    for (std::size_t i = 0; i < S.size(); ++i) {
        const double sc = S[i] - out.meanS;
        const double ac = A[i] - out.meanA;
        Sss += sc * sc;
        Sas += ac * sc;
    }
    if (Sss > 1e-12) out.a_c = Sas / Sss;

    out.valid = true;
    return out;
}

// convert a time vector of features to unified grid; if already aligned, returns raw
template <typename FeatureT, typename Getter>
std::vector<double> extract_and_resample(const std::vector<FeatureT>& feats,
                                        const std::vector<double>& dst_t,
                                        Getter get_x,
                                        bool nearest) {
    if (dst_t.empty()) return {};
    if (feats.empty()) return std::vector<double>(dst_t.size(), 0.0);

    std::vector<double> src_t(feats.size());
    std::vector<double> src_x(feats.size());
    for (std::size_t i = 0; i < feats.size(); ++i) {
        src_t[i] = feats[i].timeSeconds;
        src_x[i] = get_x(feats[i]);
    }
    if (nearest) return arcscope::resample_nearest(src_t, src_x, dst_t);
    return arcscope::resample_linear(src_t, src_x, dst_t);
}

// compute 1Hz motion from shot segments by overlap-weighted average
std::vector<double> motion_from_shots_1hz(const std::vector<arcscope::ShotSegment>& shots,
                                         const std::vector<double>& t_center) {
    std::vector<double> M(t_center.size(), 0.0);
    if (shots.empty() || t_center.empty()) return M;

    for (std::size_t k = 0; k < t_center.size(); ++k) {
        const double a = std::floor(t_center[k]);       // second start
        const double b = a + 1.0;                       // second end
        double wsum = 0.0;
        double acc  = 0.0;

        for (const auto& s : shots) {
            const double lo = std::max(a, s.startTime);
            const double hi = std::min(b, s.endTime);
            const double w  = std::max(0.0, hi - lo);
            if (w > 0.0) {
                acc  += w * s.avgMotion;
                wsum += w;
            }
        }
        M[k] = (wsum > 0.0) ? (acc / wsum) : 0.0;
    }
    return M;
}

// compute local ASL and cut density strictly from shots with overlap in window W(t)
void compute_local_asl_cd(const std::vector<arcscope::ShotSegment>& shots,
                          const std::vector<double>& t_center,
                          double windowSec,
                          std::vector<double>& ASL,
                          std::vector<double>& CD) {
    const std::size_t N = t_center.size();
    ASL.assign(N, 0.0);
    CD.assign(N, 0.0);
    if (shots.empty() || N == 0) return;

    for (std::size_t k = 0; k < N; ++k) {
        const double t = t_center[k];
        const double w0 = t - windowSec * 0.5;
        const double w1 = t + windowSec * 0.5;

        // S(t): shots overlapping window
        double dur_sum = 0.0;
        std::size_t cnt = 0;

        // cut count: number of boundaries inside window
        std::size_t cuts = 0;

        for (std::size_t i = 0; i < shots.size(); ++i) {
            const auto& s = shots[i];
            const bool overlap = !(s.endTime <= w0 || s.startTime >= w1);
            if (overlap) {
                const double d = std::max(0.0, s.endTime - s.startTime);
                dur_sum += d;
                cnt++;
            }
            // boundary at shot end (a cut) if inside window and not the last shot
            // (the last shot's endTime is film end, not a cut)
            if (i + 1 < shots.size() && s.endTime > w0 && s.endTime < w1) cuts++;
        }
        ASL[k] = (cnt > 0) ? (dur_sum / static_cast<double>(cnt)) : 0.0;
        CD[k]  = (windowSec > 1e-9) ? (static_cast<double>(cuts) / windowSec) : 0.0;
    }
}

// PCA PC1 raw score in Z-space; input columns should be RAW features (robust_z applied here)
std::vector<double> pca_pc1_raw(Matrix cols_raw, int ref_col) {
    if (cols_raw.empty()) return {};
    const std::size_t D = cols_raw.size();
    const std::size_t N = cols_raw.front().size();
    for (const auto& c : cols_raw) {
        if (c.size() != N) throw std::runtime_error("PCA input length mismatch");
    }

    // Robust standardize each column (doc 4.1.4 / 4.2.5 style)
    // This is the ONLY place where robust_z should be applied
    Matrix X;
    X.reserve(D);
    for (const auto& c : cols_raw) X.push_back(arcscope::robust_z(c));

    if (D == 1) return X[0];

    // covariance
    std::vector<double> cov(D * D, 0.0);
    const double denom = static_cast<double>((N > 1) ? (N - 1) : 1);
    for (std::size_t i = 0; i < D; ++i) {
        for (std::size_t j = i; j < D; ++j) {
            double acc = 0.0;
            for (std::size_t n = 0; n < N; ++n) acc += X[i][n] * X[j][n];
            const double v = acc / denom;
            cov[i * D + j] = v;
            cov[j * D + i] = v;
        }
    }

    // power iteration
    std::vector<double> w(D, 1.0 / std::sqrt(static_cast<double>(D)));
    for (int iter = 0; iter < 64; ++iter) {
        std::vector<double> nxt(D, 0.0);
        for (std::size_t i = 0; i < D; ++i) {
            double s = 0.0;
            for (std::size_t j = 0; j < D; ++j) s += cov[i * D + j] * w[j];
            nxt[i] = s;
        }
        double norm2 = 0.0;
        for (double v : nxt) norm2 += v * v;
        const double norm = std::sqrt(std::max(norm2, 1e-24));
        for (double& v : nxt) v /= norm;
        w.swap(nxt);
    }

    // project
    std::vector<double> score(N, 0.0);
    for (std::size_t n = 0; n < N; ++n) {
        double s = 0.0;
        for (std::size_t f = 0; f < D; ++f) s += X[f][n] * w[f];
        score[n] = s;
    }

    // sign align with reference column
    if (ref_col >= 0 && static_cast<std::size_t>(ref_col) < D) {
        const double corr_like = arcscope::dot(score, X[static_cast<std::size_t>(ref_col)]);
        if (corr_like < 0.0) for (auto& v : score) v = -v;
    }

    return score;
}

} // namespace

namespace arcscope {

std::vector<CurveSeries> CurveEngine::build_curves(const AnalysisContext& context) const {
    std::vector<CurveSeries> curves;
    curves.reserve(7);

    // Build base curves first (these are independent)
    auto pace = buildPaceCurve(context);
    auto sound = buildSoundCurve(context);
    auto color = buildColorCurve(context); // ColorMood(t)
    auto info = buildInfoCurve(context);
    auto facePresence = buildFacePresenceCurve(context);
    auto faceAffect = buildFaceAffectCurve(context);

    curves.push_back(pace);
    curves.push_back(sound);
    curves.push_back(color);
    curves.push_back(info);
    if (context.faceEnabled && !faceAffect.values.empty()) {
        curves.push_back(faceAffect);
    }
    curves.push_back(facePresence);

    // Arousal depends on the above curves (pass by reference to avoid redundant computation)
    curves.push_back(buildArousalCurve(context, pace, sound, color, info, faceAffect));

    return curves;
}

// ---------------- Pace Curve (doc 4.1 + 3.3.1) ----------------
CurveSeries CurveEngine::buildPaceCurve(const AnalysisContext& ctx) const {
    CurveSeries curve;
    curve.kind = CurveKind::Pace;
    curve.times = ctx.t_center;
    const std::size_t N = curve.times.size();
    if (N == 0 || ctx.shots.empty()) return curve;

    // motion 1hz (prefer true motion_1hz, else overlap-weighted from shots)
    std::vector<double> M = (ctx.motion_1hz.size() == N) ? ctx.motion_1hz
                                                         : motion_from_shots_1hz(ctx.shots, curve.times);

    auto compute_pace_raw_with_W = [&](double Wsec) -> std::vector<double> {
        std::vector<double> ASL, CD;
        compute_local_asl_cd(ctx.shots, curve.times, Wsec, ASL, CD);

        // Negate ASL to align direction: higher ASL => slower pace => lower raw score
        std::vector<double> neg_ASL(ASL.size());
        for (std::size_t i = 0; i < ASL.size(); ++i) neg_ASL[i] = -ASL[i];

        Matrix F;
        F.push_back(neg_ASL);  // raw -ASL
        F.push_back(CD);       // raw CD
        F.push_back(M);        // raw motion

        // align with cut density (index 1)
        return pca_pc1_raw(std::move(F), /*ref_col=*/1);
    };

    // self-consistent W iteration (doc 3.3.1)
    double W = 12.0;
    for (int iter = 0; iter < 16; ++iter) {
        auto raw = compute_pace_raw_with_W(W);
        auto z   = robust_z(raw);
        const double lambda = mean_abs_slope_1hz(z);
        const double newW = window_from_lambda(lambda);
        if (std::abs(newW - W) < 1.0) { W = newW; break; }
        W = newW;
    }

    curve.raw_score = compute_pace_raw_with_W(W);
    curve.values = robust_sigmoid01(curve.raw_score);

    // smoothing tied to lambda of raw_score
    const double lambda = mean_abs_slope_1hz(robust_z(curve.raw_score));
    const double W_s = window_from_lambda(lambda);
    const std::size_t smoothW = static_cast<std::size_t>(std::clamp(std::llround(W_s / 4.0), 1LL, 15LL));
    curve.values = smooth_series(curve.values, smoothW);

    return curve;
}

// ---------------- Sound Curve (doc 4.2) ----------------
CurveSeries CurveEngine::buildSoundCurve(const AnalysisContext& ctx) const {
    CurveSeries curve;
    curve.kind = CurveKind::Sound;
    curve.times = ctx.t_center;
    const std::size_t N = curve.times.size();
    if (N == 0 || ctx.audioFeatures.empty()) return curve;

    // resample each feature to unified 1Hz grid
    auto loud = extract_and_resample(ctx.audioFeatures, curve.times,
                                    [](const AudioFeatures& f){ return f.loudness; }, false);
    auto rhythm = extract_and_resample(ctx.audioFeatures, curve.times,
                                      [](const AudioFeatures& f){ return f.rhythmEnergy; }, false);
    auto tension = extract_and_resample(ctx.audioFeatures, curve.times,
                                       [](const AudioFeatures& f){ return f.harmonicTension; }, false);
    auto brightCentroid = extract_and_resample(ctx.audioFeatures, curve.times,
                                              [](const AudioFeatures& f){ return f.spectralCentroid; }, false);
    auto brightHFRatio = extract_and_resample(ctx.audioFeatures, curve.times,
                                             [](const AudioFeatures& f){ return f.spectralBalance; }, false);

    // SpectralBrightness(t) (doc 4.2.4): use spectral centroid or HF energy ratio.
    // Select the channel with stronger within-film variation (MAD), avoiding hard-coded thresholds.
    auto madv = [](const std::vector<double>& v) -> double {
        const double m = arcscope::median(std::vector<double>(v.begin(), v.end()));
        return arcscope::mad(v, m);
    };
    const double madCentroid = madv(brightCentroid);
    const double madHf = madv(brightHFRatio);
    const auto& bright = (madCentroid >= madHf) ? brightCentroid : brightHFRatio;

    Matrix S;
    S.push_back(loud);     // raw loudness
    S.push_back(rhythm);   // raw rhythm
    S.push_back(tension);  // raw tension
    S.push_back(bright);   // raw brightness

    curve.raw_score = pca_pc1_raw(std::move(S), /*ref_col=*/0); // align with loudness
    curve.values = robust_sigmoid01(curve.raw_score);

    const double lambda = mean_abs_slope_1hz(robust_z(curve.raw_score));
    const double W_s = window_from_lambda(lambda);
    const std::size_t smoothW = static_cast<std::size_t>(std::clamp(std::llround(W_s / 4.0), 1LL, 15LL));
    curve.values = smooth_series(curve.values, smoothW);

    return curve;
}

// ---------------- Color Curve (doc 4.3) ----------------
CurveSeries CurveEngine::buildColorCurve(const AnalysisContext& ctx) const {
    CurveSeries curve;
    curve.kind = CurveKind::Color;
    curve.times = ctx.t_center;
    const std::size_t N = curve.times.size();
    if (N == 0 || ctx.colorFeatures.empty()) return curve;

    auto warmth = extract_and_resample(ctx.colorFeatures, curve.times,
                                      [](const ColorFeatures& f){ return f.warmth; }, false);
    auto brightness = extract_and_resample(ctx.colorFeatures, curve.times,
                                          [](const ColorFeatures& f){ return f.brightness; }, false);
    auto saturation = extract_and_resample(ctx.colorFeatures, curve.times,
                                          [](const ColorFeatures& f){ return f.saturation; }, false);
    auto grayness = extract_and_resample(ctx.colorFeatures, curve.times,
                                        [](const ColorFeatures& f){ return f.grayness; }, false);
    auto harmony = extract_and_resample(ctx.colorFeatures, curve.times,
                                       [](const ColorFeatures& f){ return f.colorHarmony; }, false);
    auto contrast = extract_and_resample(ctx.colorFeatures, curve.times,
                                        [](const ColorFeatures& f){ return f.colorContrast; }, false);

    // Color–Affect Binding (doc 4.3.4): learn color->(valence, arousal) mapping per film.
    // This term is masked by FacePresence: only meaningful when dominant face is present.
    std::vector<double> align(N, 0.0);
    if (ctx.faceEnabled && !ctx.faceFeatures.empty()) {
        auto facePresence = extract_and_resample(ctx.faceFeatures, curve.times,
                                                [](const FaceFeatures& f){ return f.facePresence; }, true);
        auto faceValence = extract_and_resample(ctx.faceFeatures, curve.times,
                                               [](const FaceFeatures& f){ return f.faceValence; }, false);
        auto faceArousal = extract_and_resample(ctx.faceFeatures, curve.times,
                                               [](const FaceFeatures& f){ return f.faceArousal; }, false);

        const auto binding = fit_color_affect_binding(warmth, brightness, saturation,
                                                      facePresence, faceValence, faceArousal);
        if (binding.valid) {
            std::vector<double> vhat(N, 0.0);
            std::vector<double> ahat(N, 0.0);
            for (std::size_t k = 0; k < N; ++k) {
                vhat[k] = binding.meanV
                        + binding.a_w * (warmth[k] - binding.meanW)
                        + binding.a_b * (brightness[k] - binding.meanB);
                ahat[k] = binding.meanA
                        + binding.a_c * (saturation[k] - binding.meanS);
            }

            auto zFV = robust_z_masked(faceValence, facePresence);
            auto zFA = robust_z_masked(faceArousal, facePresence);
            auto zVh = robust_z_masked(vhat, facePresence);
            auto zAh = robust_z_masked(ahat, facePresence);

            for (std::size_t k = 0; k < N; ++k) {
                if (facePresence[k] > 0.0) {
                    const double dv = zFV[k] - zVh[k];
                    const double da = zFA[k] - zAh[k];
                    align[k] = dv * dv + da * da;
                } else {
                    align[k] = 0.0;
                }
            }
        }
    }

    // ColorMood(t): PC1 of [Warmth, J, M, Grayness, HarmonyScore, Contrast, Align] (doc 4.3.4).
    Matrix C;
    C.push_back(warmth);     // raw warmth (0..1)
    C.push_back(brightness); // raw brightness (CAM16 J)
    C.push_back(saturation); // raw saturation (CAM16 M)
    C.push_back(grayness);   // raw grayness proxy
    C.push_back(harmony);    // raw harmony score (0..1)
    C.push_back(contrast);   // raw ΔE jump (CAM16-UCS)
    if (ctx.faceEnabled) C.push_back(align); // raw Align (masked by FacePresence)

    curve.raw_score = pca_pc1_raw(std::move(C), /*ref_col=*/0); // align with warmth ("visual temperature")
    curve.values = robust_sigmoid01(curve.raw_score);

    const double lambda = mean_abs_slope_1hz(robust_z(curve.raw_score));
    const double W_s = window_from_lambda(lambda);
    const std::size_t smoothW = static_cast<std::size_t>(std::clamp(std::llround(W_s / 4.0), 1LL, 15LL));
    curve.values = smooth_series(curve.values, smoothW);

    return curve;
}

// ---------------- Info Curve (doc 4.4) ----------------
CurveSeries CurveEngine::buildInfoCurve(const AnalysisContext& ctx) const {
    CurveSeries curve;
    curve.kind = CurveKind::Info;
    curve.times = ctx.t_center;
    const std::size_t N = curve.times.size();
    if (N == 0 || ctx.infoFeatures.empty()) return curve;

    auto verbal = extract_and_resample(ctx.infoFeatures, curve.times,
                                      [](const InfoFeatures& f){ return f.verbalLoad; }, false);
    auto visual = extract_and_resample(ctx.infoFeatures, curve.times,
                                      [](const InfoFeatures& f){ return f.visualLoad; }, false);
    auto eventl = extract_and_resample(ctx.infoFeatures, curve.times,
                                     [](const InfoFeatures& f){ return f.eventLoad; }, false);
    auto expo = extract_and_resample(ctx.infoFeatures, curve.times,
                                    [](const InfoFeatures& f){ return f.exposition; }, false);

    Matrix I;
    I.push_back(verbal);  // raw verbal load
    I.push_back(visual);  // raw visual load
    I.push_back(eventl);  // raw event load
    I.push_back(expo);    // raw exposition

    curve.raw_score = pca_pc1_raw(std::move(I), /*ref_col=*/0); // align with verbal load
    curve.values = robust_sigmoid01(curve.raw_score);

    const double lambda = mean_abs_slope_1hz(robust_z(curve.raw_score));
    const double W_s = window_from_lambda(lambda);
    const std::size_t smoothW = static_cast<std::size_t>(std::clamp(std::llround(W_s / 4.0), 1LL, 15LL));
    curve.values = smooth_series(curve.values, smoothW);

    return curve;
}

// ---------------- FacePresence Curve (doc 4.5.3) ----------------
CurveSeries CurveEngine::buildFacePresenceCurve(const AnalysisContext& ctx) const {
    CurveSeries curve;
    curve.kind = CurveKind::FacePresence;
    curve.times = ctx.t_center;
    const std::size_t N = curve.times.size();

    if (N == 0) return curve;
    if (ctx.faceFeatures.empty()) {
        curve.values.assign(N, 0.0);
        return curve;
    }

    auto presence = extract_and_resample(ctx.faceFeatures, curve.times,
                                         [](const FaceFeatures& f){ return f.facePresence; }, true);
    curve.values.assign(N, 0.0);
    for (std::size_t k = 0; k < N; ++k) curve.values[k] = (presence[k] > 0.0) ? 1.0 : 0.0;

    return curve;
}

// ---------------- FaceAffect Curve (doc 4.5.3) ----------------
CurveSeries CurveEngine::buildFaceAffectCurve(const AnalysisContext& ctx) const {
    CurveSeries curve;
    curve.kind = CurveKind::FaceAffect;
    curve.times = ctx.t_center;
    const std::size_t N = curve.times.size();
    if (N == 0) return curve;

    // If face is disabled (no dominant track), the curve is absent (arcscope.md 4.5.2 persistence contract).
    if (!ctx.faceEnabled) {
        return curve;
    }
    if (ctx.faceFeatures.empty()) return curve;

    auto presence = extract_and_resample(ctx.faceFeatures, curve.times,
                                         [](const FaceFeatures& f){ return f.facePresence; }, true);
    auto expr = extract_and_resample(ctx.faceFeatures, curve.times,
                                     [](const FaceFeatures& f){ return f.expressionIntensity; }, false);
    auto ar = extract_and_resample(ctx.faceFeatures, curve.times,
                                   [](const FaceFeatures& f){ return f.faceArousal; }, false);
    auto closeUp = extract_and_resample(ctx.faceFeatures, curve.times,
                                        [](const FaceFeatures& f){ return f.closeUpWeight; }, true);

    // Robust z-score for each component (doc 4.5.3)
    auto zP = robust_z(presence);
    auto zE = robust_z(expr);
    auto zA = robust_z(ar);

    // Face_raw = w_close * (Z_P + Z_E + Z_A)
    curve.raw_score.assign(N, 0.0);
    for (std::size_t k = 0; k < N; ++k) curve.raw_score[k] = closeUp[k] * (zP[k] + zE[k] + zA[k]);

    // Mask raw_score with presence (Arousal uses FacePresence * FaceAffect; enforce here for stability)
    for (std::size_t k = 0; k < N; ++k) curve.raw_score[k] *= ((presence[k] > 0.0) ? 1.0 : 0.0);

    // FaceAffect = sigmoid(zscore(Face_raw_masked)) (doc 4.5.3)
    curve.values = robust_sigmoid01(curve.raw_score);

    // Also zero out when absent (mask enforcement)
    for (std::size_t k = 0; k < N; ++k) curve.values[k] *= ((presence[k] > 0.0) ? 1.0 : 0.0);

    const double lambda = mean_abs_slope_1hz(robust_z(curve.raw_score));
    const double W_s = window_from_lambda(lambda);
    const std::size_t smoothW = static_cast<std::size_t>(std::clamp(std::llround(W_s / 4.0), 1LL, 15LL));
    curve.values = smooth_series(curve.values, smoothW);

    return curve;
}

// ---------------- Arousal Curve (doc 4.6.3) ----------------
CurveSeries CurveEngine::buildArousalCurve(const AnalysisContext& ctx,
                                           const CurveSeries& pace,
                                           const CurveSeries& sound,
                                           const CurveSeries& color,
                                           const CurveSeries& info,
                                           const CurveSeries& faceAffect) const {
    CurveSeries curve;
    curve.kind = CurveKind::Arousal;
    curve.times = ctx.t_center;
    const std::size_t N = curve.times.size();
    if (N == 0) return curve;

    // motion 1hz
    std::vector<double> M = (ctx.motion_1hz.size() == N) ? ctx.motion_1hz
                                                         : motion_from_shots_1hz(ctx.shots, curve.times);

    std::vector<double> dM(N, 0.0);
    for (std::size_t k = 1; k < N; ++k) dM[k] = std::abs(M[k] - M[k - 1]);

    // Build PCA input matrix
    Matrix Y;
    Y.push_back(pace.raw_score);   // raw pace
    Y.push_back(sound.raw_score);  // raw sound
    Y.push_back(color.raw_score);  // raw color
    Y.push_back(info.raw_score);   // raw info

    // Only include face dimension if faceEnabled == true
    if (ctx.faceEnabled && faceAffect.raw_score.size() == N) {
        // Face dimension is already masked (raw_score *= presence)
        // This ensures presence=0 segments don't contaminate PCA with noise
        Y.push_back(faceAffect.raw_score);
    }
    // If faceEnabled == false, face dimension is completely removed from PCA

    Y.push_back(dM);  // raw motion delta

    // align sign with Sound (index 1): "high energy => high arousal"
    curve.raw_score = pca_pc1_raw(std::move(Y), /*ref_col=*/1);
    curve.values = robust_sigmoid01(curve.raw_score);

    const double lambda = mean_abs_slope_1hz(robust_z(curve.raw_score));
    const double W_s = window_from_lambda(lambda);
    const std::size_t smoothW = static_cast<std::size_t>(std::clamp(std::llround(W_s / 4.0), 1LL, 15LL));
    curve.values = smooth_series(curve.values, smoothW);

    return curve;
}

} // namespace arcscope
