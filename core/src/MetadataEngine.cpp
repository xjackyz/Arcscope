#include "arcscope/MetadataEngine.h"
#include "arcscope/Normalization.h"
#include "arcscope/CoreContract.h"

#include <algorithm>
#include <array>
#include <cmath>
#include <numeric>
#include <stdexcept>

namespace {

using Matrix = std::vector<std::vector<double>>;

double sigmoid01(double z);

std::vector<double> robust_sigmoid01_masked(const std::vector<double>& x, const std::vector<double>& mask01) {
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
    auto z = arcscope::robust_z(filtered);
    for (std::size_t j = 0; j < z.size(); ++j) z[j] = sigmoid01(z[j]);
    for (std::size_t j = 0; j < idx.size(); ++j) out[idx[j]] = z[j];
    return out;
}

std::vector<double> robust_z_masked(const std::vector<double>& x, const std::vector<double>& mask01) {
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
    auto z = arcscope::robust_z(filtered);
    for (std::size_t j = 0; j < idx.size(); ++j) out[idx[j]] = z[j];
    return out;
}

double mean_abs_slope_1hz(const std::vector<double>& z) {
    if (z.size() < 2) return 0.0;
    double acc = 0.0;
    for (std::size_t k = 1; k < z.size(); ++k) acc += std::abs(z[k] - z[k - 1]);
    return acc / static_cast<double>(z.size() - 1);
}

double window_from_lambda(double lambda) {
    return arcscope::contract::adaptive_window_size(lambda);
}

double median_stable_duration_from_mask01(const std::vector<double>& mask01, double dtSeconds) {
    if (mask01.empty()) return 0.0;
    std::vector<double> lensSec;
    std::size_t cur = 0;
    for (double x : mask01) {
        if (x > 0.0) {
            cur++;
        } else if (cur > 0) {
            lensSec.push_back(static_cast<double>(cur) * dtSeconds);
            cur = 0;
        }
    }
    if (cur > 0) lensSec.push_back(static_cast<double>(cur) * dtSeconds);
    if (lensSec.empty()) return 0.0;
    return arcscope::median(std::move(lensSec));
}

// arcscope.md 6.2 helper: median length of segments where x > Q_x(qLevel)
double stable_duration_above_quantile(const std::vector<double>& x, double qLevel, double dtSeconds) {
    if (x.empty()) return 0.0;
    const double thr = arcscope::quantile(std::vector<double>(x.begin(), x.end()), qLevel);
    std::vector<double> mask01(x.size(), 0.0);
    for (std::size_t i = 0; i < x.size(); ++i) mask01[i] = (x[i] > thr) ? 1.0 : 0.0;
    return median_stable_duration_from_mask01(mask01, dtSeconds);
}

double sigmoid01(double z) {
    if (z >= 40.0) return 1.0;
    if (z <= -40.0) return 0.0;
    return 1.0 / (1.0 + std::exp(-z));
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

std::vector<double> pca_pc1_raw(Matrix cols_raw, int ref_col) {
    if (cols_raw.empty()) return {};
    const std::size_t D = cols_raw.size();
    const std::size_t N = cols_raw.front().size();
    for (const auto& c : cols_raw) {
        if (c.size() != N) throw std::runtime_error("PCA input length mismatch");
    }

    Matrix X;
    X.reserve(D);
    for (const auto& c : cols_raw) X.push_back(arcscope::robust_z(c));
    if (D == 1) return X[0];

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

    std::vector<double> score(N, 0.0);
    for (std::size_t n = 0; n < N; ++n) {
        double s = 0.0;
        for (std::size_t f = 0; f < D; ++f) s += X[f][n] * w[f];
        score[n] = s;
    }

    if (ref_col >= 0 && static_cast<std::size_t>(ref_col) < D) {
        const double corr_like = arcscope::dot(score, X[static_cast<std::size_t>(ref_col)]);
        if (corr_like < 0.0) for (auto& v : score) v = -v;
    }

    return score;
}

double circular_mean_hue_deg_from_hist12(const std::vector<double>& hist12) {
    if (hist12.size() != 12) return 0.0;
    double sx = 0.0, sy = 0.0, wsum = 0.0;
    for (int b = 0; b < 12; ++b) {
        const double w = std::max(0.0, hist12[(std::size_t)b]);
        // bin center: (b*30 + 15) degrees
        const double ang = (double)b * 30.0 + 15.0;
        const double rad = ang * M_PI / 180.0;
        sx += w * std::cos(rad);
        sy += w * std::sin(rad);
        wsum += w;
    }
    if (!(wsum > 1e-12)) return 0.0;
    const double ang = std::atan2(sy, sx);
    double deg = ang * 180.0 / M_PI;
    if (deg < 0.0) deg += 360.0;
    return deg;
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
    }

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

} // namespace

namespace arcscope {

std::vector<MetadataOverlaySeries> MetadataEngine::build_overlays(const AnalysisContext& ctx,
                                                                  const std::vector<CurveSeries>& curves) const {
    std::vector<MetadataOverlaySeries> out;
    const std::vector<double>& t = ctx.t_center;
    const std::size_t N = t.size();
    if (N == 0) return out;

    auto add_f32 = [&](const std::string& type, const std::vector<double>& values, bool smooth, bool normalize01) {
        MetadataOverlaySeries s;
        s.type = type;
        s.fps = 1.0;
        s.channels = 1;
        s.valueType = OverlayValueType::Float32;

        std::vector<double> v = values;
        if (v.size() != N) v.resize(N, 0.0);

        if (normalize01) {
            const auto z = arcscope::robust_z(v);
            v.resize(z.size());
            for (std::size_t i = 0; i < z.size(); ++i) v[i] = sigmoid01(z[i]);
        }

        if (smooth) {
            const double lambda = mean_abs_slope_1hz(arcscope::robust_z(v));
            const double W_s = window_from_lambda(lambda);
            const std::size_t smoothW = static_cast<std::size_t>(std::clamp(std::llround(W_s / 4.0), 1LL, 15LL));
            v = arcscope::smooth_series(v, smoothW);
        }

        s.valuesF32.assign(v.begin(), v.end());
        out.push_back(std::move(s));
    };

    auto add_i32 = [&](const std::string& type, const std::vector<int32_t>& values) {
        MetadataOverlaySeries s;
        s.type = type;
        s.fps = 1.0;
        s.channels = 1;
        s.valueType = OverlayValueType::Int32;
        s.valuesI32 = values;
        if (s.valuesI32.size() != N) s.valuesI32.resize(N, 0);
        out.push_back(std::move(s));
    };

    // ---- Color overlays (doc 4.3 + 5.6/10.7) ----
    if (!ctx.colorFeatures.empty()) {
        auto warmth = extract_and_resample(ctx.colorFeatures, t, [](const ColorFeatures& f){ return f.warmth; }, false);
        add_f32("color_warmth", warmth, /*smooth=*/true, /*normalize01=*/true);

        auto brightness = extract_and_resample(ctx.colorFeatures, t, [](const ColorFeatures& f){ return f.brightness; }, false);
        auto saturation = extract_and_resample(ctx.colorFeatures, t, [](const ColorFeatures& f){ return f.saturation; }, false);
        add_f32("color_brightness", brightness, /*smooth=*/true, /*normalize01=*/true);
        add_f32("color_saturation", saturation, /*smooth=*/true, /*normalize01=*/true);

        auto grayness = extract_and_resample(ctx.colorFeatures, t, [](const ColorFeatures& f){ return f.grayness; }, false);
        add_f32("color_grayness", grayness, /*smooth=*/true, /*normalize01=*/true);

        Matrix E;
        E.push_back(brightness);
        E.push_back(saturation);
        auto energy_raw = pca_pc1_raw(std::move(E), /*ref_col=*/0);
        add_f32("color_energy", energy_raw, /*smooth=*/true, /*normalize01=*/true);

        auto harmony = extract_and_resample(ctx.colorFeatures, t, [](const ColorFeatures& f){ return f.colorHarmony; }, false);
        add_f32("color_harmony", harmony, /*smooth=*/true, /*normalize01=*/false);

        auto contrast = extract_and_resample(ctx.colorFeatures, t, [](const ColorFeatures& f){ return f.colorContrast; }, false);
        add_f32("color_contrast", contrast, /*smooth=*/true, /*normalize01=*/false);

        // Hue histogram: multi-channel 12-bin distribution per second.
        MetadataOverlaySeries hueHist;
        hueHist.type = "hue_hist12";
        hueHist.fps = 1.0;
        hueHist.channels = 12;
        hueHist.valueType = OverlayValueType::Float32;
        hueHist.valuesF32.reserve(N * 12);
        for (std::size_t i = 0; i < N; ++i) {
            if (i < ctx.colorFeatures.size() && ctx.colorFeatures[i].hueDistribution.size() == 12) {
                for (int b = 0; b < 12; ++b) hueHist.valuesF32.push_back((float)ctx.colorFeatures[i].hueDistribution[(std::size_t)b]);
            } else {
                for (int b = 0; b < 12; ++b) hueHist.valuesF32.push_back(0.0f);
            }
        }
        out.push_back(std::move(hueHist));

        // Dominant hue (degrees) derived from hue_hist12 for easy UI usage.
        std::vector<double> hueDeg(N, 0.0);
        for (std::size_t i = 0; i < N; ++i) {
            if (i < ctx.colorFeatures.size() && ctx.colorFeatures[i].hueDistribution.size() == 12) {
                std::vector<double> h(ctx.colorFeatures[i].hueDistribution.begin(), ctx.colorFeatures[i].hueDistribution.end());
                hueDeg[i] = circular_mean_hue_deg_from_hist12(h);
            }
        }
        add_f32("hue_deg", hueDeg, /*smooth=*/true, /*normalize01=*/false);

        std::vector<int32_t> colorStateRaw(N, 0);
        for (std::size_t i = 0; i < N; ++i) {
            if (i < ctx.colorFeatures.size()) colorStateRaw[i] = (int32_t)ctx.colorFeatures[i].colorState;
        }

        // Debounce/de-jitter: if a single-second spike differs from both neighbors but neighbors agree, snap it.
        std::vector<int32_t> colorState = colorStateRaw;
        if (N >= 3) {
            for (std::size_t i = 1; i + 1 < N; ++i) {
                const int32_t a = colorState[i - 1];
                const int32_t b = colorState[i];
                const int32_t c = colorState[i + 1];
                if (a == c && b != a) colorState[i] = a;
            }
        }
        add_i32("color_state", colorState);

        // StateShift(t): boundary ticks derived from debounced state (doc 4.3.5)
        std::vector<int32_t> stateShift(N, 0);
        for (std::size_t k = 1; k < N; ++k) stateShift[k] = (colorState[k] != colorState[k - 1]) ? 1 : 0;
        add_i32("state_shift", stateShift);

        // AlignNorm + Counterpoint (doc 4.3.4 / 6.3.8)
        if (ctx.faceEnabled && !ctx.faceFeatures.empty()) {
            auto facePresence = extract_and_resample(ctx.faceFeatures, t, [](const FaceFeatures& f){ return f.facePresence; }, true);
            auto faceValence = extract_and_resample(ctx.faceFeatures, t, [](const FaceFeatures& f){ return f.faceValence; }, false);
            auto faceArousal = extract_and_resample(ctx.faceFeatures, t, [](const FaceFeatures& f){ return f.faceArousal; }, false);

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

                const auto zFV = robust_z_masked(faceValence, facePresence);
                const auto zFA = robust_z_masked(faceArousal, facePresence);
                const auto zVh = robust_z_masked(vhat, facePresence);
                const auto zAh = robust_z_masked(ahat, facePresence);

                std::vector<double> align(N, 0.0);
                for (std::size_t k = 0; k < N; ++k) {
                    if (facePresence[k] > 0.0) {
                        const double dv = zFV[k] - zVh[k];
                        const double da = zFA[k] - zAh[k];
                        align[k] = dv * dv + da * da;
                    } else {
                        align[k] = 0.0;
                    }
                }

                // AlignNorm(k): robust z + sigmoid, computed ONLY on face-present seconds (mask)
                auto alignNorm = robust_sigmoid01_masked(align, facePresence);
                {
                    MetadataOverlaySeries s;
                    s.type = "align_norm";
                    s.fps = 1.0;
                    s.channels = 1;
                    s.valueType = OverlayValueType::Float32;
                    s.valuesF32.assign(alignNorm.begin(), alignNorm.end());
                    out.push_back(std::move(s));
                }

                // Counterpoint: AlignNorm > Q(0.8) on face-present seconds (adaptive, no fixed thresholds)
                std::vector<double> presentAlign;
                presentAlign.reserve(N);
                for (std::size_t k = 0; k < N; ++k) {
                    if (facePresence[k] > 0.0) presentAlign.push_back(alignNorm[k]);
                }
                double thr = presentAlign.empty() ? 1.0 : arcscope::quantile(presentAlign, 0.8);
                std::vector<int32_t> counterpoint(N, 0);
                for (std::size_t k = 0; k < N; ++k) {
                    if (facePresence[k] > 0.0 && alignNorm[k] > thr) counterpoint[k] = 1;
                }
                add_i32("counterpoint", counterpoint);
            }
        }
    }

    // ---- Emotional progress u(t) (doc 4.6.5) ----
    {
        const CurveSeries* arousal = nullptr;
        for (const auto& c : curves) {
            if (c.kind == CurveKind::Arousal) { arousal = &c; break; }
        }
        if (arousal && arousal->values.size() == N) {
            std::vector<double> u(N, 0.0);
            double total = 0.0;
            for (double a : arousal->values) total += std::max(0.0, a);
            if (!(total > 1e-12)) {
                // If Arousal is all zeros, fall back to linear progress to keep u monotonic and usable for alignment.
                for (std::size_t k = 0; k < N; ++k) u[k] = (N > 1) ? (double)k / (double)(N - 1) : 1.0;
            } else {
                double acc = 0.0;
                for (std::size_t k = 0; k < N; ++k) {
                    acc += std::max(0.0, arousal->values[k]);
                    u[k] = std::clamp(acc / total, 0.0, 1.0);
                }
                u.front() = 0.0;
                u.back() = 1.0;
            }
            add_f32("emotion_progress", u, /*smooth=*/false, /*normalize01=*/false);
        }
    }

    // ---- Audio overlays (doc 4.2 + 5.6/10.7) ----
    if (!ctx.audioFeatures.empty()) {
        auto loudness = extract_and_resample(ctx.audioFeatures, t, [](const AudioFeatures& f){ return f.loudness; }, false);
        add_f32("loudness", loudness, /*smooth=*/true, /*normalize01=*/false);

        auto rhythmDrive = extract_and_resample(ctx.audioFeatures, t, [](const AudioFeatures& f){ return f.rhythmEnergy; }, false);
        add_f32("rhythm_drive", rhythmDrive, /*smooth=*/true, /*normalize01=*/false);

        // SpectralBrightness(t) (doc 4.2.4): expose both centroid (Hz) and HF ratio; use centroid for the main overlay.
        auto spectralCentroid = extract_and_resample(ctx.audioFeatures, t, [](const AudioFeatures& f){ return f.spectralCentroid; }, false);
        add_f32("spectral_centroid_hz", spectralCentroid, /*smooth=*/true, /*normalize01=*/false);
        auto hfRatio = extract_and_resample(ctx.audioFeatures, t, [](const AudioFeatures& f){ return f.spectralBalance; }, false);
        add_f32("spectral_hf_ratio", hfRatio, /*smooth=*/true, /*normalize01=*/false);
        add_f32("spectral_brightness", spectralCentroid, /*smooth=*/true, /*normalize01=*/false);

        auto dialogueClarity = extract_and_resample(ctx.audioFeatures, t, [](const AudioFeatures& f){ return f.dialogueClarity; }, false);
        add_f32("dialogue_clarity", dialogueClarity, /*smooth=*/true, /*normalize01=*/false);

        auto speechProb = extract_and_resample(ctx.audioFeatures, t, [](const AudioFeatures& f){ return f.speechProbability; }, false);
        add_f32("speech_probability", speechProb, /*smooth=*/true, /*normalize01=*/false);

        auto bpm = extract_and_resample(ctx.audioFeatures, t, [](const AudioFeatures& f){ return f.bpm; }, false);
        add_f32("bpm", bpm, /*smooth=*/true, /*normalize01=*/false);

        auto tempoConf = extract_and_resample(ctx.audioFeatures, t, [](const AudioFeatures& f){ return f.tempoConfidence; }, false);
        add_f32("tempo_confidence", tempoConf, /*smooth=*/true, /*normalize01=*/false);

        auto tension = extract_and_resample(ctx.audioFeatures, t, [](const AudioFeatures& f){ return f.harmonicTension; }, false);
        add_f32("harmonic_tension", tension, /*smooth=*/true, /*normalize01=*/false);

        auto stereoWidth = extract_and_resample(ctx.audioFeatures, t, [](const AudioFeatures& f){ return f.stereoWidth; }, false);
        add_f32("stereo_width", stereoWidth, /*smooth=*/true, /*normalize01=*/false);

        auto reverbAmount = extract_and_resample(ctx.audioFeatures, t, [](const AudioFeatures& f){ return f.reverbAmount; }, false);
        add_f32("reverb_amount", reverbAmount, /*smooth=*/true, /*normalize01=*/false);

        // Audio role energies (doc 4.2.6): raw per-second energies (power units).
        auto dEnergy = extract_and_resample(ctx.audioFeatures, t, [](const AudioFeatures& f){ return f.dialogueEnergy; }, false);
        add_f32("dialogue_energy", dEnergy, /*smooth=*/true, /*normalize01=*/false);
        auto mEnergy = extract_and_resample(ctx.audioFeatures, t, [](const AudioFeatures& f){ return f.musicEnergy; }, false);
        add_f32("music_energy", mEnergy, /*smooth=*/true, /*normalize01=*/false);
        auto eEnergy = extract_and_resample(ctx.audioFeatures, t, [](const AudioFeatures& f){ return f.effectsEnergy; }, false);
        add_f32("effects_energy", eEnergy, /*smooth=*/true, /*normalize01=*/false);

        // SilenceScore: invert robust-sigmoid loudness (no fixed threshold).
        {
            auto loud = extract_and_resample(ctx.audioFeatures, t, [](const AudioFeatures& f){ return f.loudness; }, false);
            auto loud01 = arcscope::robust_sigmoid01(loud);
            std::vector<double> silence(loud01.size(), 0.0);
            for (std::size_t i = 0; i < loud01.size(); ++i) silence[i] = 1.0 - loud01[i];
            add_f32("silence_score", silence, /*smooth=*/true, /*normalize01=*/false);
        }

        // Dialogue/Music/Effects dominance stripes (0/1 bands): argmax over the three dominance curves.
        {
            auto d = extract_and_resample(ctx.audioFeatures, t, [](const AudioFeatures& f){ return f.dialogueDominance; }, false);
            auto m = extract_and_resample(ctx.audioFeatures, t, [](const AudioFeatures& f){ return f.musicDominance; }, false);
            auto e = extract_and_resample(ctx.audioFeatures, t, [](const AudioFeatures& f){ return f.effectsDominance; }, false);
            std::vector<int32_t> d01(N, 0), m01(N, 0), e01(N, 0);
            for (std::size_t i = 0; i < N; ++i) {
                const double dv = (i < d.size()) ? d[i] : 0.0;
                const double mv = (i < m.size()) ? m[i] : 0.0;
                const double ev = (i < e.size()) ? e[i] : 0.0;
                if (dv >= mv && dv >= ev) d01[i] = 1;
                else if (mv >= dv && mv >= ev) m01[i] = 1;
                else e01[i] = 1;
            }
            add_i32("dialogue_dominance", d01);
            add_i32("music_dominance", m01);
            add_i32("effects_dominance", e01);
        }

        // IsSpeaking(t): 0/1 band derived from speechProbability with quantile threshold and segment cleanup.
        {
            auto sp = extract_and_resample(ctx.audioFeatures, t, [](const AudioFeatures& f){ return f.speechProbability; }, false);
            if (sp.size() == N) {
                const double thr = arcscope::quantile(sp, 0.75);
                std::vector<int32_t> speak(N, 0);
                for (std::size_t i = 0; i < N; ++i) speak[i] = (sp[i] > thr) ? 1 : 0;

                // Remove very short speaking segments using the film's own stable-length statistic (no fixed seconds).
                std::vector<double> mask01(N, 0.0);
                for (std::size_t i = 0; i < N; ++i) mask01[i] = speak[i] ? 1.0 : 0.0;
                const double Tmin = std::max(1.0, stable_duration_above_quantile(sp, 0.75, 1.0));
                const std::size_t minLen = std::max<std::size_t>(1, (std::size_t)std::llround(0.5 * Tmin));
                std::size_t i = 0;
                while (i < N) {
                    while (i < N && speak[i] == 0) ++i;
                    if (i >= N) break;
                    std::size_t j = i;
                    while (j < N && speak[j] == 1) ++j;
                    if (j - i < minLen) {
                        for (std::size_t k = i; k < j; ++k) speak[k] = 0;
                    }
                    i = j;
                }
                add_i32("is_speaking", speak);
            }
        }

        std::vector<int32_t> keyPc(N, -1);
        for (std::size_t i = 0; i < N; ++i) {
            if (i < ctx.audioFeatures.size()) keyPc[i] = (int32_t)ctx.audioFeatures[i].estimatedKey;
        }
        add_i32("key_pc", keyPc);

        // KeyChangeDensity (doc 5.6): change tick smoothed with adaptive window (from Pace λ).
        {
            std::vector<double> tick(N, 0.0);
            for (std::size_t i = 1; i < N; ++i) {
                const int32_t a = keyPc[i - 1];
                const int32_t b = keyPc[i];
                tick[i] = (a >= 0 && b >= 0 && a != b) ? 1.0 : 0.0;
            }

            double windowSec = 12.0;
            for (const auto& c : curves) {
                if (c.kind == CurveKind::Pace && c.raw_score.size() >= 2) {
                    double lambda = 0.0;
                    for (std::size_t k = 1; k < c.raw_score.size(); ++k) lambda += std::abs(c.raw_score[k] - c.raw_score[k - 1]);
                    lambda /= static_cast<double>(c.raw_score.size() - 1);
                    windowSec = contract::adaptive_window_size(lambda);
                    break;
                }
            }
            const std::size_t w = std::max<std::size_t>(1, (std::size_t)std::llround(windowSec));
            std::vector<double> density(N, 0.0);
            double acc = 0.0;
            for (std::size_t i = 0; i < N; ++i) {
                acc += tick[i];
                if (i >= w) acc -= tick[i - w];
                const std::size_t den = std::min<std::size_t>(i + 1, w);
                density[i] = (den > 0) ? (acc / (double)den) : 0.0;
            }
            add_f32("key_change_density", density, /*smooth=*/true, /*normalize01=*/false);
        }

        // 12-bin chroma distribution per second (doc 4.2.x; display-only).
        MetadataOverlaySeries chroma;
        chroma.type = "chroma12";
        chroma.fps = 1.0;
        chroma.channels = 12;
        chroma.valueType = OverlayValueType::Float32;
        chroma.valuesF32.reserve(N * 12);
        for (std::size_t i = 0; i < N; ++i) {
            if (i < ctx.audioFeatures.size() && ctx.audioFeatures[i].chroma.size() == 12) {
                for (int b = 0; b < 12; ++b) chroma.valuesF32.push_back((float)ctx.audioFeatures[i].chroma[(std::size_t)b]);
            } else {
                for (int b = 0; b < 12; ++b) chroma.valuesF32.push_back(0.0f);
            }
        }
        out.push_back(std::move(chroma));
    }

    // ---- Motion/camera overlays (doc 5.6 Photography) ----
    if (!ctx.t_center.empty()) {
        if (!ctx.camera_motion_1hz.empty()) {
            auto cam = arcscope::resample_linear(ctx.t_center, ctx.camera_motion_1hz, t);
            add_f32("camera_motion", cam, /*smooth=*/true, /*normalize01=*/false);
        }
        if (!ctx.object_motion_1hz.empty()) {
            auto obj = arcscope::resample_linear(ctx.t_center, ctx.object_motion_1hz, t);
            add_f32("object_motion", obj, /*smooth=*/true, /*normalize01=*/false);
        }
        if (!ctx.camera_motion_type_1hz.empty()) {
            std::vector<double> src(ctx.camera_motion_type_1hz.size(), 0.0);
            for (std::size_t i = 0; i < ctx.camera_motion_type_1hz.size(); ++i) src[i] = (double)ctx.camera_motion_type_1hz[i];
            auto v = arcscope::resample_nearest(ctx.t_center, src, t);
            std::vector<int32_t> outType(N, 0);
            for (std::size_t i = 0; i < N && i < v.size(); ++i) outType[i] = (int32_t)std::llround(v[i]);
            add_i32("camera_motion_type", outType);
        }
        if (!ctx.camera_motion_1hz.empty() && !ctx.object_motion_1hz.empty()) {
            const std::size_t n = std::min(ctx.camera_motion_1hz.size(), ctx.object_motion_1hz.size());
            std::vector<double> delta(n, 0.0);
            for (std::size_t i = 0; i < n; ++i) delta[i] = ctx.camera_motion_1hz[i] - ctx.object_motion_1hz[i];
            auto d = arcscope::resample_linear(ctx.t_center, delta, t);
            add_f32("camera_vs_object_delta", d, /*smooth=*/true, /*normalize01=*/false);
        }
    }

    // ---- Face/Dialogue overlays (doc 4.5 + 5.6/10.7) ----
    // Note: FacePresence is persisted as a curve (curve_type='face_presence'), not as an overlay.
    if (!ctx.faceFeatures.empty()) {
        auto presence = extract_and_resample(ctx.faceFeatures, t, [](const FaceFeatures& f){ return f.facePresence; }, true);

        auto yaw = extract_and_resample(ctx.faceFeatures, t, [](const FaceFeatures& f){ return f.yaw; }, false);
        auto roll = extract_and_resample(ctx.faceFeatures, t, [](const FaceFeatures& f){ return f.roll; }, false);
        auto downRatio = extract_and_resample(ctx.faceFeatures, t, [](const FaceFeatures& f){ return f.downcastRatio; }, false);

        // FrontalScore: high when |yaw| and |roll| are small (pose close to camera).
        // Use robust normalization within this film (no fixed degrees/radians thresholds).
        {
            std::vector<double> poseMag(N, 0.0);
            for (std::size_t i = 0; i < N; ++i) {
                const double y = (i < yaw.size()) ? yaw[i] : 0.0;
                const double r = (i < roll.size()) ? roll[i] : 0.0;
                poseMag[i] = std::sqrt(y * y + r * r);
            }
            auto inv = poseMag;
            for (auto& v : inv) v = -v;
            auto frontal = arcscope::robust_sigmoid01(inv);
            add_f32("frontal_score", frontal, /*smooth=*/true, /*normalize01=*/false);
        }

        // DowncastScore: high when downcastRatio is small (nose appears closer to eyes).
        {
            auto inv = downRatio;
            for (auto& v : inv) v = -v;
            auto down = arcscope::robust_sigmoid01(inv);
            add_f32("downcast_score", down, /*smooth=*/true, /*normalize01=*/false);
        }

        // LookAtCamera: 0/1 band based on combined score (Frontal high and Downcast low), gated by FacePresence.
        {
            std::vector<double> poseMag(N, 0.0);
            for (std::size_t i = 0; i < N; ++i) {
                const double y = (i < yaw.size()) ? yaw[i] : 0.0;
                const double r = (i < roll.size()) ? roll[i] : 0.0;
                poseMag[i] = std::sqrt(y * y + r * r);
            }
            auto invPose = poseMag;
            for (auto& v : invPose) v = -v;
            const auto frontal = arcscope::robust_sigmoid01(invPose);

            auto invDown = downRatio;
            for (auto& v : invDown) v = -v;
            const auto down = arcscope::robust_sigmoid01(invDown);

            std::vector<double> score(N, 0.0);
            std::vector<double> presentScores;
            presentScores.reserve(N);
            for (std::size_t i = 0; i < N; ++i) {
                const double p = (i < presence.size()) ? presence[i] : 0.0;
                score[i] = (p > 0.0) ? (frontal[i] * (1.0 - down[i])) : 0.0;
                if (p > 0.0) presentScores.push_back(score[i]);
            }
            const double thr = presentScores.empty() ? 1.0 : arcscope::quantile(presentScores, 0.85);
            std::vector<int32_t> look(N, 0);
            for (std::size_t i = 0; i < N; ++i) {
                const double p = (i < presence.size()) ? presence[i] : 0.0;
                if (p > 0.0 && score[i] > thr) look[i] = 1;
            }
            add_i32("look_at_camera", look);
        }
    }

    // ---- Info overlays (doc 4.4 + 5.6/10.7) ----
    if (!ctx.infoFeatures.empty()) {
        auto expo = extract_and_resample(ctx.infoFeatures, t, [](const InfoFeatures& f){ return f.exposition; }, false);
        add_f32("exposition_curve", expo, /*smooth=*/true, /*normalize01=*/false);

        auto verbal = extract_and_resample(ctx.infoFeatures, t, [](const InfoFeatures& f){ return f.verbalLoad; }, false);
        auto visual = extract_and_resample(ctx.infoFeatures, t, [](const InfoFeatures& f){ return f.visualLoad; }, false);
        auto eventl = extract_and_resample(ctx.infoFeatures, t, [](const InfoFeatures& f){ return f.eventLoad; }, false);
        add_f32("verbal_load", verbal, /*smooth=*/true, /*normalize01=*/true);
        add_f32("visual_load", visual, /*smooth=*/true, /*normalize01=*/true);
        add_f32("event_load", eventl, /*smooth=*/true, /*normalize01=*/true);

        // VisualLoad subfeatures (doc 4.4.3): VisualEntropy / ShotScaleNumeric / CameraMotionComplexity / FaceChange.
        auto visualEntropy = extract_and_resample(ctx.infoFeatures, t, [](const InfoFeatures& f){ return f.visualEntropy; }, false);
        auto shotScale = extract_and_resample(ctx.infoFeatures, t, [](const InfoFeatures& f){ return f.shotScaleNumeric; }, false);
        auto camMotion = extract_and_resample(ctx.infoFeatures, t, [](const InfoFeatures& f){ return f.cameraMotionComplexity; }, false);
        auto faceChange = extract_and_resample(ctx.infoFeatures, t, [](const InfoFeatures& f){ return f.faceChange; }, false);
        add_f32("visual_entropy", visualEntropy, /*smooth=*/true, /*normalize01=*/true);
        add_f32("shot_scale_numeric", shotScale, /*smooth=*/true, /*normalize01=*/true);
        add_f32("camera_motion_complexity", camMotion, /*smooth=*/true, /*normalize01=*/true);
        add_f32("face_change", faceChange, /*smooth=*/true, /*normalize01=*/true);
    }

    // ---- Photography overlays: ShotScale band (doc 5.6) ----
    if (!ctx.faceFeatures.empty()) {
        auto presence = extract_and_resample(ctx.faceFeatures, t, [](const FaceFeatures& f){ return f.facePresence; }, true);
        auto cuw = extract_and_resample(ctx.faceFeatures, t, [](const FaceFeatures& f){ return f.closeUpWeight; }, true);

        std::vector<double> present;
        present.reserve(N);
        for (std::size_t i = 0; i < N; ++i) {
            if (presence[i] > 0.0) present.push_back(std::clamp(cuw[i], 0.0, 1.0));
        }

        // 6 bins for ELS→ECU using per-film quantiles (adaptive, no fixed thresholds).
        // If there is no dominant-face presence at all, the band stays all zeros.
        std::array<double, 5> thr{};
        if (!present.empty()) {
            thr[0] = arcscope::quantile(present, 1.0 / 6.0);
            thr[1] = arcscope::quantile(present, 2.0 / 6.0);
            thr[2] = arcscope::quantile(present, 3.0 / 6.0);
            thr[3] = arcscope::quantile(present, 4.0 / 6.0);
            thr[4] = arcscope::quantile(present, 5.0 / 6.0);
        }

        std::vector<int32_t> shotScale(N, 0);
        for (std::size_t i = 0; i < N; ++i) {
            if (!(presence[i] > 0.0) || present.empty()) { shotScale[i] = 0; continue; } // 0=no dominant face / off-screen
            const double v = std::clamp(cuw[i], 0.0, 1.0);
            int bin = 0;
            while (bin < 5 && v > thr[(std::size_t)bin]) ++bin;
            shotScale[i] = 1 + bin; // 1..6
        }
        add_i32("shot_scale", shotScale);
    }

    return out;
}

}  // namespace arcscope
