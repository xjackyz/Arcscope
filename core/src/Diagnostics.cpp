#include "arcscope/Diagnostics.h"
#include "arcscope/CoreContract.h"
#include "arcscope/Normalization.h"

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <limits>
#include <optional>
#include <string_view>
#include <unordered_map>
#include <vector>
#include <functional>

namespace {

using arcscope::CurveKind;
using arcscope::CurveSeries;
using arcscope::IssueSegment;
using arcscope::MetadataOverlaySeries;
using arcscope::OverlayValueType;

constexpr double kEps = 1e-9;

inline double clamp01(double x) { return std::min(1.0, std::max(0.0, x)); }

inline double safe_div(double num, double den) { return (std::abs(den) > kEps) ? (num / den) : 0.0; }

const CurveSeries* find_curve(const std::vector<CurveSeries>& curves, CurveKind kind) {
    for (const auto& c : curves) {
        if (c.kind == kind) return &c;
    }
    return nullptr;
}

const MetadataOverlaySeries* find_overlay(const std::vector<MetadataOverlaySeries>& overlays, std::string_view type) {
    for (const auto& o : overlays) {
        if (o.type == type) return &o;
    }
    return nullptr;
}

std::optional<std::vector<double>> overlay_f32_1ch(const MetadataOverlaySeries* o, std::size_t N) {
    if (!o) return std::nullopt;
    if (o->valueType != OverlayValueType::Float32 || o->channels != 1) return std::nullopt;
    if (o->valuesF32.size() != N) return std::nullopt;
    std::vector<double> out(N, 0.0);
    for (std::size_t i = 0; i < N; ++i) out[i] = static_cast<double>(o->valuesF32[i]);
    return out;
}

double q(const std::vector<double>& v, double p) {
    if (v.empty()) return 0.0;
    return arcscope::quantile(std::vector<double>(v.begin(), v.end()), p);
}

double estimate_dt_seconds(const std::vector<double>& times) {
    if (times.size() < 2) return 1.0;
    std::vector<double> diffs;
    diffs.reserve(times.size() - 1);
    for (std::size_t i = 1; i < times.size(); ++i) diffs.push_back(times[i] - times[i - 1]);
    double dt = arcscope::median(std::move(diffs));
    if (!std::isfinite(dt) || dt <= 1e-6) dt = 1.0;
    return dt;
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

// T_min^(X): median length of segments where X(t) > Q_X(qLevel) (doc 6.2).
double stable_duration_above_quantile(const std::vector<double>& x, double qLevel, double dtSeconds) {
    if (x.empty()) return 0.0;
    const double thr = arcscope::quantile(std::vector<double>(x.begin(), x.end()), qLevel);
    std::vector<double> mask01(x.size(), 0.0);
    for (std::size_t i = 0; i < x.size(); ++i) mask01[i] = (x[i] > thr) ? 1.0 : 0.0;
    return median_stable_duration_from_mask01(mask01, dtSeconds);
}

std::unordered_map<CurveKind, std::vector<double>> calculate_z_scores(const std::vector<CurveSeries>& curves) {
    std::unordered_map<CurveKind, std::vector<double>> z;
    for (const auto& c : curves) {
        if (c.values.empty()) continue;
        auto zv = arcscope::robust_z(c.values);
        if (zv.size() != c.values.size()) zv.resize(c.values.size(), 0.0);
        z[c.kind] = std::move(zv);
    }
    return z;
}

struct SegmentIdx {
    std::size_t i{0}; // inclusive
    std::size_t j{0}; // exclusive
    double meanSeverity{0.0};
    double startTime{0.0}; // inclusive
    double endTime{0.0};   // exclusive
};

std::vector<SegmentIdx> scan_segments(const std::vector<double>& times,
                                      const std::vector<double>& severity01,
                                      double minDurationSec) {
    std::vector<SegmentIdx> out;
    if (times.empty() || severity01.empty() || times.size() != severity01.size()) return out;
    const std::size_t N = severity01.size();
    const double dt = estimate_dt_seconds(times);
    const double minDur = std::max(minDurationSec, dt);

    std::size_t i = 0;
    while (i < N) {
        while (i < N && severity01[i] <= 0.0) ++i;
        if (i >= N) break;
        std::size_t j = i;
        while (j < N && severity01[j] > 0.0) ++j;

        const double startT = std::max(0.0, times[i] - 0.5 * dt);
        const double endT = times[j - 1] + 0.5 * dt;
        const double dur = static_cast<double>(j - i) * dt;
        if (dur >= minDur) {
            double s = 0.0;
            for (std::size_t k = i; k < j; ++k) s += severity01[k];
            s /= static_cast<double>(j - i);
            out.push_back(SegmentIdx{.i = i, .j = j, .meanSeverity = clamp01(s), .startTime = startT, .endTime = endT});
        }
        i = j;
    }
    return out;
}

std::vector<IssueSegment> emit_segments(const std::vector<double>& times,
                                        const std::vector<double>& severity01,
                                        const std::string& type,
                                        const std::string& explanation,
                                        double minDurationSec) {
    std::vector<IssueSegment> out;
    for (const auto& seg : scan_segments(times, severity01, minDurationSec)) {
        IssueSegment s;
        s.type = type;
        s.startTime = seg.startTime;
        s.endTime = seg.endTime;
        s.severity = static_cast<float>(seg.meanSeverity);
        s.explanation = explanation;
        out.push_back(std::move(s));
    }
    return out;
}

} // namespace

namespace arcscope {

std::vector<IssueSegment> DiagnosticsEngine::detect(const std::vector<CurveSeries>& curves,
                                                    const std::vector<MetadataOverlaySeries>& overlays) const {
    std::vector<IssueSegment> issues;
    if (curves.empty()) return issues;

    const CurveSeries* pace = find_curve(curves, CurveKind::Pace);
    const CurveSeries* sound = find_curve(curves, CurveKind::Sound);
    const CurveSeries* color = find_curve(curves, CurveKind::Color);
    const CurveSeries* info = find_curve(curves, CurveKind::Info);
    const CurveSeries* arousal = find_curve(curves, CurveKind::Arousal);
    const CurveSeries* face = find_curve(curves, CurveKind::FaceAffect);

    if (!pace || !sound || !arousal) return issues;
    if (pace->times.size() != pace->values.size()) return issues;
    if (pace->values.empty() || sound->values.empty() || arousal->values.empty()) return issues;

    const std::size_t N = pace->values.size();
    if (sound->values.size() != N || arousal->values.size() != N) return issues;
    if (info && info->values.size() != N) info = nullptr;
    if (color && color->values.size() != N) color = nullptr;
    if (face && face->values.size() != N) face = nullptr;

    const auto& times = pace->times;
    const double dt = estimate_dt_seconds(times);

    const auto zscores = calculate_z_scores(curves);
    if (!zscores.count(CurveKind::Pace) || !zscores.count(CurveKind::Sound) || !zscores.count(CurveKind::Arousal)) return issues;
    const auto& zP = zscores.at(CurveKind::Pace);
    const auto& zS = zscores.at(CurveKind::Sound);
    const auto& zA = zscores.at(CurveKind::Arousal);

    const std::vector<double> zI = (info && zscores.count(CurveKind::Info)) ? zscores.at(CurveKind::Info) : std::vector<double>(N, 0.0);
    const std::vector<double> zF = (face && zscores.count(CurveKind::FaceAffect)) ? zscores.at(CurveKind::FaceAffect) : std::vector<double>(N, 0.0);
    const std::vector<double> zC = (color && zscores.count(CurveKind::Color)) ? zscores.at(CurveKind::Color) : std::vector<double>(N, 0.0);

    const double QP25 = q(pace->values, 0.25);
    const double QP75 = q(pace->values, 0.75);
    const double QP85 = q(pace->values, 0.85);
    const double QP95 = q(pace->values, 0.95);

    const double QS25 = q(sound->values, 0.25);
    const double QS30 = q(sound->values, 0.30);
    const double QS75 = q(sound->values, 0.75);

    const double QA25 = q(arousal->values, 0.25);
    const double QA30 = q(arousal->values, 0.30);
    const double QA70 = q(arousal->values, 0.70);
    const double QA85 = q(arousal->values, 0.85);
    const double QA95 = q(arousal->values, 0.95);

    const double TminP = std::max(dt, stable_duration_above_quantile(pace->values, 0.5, dt));
    const double TminS = std::max(dt, stable_duration_above_quantile(sound->values, 0.5, dt));
    const double TminA = std::max(dt, stable_duration_above_quantile(arousal->values, 0.5, dt));

    // 6.3.1 LowActivity
    if (info && !info->values.empty()) {
        const double QI25 = q(info->values, 0.25);
        std::vector<double> sev(N, 0.0);
        const double dp = std::max(QP25, kEps);
        const double ds = std::max(QS25, kEps);
        const double da = std::max(QA25, kEps);
        const double di = std::max(QI25, kEps);
        for (std::size_t k = 0; k < N; ++k) {
            if (pace->values[k] < QP25 && sound->values[k] < QS25 && arousal->values[k] < QA25 && info->values[k] < QI25) {
                const double sp = safe_div(QP25 - pace->values[k], dp);
                const double ss = safe_div(QS25 - sound->values[k], ds);
                const double sa = safe_div(QA25 - arousal->values[k], da);
                const double si = safe_div(QI25 - info->values[k], di);
                sev[k] = clamp01((sp + ss + sa + si) / 4.0);
            }
        }
        auto segs = emit_segments(times, sev, "LowActivity",
                                  "LowActivity: P,S,A,I all < Q(0.25) (doc 6.3.1).",
                                  TminP);
        issues.insert(issues.end(), segs.begin(), segs.end());
    }

    // 6.3.2 OvercutFlat
    {
        std::vector<double> sev(N, 0.0);
        const double denomP = std::max(QP95 - QP75, kEps);
        const double denomA = std::max(QA25 - q(arousal->values, 0.05), kEps);
        for (std::size_t k = 0; k < N; ++k) {
            if (pace->values[k] > QP75 && arousal->values[k] < QA25) {
                const double sp = clamp01(safe_div(pace->values[k] - QP75, denomP));
                const double sa = clamp01(safe_div(QA25 - arousal->values[k], denomA));
                sev[k] = clamp01(0.5 * sp + 0.5 * sa);
            }
        }
        auto segs = emit_segments(times, sev, "OvercutFlat",
                                  "OvercutFlat: P>Q_P(0.75) while A<Q_A(0.25) (doc 6.3.2).",
                                  0.5 * TminP);
        issues.insert(issues.end(), segs.begin(), segs.end());
    }

    // 6.3.3 AudioVisualMisalign (global AV lag search ±3s, then |Z_P - Z_S| outliers)
    {
        int bestLag = 0;
        double bestScore = -std::numeric_limits<double>::infinity();
        for (int lag = -3; lag <= 3; ++lag) {
            double score = 0.0;
            for (std::size_t k = 0; k < N; ++k) {
                const long kk = static_cast<long>(k) + static_cast<long>(lag);
                if (kk < 0 || kk >= static_cast<long>(N)) continue;
                score += zP[k] * zS[static_cast<std::size_t>(kk)];
            }
            if (score > bestScore + 1e-12) {
                bestScore = score;
                bestLag = lag;
            }
        }

        std::vector<double> absd;
        absd.reserve(N);
        std::vector<double> dps(N, 0.0);
        std::vector<double> valid01(N, 0.0);
        for (std::size_t k = 0; k < N; ++k) {
            const long kk = static_cast<long>(k) + static_cast<long>(bestLag);
            if (kk < 0 || kk >= static_cast<long>(N)) continue;
            dps[k] = zP[k] - zS[static_cast<std::size_t>(kk)];
            valid01[k] = 1.0;
            absd.push_back(std::abs(dps[k]));
        }

        if (!absd.empty()) {
            const double thr = arcscope::quantile(absd, 0.8);
            const double hi = std::max(arcscope::quantile(absd, 0.95), thr + kEps);
            const double den = std::max(hi - thr, kEps);
            std::vector<double> sev(N, 0.0);
            for (std::size_t k = 0; k < N; ++k) {
                if (valid01[k] > 0.0) {
                    const double a = std::abs(dps[k]);
                    if (a > thr) sev[k] = clamp01(safe_div(a - thr, den));
                }
            }
            for (const auto& seg : scan_segments(times, sev, TminS)) {
                double meanD = 0.0;
                double meanAbsD = 0.0;
                for (std::size_t k = seg.i; k < seg.j; ++k) {
                    meanD += dps[k];
                    meanAbsD += std::abs(dps[k]);
                }
                const double denomSeg = std::max<double>(1.0, static_cast<double>(seg.j - seg.i));
                meanD /= denomSeg;
                meanAbsD /= denomSeg;

                IssueSegment s;
                s.type = "AudioVisualMisalign";
                s.startTime = seg.startTime;
                s.endTime = seg.endTime;
                s.severity = static_cast<float>(seg.meanSeverity);
                const std::string direction = (meanD >= 0.0)
                    ? "D_PS>0 (visuals fast / sound weak)"
                    : "D_PS<0 (sound strong / visuals flat)";
                s.explanation =
                    "AudioVisualMisalign: AV lag ℓ*=" + std::to_string(bestLag) +
                    "s, |D_PS|=|Z_P-Z_S(shifted)|>Q(0.8) (doc 6.3.3); " + direction + ".";
                issues.push_back(std::move(s));
            }
        }
    }

    // 6.3.4 InfoKillsEmotion (needs exposition_curve overlay)
    if (info && !info->values.empty()) {
        const auto expoOpt = overlay_f32_1ch(find_overlay(overlays, "exposition_curve"), N);
        if (expoOpt) {
            const auto& expo = *expoOpt;
            const double QI80 = q(info->values, 0.8);
            const double QExpo60 = q(expo, 0.6);
            std::vector<double> sev(N, 0.0);
            for (std::size_t k = 0; k < N; ++k) {
                if (info->values[k] > QI80 && arousal->values[k] < QA30 && expo[k] > QExpo60) {
                    const double si = clamp01(safe_div(info->values[k] - QI80, std::max(q(info->values, 0.95) - QI80, kEps)));
                    const double sa = clamp01(safe_div(QA30 - arousal->values[k], std::max(QA30 - q(arousal->values, 0.05), kEps)));
                    sev[k] = clamp01(si * sa);
                }
            }
            const double TminIHigh = std::max(dt, stable_duration_above_quantile(info->values, 0.8, dt));
            auto segs = emit_segments(times, sev, "InfoKillsEmotion",
                                      "InfoKillsEmotion: I>Q_I(0.8) and A<Q_A(0.3) while Exposition>Q_Expo(0.6) (doc 6.3.4).",
                                      TminIHigh);
            issues.insert(issues.end(), segs.begin(), segs.end());
        }
    }

    // 6.3.5 HighTensionLowPace
    {
        std::vector<double> raw(N, 0.0);
        std::vector<double> rawCond;
        rawCond.reserve(N);

        const double denomA = std::max(QA85, kEps);
        for (std::size_t k = 0; k < N; ++k) {
            if (arousal->values[k] > QA85 && pace->values[k] < QP25) {
                // doc 6.3.5: severity = mean((A - Q_A(0.85)) / Q_A(0.85)).
                raw[k] = safe_div(arousal->values[k] - QA85, denomA);
                rawCond.push_back(raw[k]);
            }
        }

        std::vector<double> sev(N, 0.0);
        if (!rawCond.empty()) {
            const double scale = std::max(arcscope::quantile(rawCond, 0.95), kEps);
            for (std::size_t k = 0; k < N; ++k) {
                if (raw[k] > 0.0) sev[k] = clamp01(safe_div(raw[k], scale));
            }
        }
        auto segs = emit_segments(times, sev, "HighTensionLowPace",
                                  "HighTensionLowPace: A>Q_A(0.85) while P<Q_P(0.25) (doc 6.3.5).",
                                  0.6 * TminP);
        issues.insert(issues.end(), segs.begin(), segs.end());
    }

    // 6.3.6 SoundOverdriveFlatEmotion
    {
        // doc 6.3.6: severity = mean(Z_S(t) + Z_LowA(t)).
        // We define LowA(t) = (Q_A(0.25) - A(t))_+ and Z_LowA via robust z-score.
        std::vector<double> lowA(N, 0.0);
        for (std::size_t k = 0; k < N; ++k) lowA[k] = std::max(0.0, QA25 - arousal->values[k]);
        const auto zLowA = arcscope::robust_z(lowA);

        std::vector<double> raw(N, 0.0);
        std::vector<double> rawCond;
        rawCond.reserve(N);
        for (std::size_t k = 0; k < N; ++k) {
            if (sound->values[k] > QS75 && arousal->values[k] < QA25) {
                raw[k] = zS[k] + zLowA[k];
                rawCond.push_back(raw[k]);
            }
        }
        if (!rawCond.empty()) {
            const double scale = std::max(arcscope::quantile(rawCond, 0.95), kEps);
            std::vector<double> sev(N, 0.0);
            for (std::size_t k = 0; k < N; ++k) {
                if (sound->values[k] > QS75 && arousal->values[k] < QA25) sev[k] = clamp01(safe_div(raw[k], scale));
            }
            auto segs = emit_segments(times, sev, "SoundOverdriveFlatEmotion",
                                      "SoundOverdriveFlatEmotion: S>Q_S(0.75) while A<Q_A(0.25) (doc 6.3.6).",
                                      TminS);
            issues.insert(issues.end(), segs.begin(), segs.end());
        }
    }

    // 6.3.7 SoundColorConflict (ColdDark = Z_gray + Z_coolHue)
    {
        const auto warmthOpt = overlay_f32_1ch(find_overlay(overlays, "color_warmth"), N);
        const auto grayOpt = overlay_f32_1ch(find_overlay(overlays, "color_grayness"), N);
        if (warmthOpt && grayOpt) {
            const auto& warmth = *warmthOpt;
            const auto& gray = *grayOpt;

            std::vector<double> cool(N, 0.0);
            for (std::size_t k = 0; k < N; ++k) cool[k] = 1.0 - warmth[k];
            const auto zGray = arcscope::robust_z(gray);
            const auto zCool = arcscope::robust_z(cool);
            std::vector<double> coldDark(N, 0.0);
            for (std::size_t k = 0; k < N; ++k) coldDark[k] = zGray[k] + zCool[k];

            const double QCold75 = q(coldDark, 0.75);
            const double TminCold = std::max(dt, stable_duration_above_quantile(coldDark, 0.5, dt));
            const double minDur = std::min(TminS, TminCold);

            // doc 6.3.7: severity = mean(Z_S(t) + Z_Cold(t)), where Z_Cold is ColdDark in z-space.
            std::vector<double> raw(N, 0.0);
            std::vector<double> rawCond;
            rawCond.reserve(N);
            for (std::size_t k = 0; k < N; ++k) {
                if (sound->values[k] > QS75 && coldDark[k] > QCold75) {
                    raw[k] = zS[k] + coldDark[k];
                    rawCond.push_back(raw[k]);
                }
            }
            if (!rawCond.empty()) {
                const double scale = std::max(arcscope::quantile(rawCond, 0.95), kEps);
                std::vector<double> sev(N, 0.0);
                for (std::size_t k = 0; k < N; ++k) {
                    if (sound->values[k] > QS75 && coldDark[k] > QCold75) sev[k] = clamp01(safe_div(raw[k], scale));
                }
                auto segs = emit_segments(times, sev, "SoundColorConflict",
                                          "SoundColorConflict: S>Q_S(0.75) and ColdDark=Z_Gray+Z_CoolHue > Q_Cold(0.75) (doc 6.3.7).",
                                          minDur);
                issues.insert(issues.end(), segs.begin(), segs.end());
            }
        }
    }

    // 6.3.8 ColorEmotionMismatch (Warm/Cold), requires AlignNorm (4.3.4) + brightness + ColdDark
    {
        const auto alignOpt = overlay_f32_1ch(find_overlay(overlays, "align_norm"), N);
        const auto warmthOpt = overlay_f32_1ch(find_overlay(overlays, "color_warmth"), N);
        const auto brightOpt = overlay_f32_1ch(find_overlay(overlays, "color_brightness"), N);
        const auto grayOpt = overlay_f32_1ch(find_overlay(overlays, "color_grayness"), N);

        if (alignOpt && warmthOpt && brightOpt && grayOpt) {
            const auto& align = *alignOpt;
            const auto& warmth = *warmthOpt;
            const auto& bright = *brightOpt;
            const auto& gray = *grayOpt;

            std::vector<double> cool(N, 0.0);
            for (std::size_t k = 0; k < N; ++k) cool[k] = 1.0 - warmth[k];
            const auto zGray = arcscope::robust_z(gray);
            const auto zCool = arcscope::robust_z(cool);
            std::vector<double> coldDark(N, 0.0);
            for (std::size_t k = 0; k < N; ++k) coldDark[k] = zGray[k] + zCool[k];

            const double QAlign80 = q(align, 0.8);
            const double QBright30 = q(bright, 0.3);
            const double QCold70 = q(coldDark, 0.7);

            const double QWarm50 = q(warmth, 0.5);
            const double QBright50 = q(bright, 0.5);

            std::vector<double> activeSev(N, 0.0);
            for (std::size_t k = 0; k < N; ++k) {
                const bool condAlign = (align[k] > QAlign80);
                const bool condBrightMismatch = (arousal->values[k] > QA70 && bright[k] < QBright30);
                const bool condColdMismatch = (arousal->values[k] < QA30 && coldDark[k] > QCold70);
                if (condAlign || condBrightMismatch || condColdMismatch) activeSev[k] = clamp01(align[k]);
            }

            for (const auto& seg : scan_segments(times, activeSev, TminA)) {
                double meanWarm = 0.0;
                double meanBright = 0.0;
                for (std::size_t k = seg.i; k < seg.j; ++k) {
                    meanWarm += warmth[k];
                    meanBright += bright[k];
                }
                const double denom = std::max<double>(1.0, static_cast<double>(seg.j - seg.i));
                meanWarm /= denom;
                meanBright /= denom;

                IssueSegment s;
                s.type = (meanWarm > QWarm50 && meanBright > QBright50) ? "ColorEmotionMismatchWarm" : "ColorEmotionMismatchCold";
                s.startTime = seg.startTime;
                s.endTime = seg.endTime;
                s.severity = static_cast<float>(seg.meanSeverity);
                s.explanation = "ColorEmotionMismatch: AlignNorm>Q_Align(0.8) (or proxy mismatch) sustained ≥Tmin(A) (doc 6.3.8).";
                issues.push_back(std::move(s));
            }
        }
    }

    // 6.3.9 CognitiveOverload
    if (info && !info->values.empty()) {
        const double QI75 = q(info->values, 0.75);
        std::vector<double> raw(N, 0.0);
        std::vector<double> rawCond;
        rawCond.reserve(N);
        for (std::size_t k = 0; k < N; ++k) {
            raw[k] = zI[k] + zP[k];
            if (info->values[k] > QI75 && pace->values[k] > QP75) rawCond.push_back(raw[k]);
        }
        if (!rawCond.empty()) {
            const double base = arcscope::quantile(rawCond, 0.5);
            const double hi = std::max(arcscope::quantile(rawCond, 0.95), base + kEps);
            const double den = std::max(hi - base, kEps);
            std::vector<double> sev(N, 0.0);
            for (std::size_t k = 0; k < N; ++k) {
                if (info->values[k] > QI75 && pace->values[k] > QP75) sev[k] = clamp01(safe_div(raw[k] - base, den));
            }
            const double TminI = std::max(dt, stable_duration_above_quantile(info->values, 0.5, dt));
            auto segs = emit_segments(times, sev, "CognitiveOverload",
                                      "CognitiveOverload: I>Q_I(0.75) and P>Q_P(0.75) (doc 6.3.9).",
                                      std::min(TminI, TminP));
            issues.insert(issues.end(), segs.begin(), segs.end());
        }
    }

    // 6.3.10 FacialAffectSuppressed
    if (face && !face->values.empty()) {
        const double QF75 = q(face->values, 0.75);
        std::vector<double> raw(N, 0.0);
        std::vector<double> rawCond;
        rawCond.reserve(N);
        for (std::size_t k = 0; k < N; ++k) {
            raw[k] = zF[k] - zA[k];
            if (face->values[k] > QF75 && arousal->values[k] < QA25) rawCond.push_back(raw[k]);
        }
        if (!rawCond.empty()) {
            const double base = arcscope::quantile(rawCond, 0.5);
            const double hi = std::max(arcscope::quantile(rawCond, 0.95), base + kEps);
            const double den = std::max(hi - base, kEps);
            std::vector<double> sev(N, 0.0);
            for (std::size_t k = 0; k < N; ++k) {
                if (face->values[k] > QF75 && arousal->values[k] < QA25) sev[k] = clamp01(safe_div(raw[k] - base, den));
            }
            const double TminF = std::max(dt, stable_duration_above_quantile(face->values, 0.5, dt));
            auto segs = emit_segments(times, sev, "FacialAffectSuppressed",
                                      "FacialAffectSuppressed: F>Q_F(0.75) while A<Q_A(0.25) (doc 6.3.10).",
                                      TminF);
            issues.insert(issues.end(), segs.begin(), segs.end());
        }
    }

    // 6.3.11 NonCharacterDrivenPeak
    if (face && !face->values.empty()) {
        const double QF25 = q(face->values, 0.25);
        std::vector<double> raw(N, 0.0);
        std::vector<double> rawCond;
        rawCond.reserve(N);
        for (std::size_t k = 0; k < N; ++k) {
            raw[k] = zA[k] - zF[k];
            if (arousal->values[k] > QA85 && face->values[k] < QF25) rawCond.push_back(raw[k]);
        }
        if (!rawCond.empty()) {
            const double base = arcscope::quantile(rawCond, 0.5);
            const double hi = std::max(arcscope::quantile(rawCond, 0.95), base + kEps);
            const double den = std::max(hi - base, kEps);
            std::vector<double> sev(N, 0.0);
            for (std::size_t k = 0; k < N; ++k) {
                if (arousal->values[k] > QA85 && face->values[k] < QF25) sev[k] = clamp01(safe_div(raw[k] - base, den));
            }
            auto segs = emit_segments(times, sev, "NonCharacterDrivenPeak",
                                      "NonCharacterDrivenPeak: A>Q_A(0.85) while F<Q_F(0.25) (doc 6.3.11).",
                                      0.5 * TminA);
            issues.insert(issues.end(), segs.begin(), segs.end());
        }
    }

    // 6.3.12 SilentHighTension (uses harmonic_tension overlay as Tension(t))
    {
        const auto tensionOpt = overlay_f32_1ch(find_overlay(overlays, "harmonic_tension"), N);
        if (tensionOpt) {
            const auto& tension = *tensionOpt;
            const double QT70 = q(tension, 0.7);
            const auto zT = arcscope::robust_z(tension);

            std::vector<double> raw(N, 0.0);
            std::vector<double> rawCond;
            rawCond.reserve(N);
            for (std::size_t k = 0; k < N; ++k) {
                if (sound->values[k] < QS30 && tension[k] > QT70) {
                    // doc 6.3.12: severity = mean(Z_Tension(t) - Z_S(t)).
                    raw[k] = zT[k] - zS[k];
                    rawCond.push_back(raw[k]);
                }
            }
            if (!rawCond.empty()) {
                const double scale = std::max(arcscope::quantile(rawCond, 0.95), kEps);
                std::vector<double> sev(N, 0.0);
                for (std::size_t k = 0; k < N; ++k) {
                    if (sound->values[k] < QS30 && tension[k] > QT70) sev[k] = clamp01(safe_div(raw[k], scale));
                }
                const double TminT = std::max(dt, stable_duration_above_quantile(tension, 0.5, dt));
                auto segs = emit_segments(times, sev, "SilentHighTension",
                                          "SilentHighTension: S<Q_S(0.3) while Tension>Q_T(0.7) (doc 6.3.12).",
                                          TminT);
                issues.insert(issues.end(), segs.begin(), segs.end());
            }
        }
    }

    // 6.3.13 ModeCounterpoint (window W from Pace, doc 3.3.1)
    {
        double lambda = 0.0;
        if (zP.size() >= 2) {
            double acc = 0.0;
            for (std::size_t k = 1; k < N; ++k) acc += std::abs(zP[k] - zP[k - 1]);
            lambda = acc / static_cast<double>(N - 1);
        }
        const double Wsec = arcscope::contract::adaptive_window_size(lambda);
        const std::size_t W = static_cast<std::size_t>(std::clamp(std::llround(Wsec / dt), 1LL, (long long)N));
        const double minDur = 0.5 * Wsec;

        auto sign = [](double x) -> int {
            if (x > 1e-6) return 1;
            if (x < -1e-6) return -1;
            return 0;
        };

        std::vector<double> conflict01(N, 0.0);
        for (std::size_t k = 0; k < N; ++k) {
            const int sA = sign(zA[k]);
            if (sA == 0) continue;
            int conflicts = 0;
            if (sign(zP[k]) == -sA) conflicts++;
            if (sign(zS[k]) == -sA) conflicts++;
            if (color && sign(zC[k]) == -sA) conflicts++;
            if (conflicts >= 2) conflict01[k] = 1.0;
        }

        std::vector<double> rolling(N, 0.0);
        if (W > 0) {
            double acc = 0.0;
            for (std::size_t k = 0; k < N; ++k) {
                acc += conflict01[k];
                if (k >= W) acc -= conflict01[k - W];
                const std::size_t denom = std::min<std::size_t>(k + 1, W);
                rolling[k] = (denom > 0) ? clamp01(acc / static_cast<double>(denom)) : 0.0;
            }
        }

        std::vector<double> sev(N, 0.0);
        for (std::size_t k = 0; k < N; ++k) sev[k] = (conflict01[k] > 0.0) ? rolling[k] : 0.0;
        auto segs = emit_segments(times, sev, "ModeCounterpoint",
                                  "ModeCounterpoint: ≥2 of {P,S,C} have sign opposite to A for ≥0.5W (doc 6.3.13).",
                                  minDur);
        issues.insert(issues.end(), segs.begin(), segs.end());
    }

    return issues;
}

} // namespace arcscope
