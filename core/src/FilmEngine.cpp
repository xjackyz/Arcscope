#include "arcscope/FilmEngine.h"
#include "arcscope/Utility.h"
#include "arcscope/Normalization.h"
#include "arcscope/CoreContract.h"

#include <algorithm>
#include <cmath>
#include <cstdint>
#include <numeric>

namespace arcscope {

namespace {

const CurveSeries* findCurve(const std::vector<CurveSeries>& curves, CurveKind kind) {
    for (const auto& c : curves) {
        if (c.kind == kind) return &c;
    }
    return nullptr;
}

std::vector<double> emotion_progress_boundaries_from_arousal(const std::vector<double>& arousal01) {
    const std::size_t N = arousal01.size();
    std::vector<double> u(N + 1, 0.0);
    double total = 0.0;
    for (double a : arousal01) total += std::max(0.0, a);
    if (!std::isfinite(total) || total <= 1e-12) {
        for (std::size_t k = 0; k <= N; ++k) u[k] = (N > 0) ? (static_cast<double>(k) / static_cast<double>(N)) : 0.0;
        return u;
    }
    double acc = 0.0;
    for (std::size_t k = 0; k < N; ++k) {
        acc += std::max(0.0, arousal01[k]);
        u[k + 1] = acc / total;
    }
    u.front() = 0.0;
    u.back() = 1.0;
    return u;
}

std::vector<double> centers_from_boundaries(const std::vector<double>& boundaries) {
    if (boundaries.size() < 2) return {};
    std::vector<double> centers(boundaries.size() - 1, 0.0);
    for (std::size_t k = 0; k + 1 < boundaries.size(); ++k) centers[k] = 0.5 * (boundaries[k] + boundaries[k + 1]);
    return centers;
}

std::vector<double> uniform_boundaries(std::size_t N) {
    std::vector<double> u(N + 1, 0.0);
    const double den = std::max<std::size_t>(1, N);
    for (std::size_t k = 0; k <= N; ++k) u[k] = static_cast<double>(k) / den;
    return u;
}

std::vector<double> clamp_and_dedup_sorted(std::vector<double> v, double lo, double hi) {
    for (auto& x : v) x = std::clamp(x, lo, hi);
    std::sort(v.begin(), v.end());
    v.erase(std::unique(v.begin(), v.end(), [](double a, double b) { return std::abs(a - b) < 1e-9; }), v.end());
    return v;
}

std::vector<CurveSeries> reparameterize_curves_on_u(const std::vector<CurveSeries>& curves,
                                                    const std::vector<double>& uCentersSrc,
                                                    const std::vector<double>& uCentersDst) {
    std::vector<CurveSeries> out;
    out.reserve(curves.size());
    for (const auto& c : curves) {
        CurveSeries r;
        r.kind = c.kind;
        r.times = uCentersDst;
        r.values = resample_linear(uCentersSrc, c.values, uCentersDst);
        if (!c.raw_score.empty() && c.raw_score.size() == c.values.size()) {
            r.raw_score = resample_linear(uCentersSrc, c.raw_score, uCentersDst);
        }
        out.push_back(std::move(r));
    }
    return out;
}

} // namespace

FilmEngine::FilmEngine(std::string databasePath)
    : store_(std::make_unique<SQLiteStore>(std::move(databasePath))) {
    store_->initialize();
}

FilmMetrics FilmEngine::compute_metrics(const std::vector<CurveSeries>& curves) const {
    FilmMetrics metrics;
    for (const auto& curve : curves) {
        if (curve.values.empty()) {
            continue;
        }
        const double sum = std::accumulate(curve.values.begin(), curve.values.end(), 0.0);
        const double mean = sum / static_cast<double>(curve.values.size());
        const double maxVal = *std::max_element(curve.values.begin(), curve.values.end());
        const double minVal = *std::min_element(curve.values.begin(), curve.values.end());
        const std::string prefix = curve_kind_to_string(curve.kind);
        metrics.scalars[prefix + ".mean"] = mean;
        metrics.scalars[prefix + ".max"] = maxVal;
        metrics.scalars[prefix + ".min"] = minVal;
    }
    return metrics;
}

AnalysisContext FilmEngine::buildAnalysisContext(const std::vector<FeatureSample>& samples,
                                                 const std::vector<FaceTrack>& tracks,
                                                 const std::vector<double>& cutTimesSec) const {
    AnalysisContext context;

    // Step 0: Setup unified 1Hz timeline
    // Contract: samples[k].timeSeconds = k + 0.5, sample_count = duration_seconds
    // So runtime_sec = samples.size() (strict 1Hz, no need for ceil)
    if (!samples.empty()) {
        const std::size_t N = samples.size();
        context.runtime_sec = static_cast<double>(N);
        context.t_center.resize(N);
        for (std::size_t k = 0; k < N; ++k) {
            context.t_center[k] = static_cast<double>(k) + 0.5;
        }
    }

    // Step 1: Run feature analyzers to expand simple samples into detailed features
    context.audioFeatures = soundAnalyzer_.analyze(samples);
    context.colorFeatures = colorAnalyzer_.analyze(samples);
    context.infoFeatures = infoAnalyzer_.analyze(samples);

    // Step 2: Build shot segments (doc 4.1.2 / 5.5): driven by cut detection timestamps (not motion peaks).
    context.shots = motionAnalyzer_.detectShots(samples, cutTimesSec);

    // Step 3: Build face features using real tracks from Bridge layer
    auto faceOut = faceAnalyzer_.analyze(samples, tracks);
    context.faceFeatures = faceOut.perSecond;
    context.dominantTrackId = faceOut.dominantTrackId;
    context.faceEnabled = faceOut.faceEnabled;

    // Step 3.5: Fill Info visual subfeatures that depend on dominant face (doc 4.4.3).
    // - ShotScaleNumeric: use dominant close-up weight when face present; otherwise 0 (wide/subject absent).
    // - FaceChange: 0/1 transition rate from face presence toggles.
    {
        const std::size_t N = std::min({context.infoFeatures.size(), context.faceFeatures.size(), context.t_center.size()});
        double prevPresence = 0.0;
        for (std::size_t k = 0; k < N; ++k) {
            const double presence = context.faceFeatures[k].facePresence;
            context.infoFeatures[k].shotScaleNumeric = (presence > 0.0) ? std::clamp(context.faceFeatures[k].closeUpWeight, 0.0, 1.0) : 0.0;
            context.infoFeatures[k].faceChange = (k == 0) ? 0.0 : std::abs(presence - prevPresence);
            prevPresence = presence;
        }
    }

    // Step 4: Build 1Hz motion vector from samples (O(N) direct mapping)
    // Contract: samples[k].timeSeconds = k + 0.5, so direct index mapping is safe
    if (!samples.empty()) {
        const std::size_t N = std::min(context.t_center.size(), samples.size());
        context.motion_1hz.resize(N, 0.0);
        context.camera_motion_1hz.resize(N, 0.0);
        context.object_motion_1hz.resize(N, 0.0);
        context.camera_motion_type_1hz.resize(N, 0);
        for (std::size_t k = 0; k < N; ++k) {
            context.motion_1hz[k] = samples[k].motionAmplitude;
            context.camera_motion_1hz[k] = samples[k].cameraMotion;
            context.object_motion_1hz[k] = samples[k].objectMotion;
            context.camera_motion_type_1hz[k] = static_cast<int32_t>(samples[k].cameraMotionType);
        }
    }

    return context;
}

FilmAnalysisResult FilmEngine::analyze_and_store(const FilmAnalysisRequest& request) {
    FilmAnalysisResult result;
    result.filmId = request.filmId;
    result.filmTitle = request.filmTitle;
    result.filePath = request.filePath;      // NEW: actual file path
    result.imdbId = request.imdbId;          // NEW: IMDb ID (optional)
    result.source = request.source;          // NEW: data source
    result.year = request.year;
    result.director = request.director;
    result.durationSeconds = request.durationSeconds;

    // Store face tracking info
    result.faceTracks = request.tracks;

    // Step 1: Build complete AnalysisContext from FeatureSamples and FaceTracks
    AnalysisContext context = buildAnalysisContext(request.samples, request.tracks, request.cutTimes);

    // Extract face analysis metadata from context (already computed in buildAnalysisContext)
    result.dominantTrackId = context.dominantTrackId;
    result.hasDominantFace = context.faceEnabled;

    // Step 2: Build curves using complete context
    result.curves = curveEngine_.build_curves(context);

    // Step 2.5: Build metadata overlays (display-only, never enters PCA/Diagnostics)
    result.overlays = metadataEngine_.build_overlays(context, result.curves);

    // Step 2.75: Persist shot segments (doc 5.5) for Shot grid visualization.
    result.shots = context.shots;

    // Step 3: Run diagnostics
    result.issues = diagnostics_.detect(result.curves, result.overlays);

    // Step 4: Detect multi-scale structure with adaptive window
    double windowSize = 10.0;  // fallback default
    // Compute adaptive window from Pace or Arousal raw_score lambda
    // λ = mean absolute slope of raw_score (Z-space), then W = adaptive_window_size(λ)
    for (const auto& curve : result.curves) {
        if (curve.kind == CurveKind::Arousal || curve.kind == CurveKind::Pace) {
            if (curve.raw_score.size() >= 2) {
                double lambda = 0.0;
                for (std::size_t k = 1; k < curve.raw_score.size(); ++k) {
                    lambda += std::abs(curve.raw_score[k] - curve.raw_score[k - 1]);
                }
                lambda /= static_cast<double>(curve.raw_score.size() - 1);
                windowSize = contract::adaptive_window_size(lambda);
            }
            break;  // Use first available curve (prefer Arousal, fallback to Pace)
        }
    }

    // arcscope.md 5.6: optional emotion-progress reparameterization u(t) for scene/sequence slicing.
    if ((request.options & StructureUseEmotionProgressReparam) != 0u) {
        const auto* arousal = findCurve(result.curves, CurveKind::Arousal);
        if (arousal && !arousal->values.empty()) {
            const std::size_t N = arousal->values.size();
            const double duration = std::max(1.0, result.durationSeconds);

            const auto uBoundSrc = emotion_progress_boundaries_from_arousal(arousal->values);  // length N+1
            const auto uCenterSrc = centers_from_boundaries(uBoundSrc);                        // length N

            const auto uBoundDst = uniform_boundaries(N);
            const auto uCenterDst = centers_from_boundaries(uBoundDst);

            const auto curvesU = reparameterize_curves_on_u(result.curves, uCenterSrc, uCenterDst);
            const double windowU = windowSize / duration; // W seconds -> W in u-axis (u in [0,1])

            const auto boundariesU = structureAnalyzer_.detectSceneBoundaries(curvesU, windowU);

            std::vector<double> tBoundSrc(N + 1, 0.0);
            for (std::size_t k = 0; k <= N; ++k) tBoundSrc[k] = static_cast<double>(k);

            // Map u-boundaries back to seconds boundaries.
            auto boundariesSec = resample_linear(uBoundSrc, tBoundSrc, boundariesU);
            boundariesSec = clamp_and_dedup_sorted(std::move(boundariesSec), 0.0, duration);
            result.scenes = structureAnalyzer_.scenesFromBoundaries(result.curves, boundariesSec);
        } else {
            result.scenes = structureAnalyzer_.detectScenes(result.curves, windowSize);
        }
    } else {
        result.scenes = structureAnalyzer_.detectScenes(result.curves, windowSize);
    }

    result.sequences = structureAnalyzer_.groupSequences(result.scenes, -1);
    result.acts = structureAnalyzer_.groupActs(result.sequences, result.durationSeconds);

    // Step 5: Compute metrics
    result.metrics = compute_metrics(result.curves);

    // Step 6: Store results
    store_->save_analysis(result);

    return result;
}

}  // namespace arcscope
