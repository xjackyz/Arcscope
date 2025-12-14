#include "arcscope/features/ColorFeatures.h"
#include "arcscope/Normalization.h"

#include <algorithm>
#include <array>
#include <cmath>
#include <cstddef>
#include <limits>
#include <numeric>
#include <vector>

namespace arcscope {
namespace features {

namespace {

// ---------- small utils ----------
inline double clamp01(double x) {
    if (!std::isfinite(x)) return 0.0;
    return std::clamp(x, 0.0, 1.0);
}

inline double safe_div(double a, double b, double fallback = 0.0) {
    return (std::abs(b) > 1e-12) ? (a / b) : fallback;
}

double warm_ratio_from_cam16_hist(const std::array<double,12>& hist) {
    double sum = 0.0;
    for (double v : hist) sum += v;
    if (!(sum > 1e-12)) return 0.0;

    // 12 bins, 30° each. Warm hues ≈ red/orange/yellow (330..90°).
    // This is a deterministic bin partition, not a tunable threshold.
    const int warmBins[] = {11, 0, 1, 2};
    double warm = 0.0;
    for (int b : warmBins) warm += hist[(size_t)b];
    return clamp01(warm / sum);
}

double grayness_from_cam16_m(double M) {
    return safe_div(1.0, std::max(1e-6, M), 0.0);
}

inline double l2_sq3(const std::array<double,3>& a, const std::array<double,3>& b) {
    double dx = a[0]-b[0];
    double dy = a[1]-b[1];
    double dz = a[2]-b[2];
    return dx*dx + dy*dy + dz*dz;
}

// Gentle temporal de-jitter for discrete states
// Rule: if i is different from both neighbors but neighbors agree -> snap to neighbor
void temporal_majority_smooth(std::vector<int>& states) {
    if (states.size() < 3) return;
    for (size_t i = 1; i + 1 < states.size(); ++i) {
        int a = states[i-1];
        int b = states[i];
        int c = states[i+1];
        if (a == c && b != a) {
            states[i] = a;
        }
    }
}

// Viterbi algorithm for temporal smoothing with transition costs
// observations: [N x 3] array as (warmth, brightness, saturation)
// centroids: [K x 3] cluster centers
// transition_penalty: cost for switching states (higher = prefer staying in same state)
std::vector<int> viterbi_smooth(const std::vector<std::array<double,3>>& observations,
                                 const std::vector<std::array<double,3>>& centroids,
                                 double transition_penalty = 0.15) {
    const int N = static_cast<int>(observations.size());
    const int K = static_cast<int>(centroids.size());

    if (N == 0 || K == 0) return {};
    if (K == 1) return std::vector<int>(N, 0);

    // DP tables: cost[t][k] = min cost to reach state k at time t
    std::vector<std::vector<double>> cost(N, std::vector<double>(K, 0.0));
    std::vector<std::vector<int>> backptr(N, std::vector<int>(K, -1));

    // Initialize: t=0, cost = emission cost only
    for (int k = 0; k < K; ++k) {
        cost[0][k] = l2_sq3(observations[0], centroids[k]);
    }

    // Forward pass
    for (int t = 1; t < N; ++t) {
        for (int k = 0; k < K; ++k) {
            double emission = l2_sq3(observations[t], centroids[k]);
            double best_cost = std::numeric_limits<double>::infinity();
            int best_prev = 0;

            for (int prev_k = 0; prev_k < K; ++prev_k) {
                // Transition cost: penalty if switching states
                double transition = (prev_k == k) ? 0.0 : transition_penalty;
                double total = cost[t-1][prev_k] + transition + emission;

                if (total < best_cost) {
                    best_cost = total;
                    best_prev = prev_k;
                }
            }

            cost[t][k] = best_cost;
            backptr[t][k] = best_prev;
        }
    }

    // Backtrack to find optimal path
    std::vector<int> path(N);

    // Find best final state
    int best_final = 0;
    double best_final_cost = cost[N-1][0];
    for (int k = 1; k < K; ++k) {
        if (cost[N-1][k] < best_final_cost) {
            best_final_cost = cost[N-1][k];
            best_final = k;
        }
    }

    path[N-1] = best_final;
    for (int t = N-2; t >= 0; --t) {
        path[t] = backptr[t+1][path[t+1]];
    }

    return path;
}

} // namespace

// ---------------- ColorAnalyzer ----------------

std::vector<ColorFeatures> ColorAnalyzer::analyze(const std::vector<FeatureSample>& samples) const {
    std::vector<ColorFeatures> result;
    if (samples.empty()) return result;
    result.reserve(samples.size());

    const int N = static_cast<int>(samples.size());

    // ColorState: KMeans in [Warmth, Brightness(J), Saturation(M)] space (doc 4.3.5)
    int K = std::clamp((int)std::lround(std::sqrt((double)N) / 3.0 + 2.0), 2, 6);
    K = std::min(K, N);
    std::vector<int> colorStates = clusterColorStates(samples, K);
    if (colorStates.size() != samples.size()) colorStates.assign(samples.size(), 0);

    for (size_t i = 0; i < samples.size(); ++i) {
        ColorFeatures cf;
        cf.timeSeconds = samples[i].timeSeconds;
        const bool hasCam16 = std::isfinite(samples[i].cam16J) && samples[i].cam16J > 0.0;

        cf.brightness = hasCam16 ? samples[i].cam16J : 0.0;
        cf.saturation = hasCam16 ? samples[i].cam16M : 0.0;
        cf.warmth     = hasCam16 ? warm_ratio_from_cam16_hist(samples[i].cam16HueHist) : 0.0;
        cf.grayness   = hasCam16 ? grayness_from_cam16_m(samples[i].cam16M) : 0.0;

        cf.colorHarmony  = hasCam16 ? clamp01(samples[i].cam16Harmony) : 0.0;
        cf.colorContrast = hasCam16 ? std::max(0.0, samples[i].cam16DeltaE) : 0.0;

        // Hue distribution is supplied by Bridge; it is normalized there.
        cf.hueDistribution.assign(samples[i].cam16HueHist.begin(), samples[i].cam16HueHist.end());

        // Discrete state
        cf.colorState = colorStates[i];

        result.push_back(std::move(cf));
    }

    return result;
}

std::vector<int> ColorAnalyzer::clusterColorStates(const std::vector<FeatureSample>& samples, int K) const {
    const int N = static_cast<int>(samples.size());
    if (N <= 0) return {};

    // Don't cluster if sample size is too small - insufficient data for meaningful states
    // Minimum ~12 samples ensures at least a few seconds per state
    if (N < 12) return std::vector<int>(N, 0);

    K = std::clamp(K, 1, N);

    // Features (doc 4.3.5): [Warmth, Brightness(J), Saturation(M)] in CAM16 space.
    // We robust-standardize columns to remove unit/scale effects without hard-coded weights.
    std::vector<double> warmths; warmths.reserve(N);
    std::vector<double> brs;     brs.reserve(N);
    std::vector<double> sats;    sats.reserve(N);

    for (const auto& s : samples) {
        const bool hasCam16 = std::isfinite(s.cam16J) && s.cam16J > 0.0;
        warmths.push_back(hasCam16 ? warm_ratio_from_cam16_hist(s.cam16HueHist) : 0.0);
        brs.push_back(hasCam16 ? s.cam16J : 0.0);
        sats.push_back(hasCam16 ? s.cam16M : 0.0);
    }

    auto zW = arcscope::robust_z(warmths);
    auto zB = arcscope::robust_z(brs);
    auto zS = arcscope::robust_z(sats);

    std::vector<std::array<double,3>> observations;
    observations.reserve(N);
    for (int i = 0; i < N; ++i) observations.push_back({zW[(size_t)i], zB[(size_t)i], zS[(size_t)i]});

    // Initialize centroids by quantiles (deterministic)
    std::vector<std::array<double,3>> centroids;
    centroids.reserve(K);

    for (int k = 0; k < K; ++k) {
        double q = static_cast<double>(k + 1) / static_cast<double>(K + 1);
        std::array<double,3> c = { quantile(zW, q), quantile(zB, q), quantile(zS, q) };
        centroids.push_back(c);
    }

    // De-duplicate centroids lightly (if quantiles collapse)
    for (int i = 1; i < K; ++i) {
        if (l2_sq3(centroids[i], centroids[i-1]) < 1e-6) {
            centroids[i][0] = clamp01(centroids[i][0] + 0.01 * (double)i);
            centroids[i][1] = clamp01(centroids[i][1] + 0.005 * (double)i);
            centroids[i][2] = clamp01(centroids[i][2] + 0.008 * (double)i);
        }
    }

    std::vector<int> assign(N, 0);

    const int maxIter = 25;
    for (int iter = 0; iter < maxIter; ++iter) {
        // Assignment
        for (int i = 0; i < N; ++i) {
            const std::array<double,3> x = observations[(size_t)i];
            double best = std::numeric_limits<double>::infinity();
            int bestK = 0;
            for (int k = 0; k < K; ++k) {
                double d = l2_sq3(x, centroids[k]);
                if (d < best) {
                    best = d;
                    bestK = k;
                }
            }
            assign[i] = bestK;
        }

        // Update
        std::vector<std::array<double,3>> newC(K, {0.0,0.0,0.0});
        std::vector<int> cnt(K, 0);

        for (int i = 0; i < N; ++i) {
            int k = assign[i];
            newC[k][0] += observations[(size_t)i][0];
            newC[k][1] += observations[(size_t)i][1];
            newC[k][2] += observations[(size_t)i][2];
            cnt[k] += 1;
        }

        // Handle empty clusters: reinit to farthest point from its assigned centroid
        for (int k = 0; k < K; ++k) {
            if (cnt[k] > 0) continue;

            // Find the point with largest error to its centroid
            double worst = -1.0;
            int worstIdx = 0;
            for (int i = 0; i < N; ++i) {
                double d = l2_sq3(observations[(size_t)i], centroids[assign[i]]);
                if (d > worst) {
                    worst = d;
                    worstIdx = i;
                }
            }
            newC[k] = observations[(size_t)worstIdx];
            cnt[k] = 1;
        }

        // Compute shift for convergence
        double shift = 0.0;
        for (int k = 0; k < K; ++k) {
            std::array<double,3> updated = {
                newC[k][0] / (double)cnt[k],
                newC[k][1] / (double)cnt[k],
                newC[k][2] / (double)cnt[k]
            };
            shift += l2_sq3(updated, centroids[k]);
            centroids[k] = updated;
        }

        if (shift < 1e-6) {
            break;
        }
    }

    // Adaptive Viterbi smoothing: compute transition penalty based on data scale (z-space)
    std::vector<double> emissions;
    emissions.reserve(N);
    for (int i = 0; i < N; ++i) {
        double dist = l2_sq3(observations[i], centroids[assign[i]]);
        emissions.push_back(dist);
    }

    // Use median emission distance as scale for transition penalty
    std::vector<double> emissionsSorted = emissions;
    std::nth_element(emissionsSorted.begin(),
                     emissionsSorted.begin() + N/2,
                     emissionsSorted.end());
    double medianEmission = emissionsSorted[N/2];

    // Penalty = alpha * scale, where alpha controls preference for temporal stability
    // alpha ~ 0.8 means switching costs about 80% of typical emission distance
    const double alpha = 0.8;
    double transition_penalty = alpha * std::max(medianEmission, 0.01); // prevent zero penalty

    assign = viterbi_smooth(observations, centroids, transition_penalty);

    // Final pass: temporal majority smoothing to remove single-frame jitter
    // This is cheap and catches isolated single-frame state flips that Viterbi might miss
    temporal_majority_smooth(assign);

    return assign;
}

}  // namespace features
}  // namespace arcscope
