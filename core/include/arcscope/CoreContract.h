#pragma once

/**
 * CoreContract.h - ArcScope Core System Constants
 *
 * This file defines all contract-level constants for the ArcScope system.
 * These constants are PART OF THE SYSTEM CONTRACT and should NOT be changed
 * without understanding the implications for:
 *   - Historical database compatibility
 *   - Curve reproducibility across versions
 *   - UI timeline consistency
 *
 * VERSION: 1.0.0
 * Last Updated: 2025-12-13
 *
 * See: ArcScope_CoreContract.md for full specification
 */

namespace arcscope {
namespace contract {

// ============================================================================
// Time Invariant Constants (§ I.R3)
// ============================================================================

/**
 * TIME_SAMPLE_RATE_HZ - Unified timeline sample rate
 *
 * All curves and temporal data use 1Hz sampling:
 *   - Sample times: t[i] = i + 0.5 (center-aligned)
 *   - Sample count: N = ceil(runtime_seconds)
 *
 * Changing this breaks database compatibility and UI synchronization.
 */
constexpr double TIME_SAMPLE_RATE_HZ = 1.0;

/**
 * TIME_CENTER_OFFSET - Center-alignment offset for 1-second intervals
 *
 * Each sample represents interval [k, k+1) with center at k+0.5
 */
constexpr double TIME_CENTER_OFFSET = 0.5;

// ============================================================================
// Adaptive Window Constants (§ 3.3.1)
// ============================================================================

/**
 * Window computation for adaptive smoothing and local statistics.
 *
 * FORMULA: W* = W_BASE * exp(-W_LAMBDA_SCALE * lambda)
 *
 * Where:
 *   - lambda = mean_abs_slope_1hz(robust_z(raw_score))
 *   - W* is clamped to [W_MIN, W_MAX]
 *
 * RATIONALE:
 *   - High lambda (rapid changes) => smaller window (preserve detail)
 *   - Low lambda (smooth curves) => larger window (noise reduction)
 *
 * These constants determine the system's responsiveness to temporal dynamics.
 * Changing them affects curve smoothness and peak detection sensitivity.
 */

constexpr double W_BASE = 12.0;           // Base window size (seconds)
constexpr double W_LAMBDA_SCALE = 2.0;    // Lambda decay rate
constexpr double W_MIN = 4.0;             // Minimum window size (seconds)
constexpr double W_MAX = 30.0;            // Maximum window size (seconds)

/**
 * Compute adaptive window size from lambda
 *
 * @param lambda Mean absolute slope of robust_z(raw_score)
 * @return Window size in seconds, clamped to [W_MIN, W_MAX]
 */
inline double adaptive_window_size(double lambda) {
    const double w_star = W_BASE * std::exp(-W_LAMBDA_SCALE * lambda);
    return (w_star < W_MIN) ? W_MIN : (w_star > W_MAX) ? W_MAX : w_star;
}

// ============================================================================
// PCA and Normalization Constants (§ 4.1.4, 4.2.5)
// ============================================================================

/**
 * ROBUST_Z_MAD_SCALE - Robust standardization scaling factor
 *
 * arcscope.md §3.2 defines:
 *   Z(k) = (X_k - median(X)) / (MAD(X) + 1e-6)
 * with MAD(X) = median(|X - median(X)|).
 *
 * No sigma-alignment factor (e.g., 1.4826) is applied, by design.
 */
constexpr double ROBUST_Z_MAD_SCALE = 1.0;

/**
 * PCA_POWER_ITERATION_MAX - Maximum iterations for PC1 extraction
 *
 * Power iteration converges quickly for well-conditioned covariance matrices.
 * 64 iterations is overkill but guarantees convergence.
 */
constexpr int PCA_POWER_ITERATION_MAX = 64;

/**
 * PCA_CONVERGENCE_EPS - Numerical epsilon for PCA stability
 *
 * Used to prevent division by zero in covariance eigenvalue computation.
 */
constexpr double PCA_CONVERGENCE_EPS = 1e-24;

// ============================================================================
// Curve Gating and Threshold Constants
// ============================================================================

/**
 * PRESENCE_GATE_QUANTILE - Adaptive gate for face presence
 *
 * Frames with presence < quantile(presence, PRESENCE_GATE_QUANTILE) are zeroed.
 * This replaces the hard-coded 0.5 threshold with a scale-invariant gate.
 *
 * Range: [0.5, 0.8]
 *   - 0.6: Moderate gate (keeps ~40% of detected faces)
 *   - 0.7: Conservative gate (keeps ~30% of high-confidence faces)
 */
constexpr double PRESENCE_GATE_QUANTILE = 0.6;

/**
 * SHOT_CUT_QUANTILE - Adaptive threshold for shot segmentation
 *
 * Cut density values above quantile(cut_density, SHOT_CUT_QUANTILE) trigger cuts.
 */
constexpr double SHOT_CUT_QUANTILE = 0.70;

/**
 * SHOT_MOTION_JUMP_QUANTILE - Adaptive threshold for motion-based cuts
 *
 * Motion jumps above quantile(motion_jump, SHOT_MOTION_JUMP_QUANTILE) trigger cuts.
 */
constexpr double SHOT_MOTION_JUMP_QUANTILE = 0.75;

/**
 * SHOT_MIN_DURATION_SEC - Minimum shot duration to prevent spurious cuts
 */
constexpr double SHOT_MIN_DURATION_SEC = 1.0;

// ============================================================================
// Sigmoid and Output Layer Constants (§ 4.1.5)
// ============================================================================

/**
 * ROBUST_SIGMOID_CENTER - Center point for robust sigmoid mapping
 *
 * Maps Z-score 0 (median) to output 0.5.
 */
constexpr double ROBUST_SIGMOID_CENTER = 0.0;

/**
 * ROBUST_SIGMOID_STEEPNESS - Steepness of sigmoid transition
 *
 * Controls how sharply raw scores map to [0,1] range.
 * Higher values => sharper transitions (more binary-like)
 * Lower values => smoother transitions (more gradual)
 */
constexpr double ROBUST_SIGMOID_STEEPNESS = 1.0;

// ============================================================================
// Version Tracking
// ============================================================================

/**
 * CORE_CONTRACT_VERSION - Semantic version of this contract
 *
 * Increment when:
 *   - Major: Breaking changes to constants (database incompatible)
 *   - Minor: New constants added (backward compatible)
 *   - Patch: Documentation or non-functional changes
 *
 * This version should be stored in the database schema_version table.
 */
constexpr const char* CORE_CONTRACT_VERSION = "1.0.0";

} // namespace contract
} // namespace arcscope
