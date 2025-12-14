#pragma once

#include <stdint.h>
#include <stddef.h>

#ifdef __cplusplus
extern "C" {
#endif

// ============================================================================
// Analysis options (arcscope.md 5.6)
// ============================================================================
// These flags affect FilmEngine behavior during analysis only.
// They are stored in FilmAnalysisRequest.options (core) and do not change the 1Hz contract.
//
// NOTE: Keep values in sync with core/include/arcscope/FilmTypes.h (FilmAnalysisOptions).
#define ARCSCOPE_ANALYZE_OPT_STRUCTURE_U_REPARAM (1u << 0)

// Progress callback invoked during long running operations (e.g., video analysis).
typedef void (*ArcScopeProgressCallback)(float progress, void* user_data);

// ============================================================================
// ArcScopeBridgeSample - Raw Observation Layer
// ============================================================================
//
// LEGAL CONTRACT:
// This struct represents ONLY raw per-second observations extracted from
// video/audio/subtitle streams. It MUST NOT contain:
//   - Fused/derived quantities (e.g., arousalProxy, dialogueDensity)
//   - Semantic interpretations (e.g., sceneType, tension)
//   - Normalized/scaled values (FilmEngine responsibility)
//
// TIME INVARIANT:
//   - Each sample represents exactly 1 second: interval [k, k+1)
//   - time_seconds = k + 0.5 (canonical mid-point)
//   - Sample count = film duration in seconds (integer)
//
// DATA SEMANTICS:
//   - All fields are measurements of the physical world
//   - Values are unnormalized (preserve original scale/units)
//   - Missing/invalid data represented as 0.0 or -1.0 (see field docs)
//
// See: ArcScope_CoreContract.md § II.R4 for enforcement rules
// ============================================================================

typedef struct {
    // -------------------------------------------------------------------------
    // Temporal Coordinate (required)
    // -------------------------------------------------------------------------
    double time_seconds;              // Mid-point of 1-second interval [k, k+1)

    // -------------------------------------------------------------------------
    // Visual Observations (from AVFoundation video frames)
    // -------------------------------------------------------------------------
    // Number of detected hard cuts whose timestamp t satisfies floor(t) == k.
    // (Windowed CutDensity CD(t) is derived in Core from shot segments + adaptive window W.)
    double cut_density;
    double motion_amplitude;          // Average pixel displacement (0.0–1.0+)
    double brightness;                // Mean luminance (0.0–1.0)
    double saturation;                // Mean chroma intensity (0.0–1.0)
    double warmth;                    // Red/blue ratio (-1.0=cool, +1.0=warm)
    double grayness;                  // Chroma absence (0.0=vivid, 1.0=gray)
    double avg_hue;                   // Dominant hue angle [0,360), or -1.0 if invalid

    // -------------------------------------------------------------------------
    // Color Science Pipeline (Rec.2020 linear -> CAM16-UCS)
    // -------------------------------------------------------------------------
    // Per-second aggregates computed from a fixed analysis frame rate (10fps) with
    // uniform downsampling to fixed height (720p), per arcscope.md §3.4.
    //
    // These are observation-level descriptors in a perceptual uniform space:
    //   - CAM16 lightness J (0..100)
    //   - CAM16 colorfulness M (0..~)
    //   - hue angle h (0..360)
    //   - CAM16-UCS coordinates (J', a', b') for ΔE and clustering
    //   - ΔE jump between consecutive seconds in CAM16-UCS
    //   - 12-bin hue histogram (in CAM16 hue angle), normalized to sum=1
    //   - harmony score (K-medians in CAM16-UCS), normalized to [0,1]
    double cam16_j;
    double cam16_m;
    double cam16_h;
    double cam16_jp;
    double cam16_ap;
    double cam16_bp;
    double cam16_delta_e;
    double cam16_hue_hist[12];
    double cam16_harmony;

    // -------------------------------------------------------------------------
    // Audio Observations (from AVFoundation audio frames)
    // -------------------------------------------------------------------------
    double audio_rms;                 // RMS amplitude (0.0–1.0)
    double audio_transient;           // Onset density (attacks/sec, 0.0–10.0+)
    double spectral_balance;          // High/low frequency ratio (0.0–1.0)

    // -------------------------------------------------------------------------
    // Audio Analysis (FFT-derived, still observation-level)
    // -------------------------------------------------------------------------
    // These are computed in Bridge (AudioAnalysisEngine) and passed through as
    // physical measurements / signal descriptors (NOT semantic FilmEngine outputs).
    double audio_loudness;            // LUFS-like (dB), A-weighted spectrum integration
    double audio_spectral_flux;       // Spectral flux (0.0–1.0+)
    double audio_spectral_centroid;   // Spectral centroid (Hz)
    double audio_zero_crossing_rate;  // ZCR (0.0–1.0)
    double audio_bpm;                 // Tempo (BPM)
    double audio_tempo_confidence;    // [0,1]
    double audio_rhythm_energy;       // Rhythm drive proxy [0,1]
    int audio_estimated_key;          // 0=C..11=B, -1 unknown
    double audio_harmonic_tension;    // Harmonic tension [0,1]
    double audio_dialogue_clarity;    // Speech-band SNR proxy [0,1]
    double audio_speech_probability;  // Speech presence [0,1]
    double audio_chroma[12];          // 12-bin chroma vector (sum≈1), or all zeros if unavailable

    // -------------------------------------------------------------------------
    // Audio role descriptors (arcscope.md 4.2.6 / 5.6)
    // -------------------------------------------------------------------------
    double audio_stereo_width;        // Mid/side energy ratio [0,1]
    double audio_reverb_amount;       // Reverb proxy [0,1]
    double audio_dialogue_energy;     // Speech-band power (arcscope.md 4.2.6), raw units
    double audio_music_energy;        // Non-speech tonal power (arcscope.md 4.2.6), raw units
    double audio_effects_energy;      // Residual non-speech power (arcscope.md 4.2.6), raw units
    double audio_dialogue_dominance;  // DialogueEnergy / total [0,1]
    double audio_music_dominance;     // MusicEnergy / total [0,1]
    double audio_effects_dominance;   // EffectsEnergy / total [0,1]

    // -------------------------------------------------------------------------
    // Camera vs object motion (arcscope.md 5.6 Photography/Color)
    // -------------------------------------------------------------------------
    double camera_motion;             // Global motion magnitude proxy
    double object_motion;             // Residual / local motion magnitude proxy
    int camera_motion_type;           // 0=none,1=pan,2=tilt,3=push,4=pull,5=track

    // -------------------------------------------------------------------------
    // Subtitle/Text Observations (from SRT/VTT parsing)
    // -------------------------------------------------------------------------
    double subtitle_density;          // Words per second (0.0–10.0+)
    double exposition_probability;    // Heuristic: likely expository dialogue (0.0–1.0)

    // -------------------------------------------------------------------------
    // NLP Features (from Apple NaturalLanguage processing)
    // -------------------------------------------------------------------------
    // These fields are populated by SubtitleNLPProcessor in Bridge layer.
    // If no subtitles or NLP unavailable, all fields will be 0.
    // Core's InfoAnalyzer will fall back to proxy features if these are 0.
    int nlp_word_count;               // True word count from NL tokenizer (0+)
    double nlp_sentence_complexity;   // Sentence complexity score (0.0–1.0)
    double nlp_new_concept_ratio;     // Ratio of new concepts (0.0–1.0)
    double nlp_entity_density;        // Named entity density (0.0–1.0)

    // -------------------------------------------------------------------------
    // Information Observations (Bridge-level coarse estimate - USE WITH CAUTION)
    // -------------------------------------------------------------------------
    // ⚠️ WARNING: info_density is a COARSE BRIDGE-LEVEL ESTIMATE.
    //    It should be treated as METADATA ONLY, NOT used for Core's PCA/InfoAnalyzer.
    //    Using it in Core will cause:
    //      - Non-reproducible results across different macOS/model versions
    //      - PCA axis drift when Bridge's estimation logic changes
    //      - Violation of the "atomic observations only" contract
    //
    // RECOMMENDED USAGE:
    //    - Bridge: Populate as a rough diagnostic hint (optional)
    //    - Core: DO NOT use in PCA, InfoAnalyzer, or any normalization
    //    - Core: Calculate true InfoDensity from atomic features (subtitle_density,
    //            visual_entropy, audio_transient, etc.) via PCA
    //
    // If you need info density in Core, derive it from:
    //    subtitle_density + visual_motion + audio_complexity
    //    (see § IV.R8 of ArcScope_CoreContract.md)
    double info_density;              // Bridge-level info density estimate (0.0–1.0)

    // -------------------------------------------------------------------------
    // Face/Expression Observations (from CoreML Vision + AU model)
    // -------------------------------------------------------------------------
    double face_presence;             // Total face area / frame area (0.0–1.0)
    double face_arousal;              // Facial action unit intensity (0.0–1.0)
    double face_valence;              // Facial valence (negative/positive) (-1.0–1.0)
    double face_confidence;           // Detection confidence (0.0–1.0)

    // -------------------------------------------------------------------------
    // FORBIDDEN FIELDS (must never appear here):
    // -------------------------------------------------------------------------
    // ❌ dialogue_density      → FilmEngine fusion (audio + subtitle)
    // ❌ arousal_proxy         → FilmEngine fusion (motion + audio + face)
    // ❌ pace                  → FilmEngine PCA output
    // ❌ tension               → FilmEngine derived metric
    // ❌ scene_type            → FilmEngine semantic classification
    // ❌ verbal_load           → FilmEngine InfoAnalyzer output
    // ❌ visual_load           → FilmEngine InfoAnalyzer output
    // ❌ event_load            → FilmEngine InfoAnalyzer output
    //
    // Note: info_density is ALLOWED as a Bridge-level coarse estimate,
    //       but VerbalLoad/VisualLoad/EventLoad are FilmEngine-only.
    //
    // Adding forbidden fields violates § II.R5 of ArcScope_CoreContract.md
    // -------------------------------------------------------------------------

} ArcScopeBridgeSample;

// ============================================================================
// ArcScopeBridgeFaceTrack - Raw Face Tracking Data
// ============================================================================
//
// LEGAL CONTRACT:
// This struct represents continuous face tracking across a film's timeline.
// It is a RAW OBSERVATION from CoreML Vision + Facial AU models.
//
// MEMORY OWNERSHIP:
//   - Pointers (time_points, face_areas, etc.) are allocated by Bridge
//   - When used internally by arcscope_analyze_film():
//     * Bridge allocates, Core copies, Bridge deallocates automatically
//   - When used with arcscope_bridge_analyze():
//     * CALLER MUST use arcscope_bridge_free_tracks() to free memory
//     * Failure to call free_tracks() will cause memory leaks
//
// LIFECYCLE (arcscope_analyze_film):
//   1. Bridge allocates tracks during video analysis
//   2. Core copies track data in analyze_and_store()
//   3. Bridge deallocates tracks immediately after Core returns
//
// LIFECYCLE (arcscope_bridge_analyze with external tracks):
//   1. Caller allocates tracks (e.g., via ConvertObjCTracksToC)
//   2. Caller passes tracks to arcscope_bridge_analyze()
//   3. Core copies track data
//   4. Caller MUST call arcscope_bridge_free_tracks() to deallocate
//
// DATA SEMANTICS:
//   - Each track represents a single person's continuous appearance
//   - time_points[i] = when this face was observed (seconds)
//   - All arrays have length = point_count
//
// See: ArcScope_CoreContract.md § III.R7 for memory rules
// ============================================================================

typedef struct {
    int track_id;                     // Unique ID for this face track
    double appearance_duration;       // Total time this face appears (seconds)
    double avg_face_area;             // Mean face area / frame area (0.0–1.0)
    double avg_confidence;            // Mean detection confidence (0.0–1.0)
    double score;                     // Track quality score (higher = better)

    // Per-frame arrays (all of length point_count):
    double* time_points;              // Observation timestamps (seconds)
    double* face_areas;               // Face area at each time (0.0–1.0)
    double* confidences;              // Detection confidence at each time (0.0–1.0)
    double* valences;                 // Facial valence at each time (-1.0–1.0)
    double* arousals;                 // Facial arousal at each time (0.0–1.0)
    double* yaws;                     // Head yaw (radians), 0 = frontal
    double* rolls;                    // Head roll (radians), 0 = upright
    double* downcast_ratios;          // Landmark proxy in face coords: (nose-eye)/(mouth-eye)
    size_t point_count;               // Number of observations in this track
} ArcScopeBridgeFaceTrack;

// ============================================================================
// ArcScopeBridgeIssue - Diagnostic/Semantic Output (FilmEngine → Bridge)
// ============================================================================
//
// LEGAL CONTRACT:
// This struct is SEMANTIC OUTPUT from FilmEngine, not Bridge input.
// It represents a diagnosed issue/anomaly in the film's structure.
//
// DATA FLOW:
//   FilmEngine → SQLiteStore → Bridge → UI
//   (NOT: Bridge → FilmEngine)
//
// MEMORY OWNERSHIP:
//   - Strings (type, explanation) are allocated by Bridge on read
//   - Caller must use bridge_free_detail() to deallocate
//
// See: ArcScope_CoreContract.md § II (Semantic Boundary)
// ============================================================================

typedef struct {
    const char* type;                 // Issue category (e.g., "pacing_stall")
    double start_time;                // Issue start (seconds)
    double end_time;                  // Issue end (seconds)
    double severity;                  // Issue severity (0.0–1.0)
    const char* explanation;          // Human-readable description
} ArcScopeBridgeIssue;

// ============================================================================
// ArcScopeBridgeCurve - Time-Series Data (Raw or Derived)
// ============================================================================
//
// LEGAL CONTRACT:
// This struct represents a time-series curve stored in the database.
// It can contain EITHER raw observations OR derived FilmEngine outputs.
//
// TIME INVARIANT:
//   - count MUST equal film duration in seconds (integer)
//   - times[i] = i + 0.5 (for i = 0 to count-1)
//   - fps is always 1.0 (enforced by SQLiteStore)
//
// MEMORY OWNERSHIP:
//   - times/values arrays are allocated by Bridge on read
//   - Caller must use bridge_free_curve() to deallocate
//
// CURVE TYPES:
//   Raw observations: motion, audio_rms, brightness, etc.
//   Derived quantities: pace, tension, arousal_proxy, etc.
//
// See: ArcScope_CoreContract.md § I.R3 (Time Invariant)
// ============================================================================

typedef struct {
    const char* kind;                 // Curve identifier (e.g., "motion", "pace")
    float* times;                     // Timestamps (seconds, always i+0.5)
    float* values;                    // Curve values (range depends on kind)
    size_t count;                     // Number of points = film duration (seconds)
} ArcScopeBridgeCurve;

// ============================================================================
// Metadata Overlays - Display-only Tracks (NEVER enter PCA/Diagnostics)
// ============================================================================

typedef enum {
    ARCSCOPE_OVERLAY_F32 = 0,
    ARCSCOPE_OVERLAY_I32 = 1
} ArcScopeOverlayValueType;

typedef struct {
    const char* kind;                 // Overlay identifier (e.g., "color_warmth", "bpm", "hue_hist12")
    float* times;                     // Timestamps (seconds, always i+0.5)
    size_t count;                     // Number of points = film duration (seconds)
    int channels;                     // 1 for scalar; >1 for multi-channel
    ArcScopeOverlayValueType value_type;
    float* values_f32;                // length = count * channels (if value_type==F32)
    int32_t* values_i32;              // length = count * channels (if value_type==I32)
} ArcScopeBridgeOverlay;

// ============================================================================
// ArcScopeFilmSummary - Film Metadata (User-Facing)
// ============================================================================
//
// LEGAL CONTRACT:
// This struct contains user-facing film metadata for lists/selection UI.
//
// ID SEMANTICS:
//   - film_id is the LOGICAL ID (stable across re-imports)
//   - Example: "blade_runner_1982", "2001_a_space_odyssey"
//   - This is NOT the database primary key (see § III.R6)
//
// MEMORY OWNERSHIP:
//   - Strings are allocated by Bridge on read
//   - Caller must use bridge_free_film_list() to deallocate
//
// See: ArcScope_CoreContract.md § III.R6 (ID Duality)
// ============================================================================

typedef struct {
    const char* film_id;              // Logical ID (permanent, user-visible)
    const char* title;                // Film title
    int year;                         // Release year
    const char* director;             // Director name
    double duration_seconds;          // Film duration
} ArcScopeFilmSummary;

typedef struct {
    ArcScopeFilmSummary* items;
    size_t count;
} ArcScopeFilmList;

// ============================================================================
// Structure Segments - 4-Layer Hierarchical Structure
// ============================================================================

typedef struct {
    double start_time;
    double end_time;
    double duration;
    double avg_motion;
} ArcScopeShotSegment;

typedef struct {
    double start_time;
    double end_time;
    const char* label;              // Scene type (e.g., "BuildUp", "Peak")
} ArcScopeSceneSegment;

typedef struct {
    double start_time;
    double end_time;
} ArcScopeSequenceSegment;

typedef struct {
    int act_number;
    double start_time;
    double end_time;
} ArcScopeActSegment;

typedef struct {
    const char* film_id;
    const char* title;
    int year;
    const char* director;
    const char* file_path;           // May be NULL if unknown
    double duration_seconds;
    ArcScopeBridgeCurve* curves;
    size_t curve_count;
    ArcScopeBridgeOverlay* overlays;
    size_t overlay_count;
    ArcScopeBridgeIssue* issues;
    size_t issue_count;
    ArcScopeShotSegment* shots;
    size_t shot_count;
    ArcScopeSceneSegment* scenes;
    size_t scene_count;
    ArcScopeSequenceSegment* sequences;
    size_t sequence_count;
    ArcScopeActSegment* acts;
    size_t act_count;
} ArcScopeFilmDetail;

// Analyzer lifecycle -------------------------------------------------------

typedef struct ArcScopeBridgeAnalyzer ArcScopeBridgeAnalyzer;

ArcScopeBridgeAnalyzer* arcscope_bridge_create(const char* database_path);
void arcscope_bridge_destroy(ArcScopeBridgeAnalyzer* analyzer);

// Video analysis -----------------------------------------------------------

int arcscope_analyze_film(const char* input_path,
                          const char* db_path,
                          const char* film_id,
                          const char* title,
                          int year,
                          const char* director,
                          ArcScopeProgressCallback progress_cb,
                          void* user_data);

int arcscope_analyze_film_with_options(const char* input_path,
                                      const char* db_path,
                                      const char* film_id,
                                      const char* title,
                                      int year,
                                      const char* director,
                                      uint32_t options,
                                      ArcScopeProgressCallback progress_cb,
                                      void* user_data);

int arcscope_bridge_analyze(ArcScopeBridgeAnalyzer* analyzer,
                            const char* film_id,
                            const char* film_title,
                            const ArcScopeBridgeSample* samples,
                            size_t sample_count,
                            const ArcScopeBridgeFaceTrack* tracks,
                            size_t track_count);

int arcscope_bridge_analyze_with_options(ArcScopeBridgeAnalyzer* analyzer,
                                        const char* film_id,
                                        const char* film_title,
                                        const ArcScopeBridgeSample* samples,
                                        size_t sample_count,
                                        const ArcScopeBridgeFaceTrack* tracks,
                                        size_t track_count,
                                        uint32_t options);

// Data access --------------------------------------------------------------

int arcscope_load_film_detail(const char* db_path,
                              const char* film_id,
                              ArcScopeFilmDetail* out_detail);

int arcscope_list_films(const char* db_path, ArcScopeFilmList* out_list);

// Reference library --------------------------------------------------------

// List films that are marked as references (references_library table).
int arcscope_list_reference_films(const char* db_path, ArcScopeFilmList* out_list);

// Mark/unmark a film as reference.
// - film_id: logical_id (e.g., "blade-runner-1982")
// - is_reference: 1 to add/update, 0 to remove
// - note: optional note stored in references_library.note (may be NULL)
int arcscope_set_reference_film(const char* db_path,
                                const char* film_id,
                                int is_reference,
                                const char* note);

void arcscope_bridge_free_film_list(ArcScopeFilmList* list);
void arcscope_bridge_free_detail(ArcScopeFilmDetail* detail);

// Memory management for face tracks ----------------------------------------

// Free memory for face tracks (if allocated externally)
// ONLY use this if you created tracks outside of arcscope_analyze_film()
// (e.g., via ConvertObjCTracksToC or custom allocation)
void arcscope_bridge_free_tracks(ArcScopeBridgeFaceTrack* tracks, size_t track_count);

#ifdef __cplusplus
}
#endif
