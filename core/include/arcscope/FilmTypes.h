#pragma once

#include <string>
#include <unordered_map>
#include <array>
#include <cstdint>
#include <vector>
#include <optional>

namespace arcscope {

// ========== Basic Feature Sample (Bridge layer interface) ==========
struct FeatureSample {
    double timeSeconds{0.0};
    // Cut events aggregated in Bridge at analysis FPS (e.g., 10fps), then bucketed into 1-second bins.
    // Semantics: number of detected hard cuts whose timestamp t satisfies floor(t) == k.
    // (This is NOT the windowed CutDensity CD(t); Core derives CD(t) from shot segments + adaptive window W.)
    double cutDensity{0.0};
    double motionAmplitude{0.0};  // normalized optical flow magnitude
    double cameraMotion{0.0};     // global motion component (normalized)
    double objectMotion{0.0};     // residual/local motion component (normalized)
    int cameraMotionType{0};      // 0=none,1=pan,2=tilt,3=push,4=pull,5=track (observation-level)
    double audioRms{0.0};
    double audioTransient{0.0};
    double spectralBalance{0.0};

    // FFT-derived audio observations (Bridge AudioAnalysisEngine)
    // Psychoacoustic loudness proxy in LUFS-like log units (dB), computed from A-weighted spectrum.
    // Stored as a raw descriptor; Core performs robust z + sigmoid normalization.
    double audioLoudness{0.0};
    double audioSpectralFlux{0.0};       // spectral flux (0..1+)
    double audioSpectralCentroid{0.0};   // spectral centroid (Hz)
    double audioZeroCrossingRate{0.0};   // ZCR (0..1)
    double audioBpm{0.0};                // BPM
    double audioTempoConfidence{0.0};    // [0,1]
    double audioRhythmEnergy{0.0};       // [0,1]
    int audioEstimatedKey{-1};           // 0..11, -1 unknown
    double audioHarmonicTension{0.0};    // [0,1]
    double audioDialogueClarity{0.0};    // [0,1]
    double audioSpeechProbability{0.0};  // [0,1]
    std::vector<double> audioChroma;     // 12-bin chroma vector, or empty if unavailable
    double audioStereoWidth{0.0};        // [0,1]
    double audioReverbAmount{0.0};       // [0,1]
    double audioDialogueEnergy{0.0};     // speech-band power (raw units)
    double audioMusicEnergy{0.0};        // non-speech tonal power (raw units)
    double audioEffectsEnergy{0.0};      // residual non-speech power (raw units)
    double audioDialogueDominance{0.0};  // [0,1]
    double audioMusicDominance{0.0};     // [0,1]
    double audioEffectsDominance{0.0};   // [0,1]

    // Color science observations (Rec.2020 linear -> CAM16-UCS)
    double cam16J{0.0};                 // J in [0,100]
    double cam16M{0.0};                 // M (colorfulness)
    double cam16H{0.0};                 // hue angle [0,360), -1 if invalid
    double cam16Jp{0.0};                // CAM16-UCS J'
    double cam16Ap{0.0};                // CAM16-UCS a'
    double cam16Bp{0.0};                // CAM16-UCS b'
    double cam16DeltaE{0.0};            // ΔE in CAM16-UCS between consecutive seconds
    std::array<double, 12> cam16HueHist{}; // 12-bin hue histogram (sum≈1), or all zeros
    double cam16Harmony{0.0};           // [0,1] harmony score from K-medians in CAM16-UCS

    double brightness{0.0};
    double saturation{0.0};
    double warmth{0.0};
    double grayness{0.0};
    double avgHue{0.0};           // Average hue angle [0, 360) from CAM16/Lab, -1 if invalid
    double subtitleDensity{0.0};
    double expositionProbability{0.0};
    double infoDensity{0.0};      // Bridge-level coarse estimate

    // NLP features (from Bridge's SubtitleNLPProcessor)
    // Level 0 (no NLP): all zeros, InfoAnalyzer uses fallback proxy
    // Level 1 (has NLP): populated from Apple NaturalLanguage
    int nlpWordCount{0};                   // True word count from tokenizer
    double nlpSentenceComplexity{0.0};     // [0,1] sentence complexity
    double nlpNewConceptRatio{0.0};        // [0,1] new concept ratio
    double nlpEntityDensity{0.0};          // [0,1] named entity density

    double facePresence{0.0};
    double faceValence{0.0};
    double faceArousal{0.0};
    double faceConfidence{0.0};   // Detection confidence
};

// ========== Detailed Audio Features ==========
struct AudioFeatures {
    double timeSeconds{0.0};            // center time k+0.5
    double rms{0.0};                    // Root mean square energy
    double loudness{0.0};               // Psychoacoustic loudness (LUFS-style / A-weighted)
    double spectralFlux{0.0};           // Frame-to-frame spectral change
    double spectralCentroid{0.0};       // Spectral brightness (raw)
    double spectralBalance{0.0};        // High-freq vs low-freq ratio
    double zeroCrossingRate{0.0};       // Transient indicator
    double bpm{0.0};                    // Estimated beats per minute
    double tempoConfidence{0.0};        // Tempo detection confidence [0,1]
    double rhythmEnergy{0.0};           // RhythmDrive proxy
    std::vector<double> chroma;         // 12-bin chroma vector
    int estimatedKey{-1};               // 0=C, 1=C#, ..., 11=B, -1=unknown
    double harmonicTension{0.0};        // HarmonicTension
    double dialogueClarity{0.0};        // SNR in speech band (300-3400Hz)
    double speechProbability{0.0};      // VAD-style speech presence [0,1]
    double stereoWidth{0.0};            // Mid/side energy ratio [0,1]
    double reverbAmount{0.0};           // Reverb proxy [0,1]
    double dialogueEnergy{0.0};         // speech-band power (raw units)
    double musicEnergy{0.0};            // non-speech tonal power (raw units)
    double effectsEnergy{0.0};          // residual non-speech power (raw units)
    double dialogueDominance{0.0};      // DialogueEnergy / total [0,1]
    double musicDominance{0.0};         // MusicEnergy / total [0,1]
    double effectsDominance{0.0};       // EffectsEnergy / total [0,1]
};

// ========== Detailed Color Features ==========
struct ColorFeatures {
    double timeSeconds{0.0};
    double brightness{0.0};             // CAM16 J (0..100)
    double saturation{0.0};             // CAM16 M (colorfulness)
    double warmth{0.0};                 // Warm-hue mass ratio from CAM16 hue histogram [0,1]
    double grayness{0.0};               // Grayness proxy (high=gray), derived from CAM16 M
    double colorHarmony{0.0};           // HarmonyScore from K-medians in CAM16-UCS [0,1]
    double colorContrast{0.0};          // ΔE jump in CAM16-UCS between consecutive seconds
    int colorState{0};                  // Discrete palette mode (KMeans cluster ID)
    std::vector<double> hueDistribution; // 12-bin hue histogram
};

// ========== Info/NLP Features ==========
struct InfoFeatures {
    double timeSeconds{0.0};
    int wordCount{0};                       // Words spoken/shown per second
    double speechDuty{0.0};                 // Percentage of time with speech [0,1]
    double conceptDensity{0.0};             // Unique concepts per second
    double newConceptRatio{0.0};            // New vs repeated concepts [0,1]
    double sentenceComplexity{0.0};         // Avg clause depth/complexity [0,1]
    double dialogueClarity{0.0};            // Speech clarity proxy [0,1] (doc 4.4.2)
    double expositionScore{0.0};            // Technical/expository content [0,1]
    double visualEntropy{0.0};              // Frame visual complexity
    double shotScaleNumeric{0.0};           // Shot scale: 0=ELS, 1=ECU
    double cameraMotionComplexity{0.0};     // Camera movement complexity
    double faceChange{0.0};                 // Face appearance/disappearance rate
    double soundEvents{0.0};                // Non-speech sound event density
    double musicChange{0.0};                // Music transition density
    double cutDensity{0.0};                 // CutDensity (doc 4.4.4), cuts per second
    double objectMotionJump{0.0};           // Object motion discontinuity
    // Aggregated per second (as per doc 4.4)
    double verbalLoad{0.0};
    double visualLoad{0.0};
    double eventLoad{0.0};
    double exposition{0.0};                 // ExpositionCurve/ExpoScore
};

// ========== Face Tracking ==========
struct FaceTrack {
    int trackId{-1};
    double appearanceDuration{0.0};     // Total seconds on screen
    double avgFaceArea{0.0};            // Mean face/frame ratio
    double avgConfidence{0.0};          // Mean detection confidence
    double score{0.0};                  // T_i * S_i * C_i
    std::vector<double> timePoints;     // Seconds where face appears
    std::vector<double> faceAreas;      // Face area at each time point
    std::vector<double> confidences;    // Detection confidence at each time point
    std::vector<double> valences;       // Valence at each time point
    std::vector<double> arousals;       // Arousal at each time point
    std::vector<double> yaws;           // Head yaw (radians)
    std::vector<double> rolls;          // Head roll (radians)
    std::vector<double> downcastRatios; // (nose-eye)/(mouth-eye) in face coords
};

struct FaceFeatures {
    double timeSeconds{0.0};
    int dominantTrackId{-1};            // ID of main character track (-1 if none)
    double facePresence{0.0};           // 0/1 or probability; used as mask
    double faceArea{0.0};               // Face area / frame area
    double faceValence{0.0};            // Emotional valence [-1,1]
    double faceArousal{0.0};            // A_* (model output)
    double expressionIntensity{0.0};    // E_*
    double faceAreaRatio{0.0};          // S_* = face_pixels / frame_pixels
    double closeUpWeight{0.0};          // Shot scale weighting [0,1]
    double confidence{1.0};             // optional quality
    double yaw{0.0};                    // head yaw (radians)
    double roll{0.0};                   // head roll (radians)
    double downcastRatio{0.5};          // landmark proxy ratio (0..1)
};

// ========== Shot/Scene/Sequence/Act Segments ==========
struct ShotSegment {
    double startTime{0.0};
    double endTime{0.0};
    double duration{0.0};
    double avgMotion{0.0};              // Average optical flow in shot
};

enum class SceneType {
    Unknown,
    LowEnergyBasin,
    BuildUp,
    HighEnergyPeak,
    EmotionalReversal,
    MisalignedTension,
    Release
};

struct SceneSegment {
    double startTime{0.0};
    double endTime{0.0};
    SceneType type{SceneType::Unknown};
    double avgPace{0.0};
    double avgSound{0.0};
    double avgArousal{0.0};
    double avgArousalSlope{0.0};        // A'(t)
    double avgPaceSoundDivergence{0.0}; // D_PA(t)
    int colorState{0};                  // Dominant color mode
};

struct SequenceSegment {
    double startTime{0.0};
    double endTime{0.0};
    std::vector<int> sceneIndices;      // Indices into SceneSegment array
};

struct ActSegment {
    double startTime{0.0};
    double endTime{0.0};
    int actNumber{1};                   // 1, 2, 3, ...
    std::vector<int> sequenceIndices;   // Indices into SequenceSegment array
};

enum class CurveKind {
    Pace,
    Sound,
    Color,
    Info,
    Arousal,
    FaceAffect,
    FacePresence
};

struct CurveSeries {
    CurveKind kind{CurveKind::Pace};
    std::vector<double> times;
    std::vector<double> values;        // Final curve values in (0,1) after robust_sigmoid01
    std::vector<double> raw_score;     // PC1 raw score in Z-space before sigmoid (for λ, for Arousal PCA)
};

enum class OverlayValueType {
    Float32,
    Int32
};

// Metadata Overlays: displayed as optional tracks (arcscope.md 5.6 / 10.7).
// Note: some Diagnostics rules (arcscope.md 6.x) intentionally consume derived overlays (e.g. AlignNorm/Tension).
struct MetadataOverlaySeries {
    std::string type;                  // e.g. "color_warmth", "bpm", "hue_hist12"
    double fps{1.0};                   // fixed to 1.0 for unified timeline
    int channels{1};                   // 1 for scalar tracks; >1 for multi-channel (e.g. hue hist)
    OverlayValueType valueType{OverlayValueType::Float32};
    std::vector<float> valuesF32;      // length = lengthSeconds * channels
    std::vector<int32_t> valuesI32;    // length = lengthSeconds * channels
};

struct IssueSegment {
    std::string type;          // e.g. "OvercutFlat", "AudioVisualMisalign"
    double startTime{0.0};
    double endTime{0.0};
    double severity{0.0};     // [0,1]
    std::string explanation;
};

struct FilmMetrics {
    std::unordered_map<std::string, double> scalars;
};

struct FilmAnalysisResult {
    std::string filmId;                   // logical_id (e.g., "blade-runner-1982")
    std::string filmTitle;                // Display title
    std::string filePath;                 // Actual file path (e.g., "/Users/foo/Videos/blade_runner.mp4")
    std::string imdbId;                   // IMDb ID (e.g., "tt0083658"), optional
    std::string source;                   // Data source: "local" / "imdb" / "omdb" / "tmdb"
    int year{0};
    std::string director;
    double durationSeconds{0.0};
    std::vector<CurveSeries> curves;
    std::vector<MetadataOverlaySeries> overlays;
    std::vector<IssueSegment> issues;
    FilmMetrics metrics;

    // Multi-scale structure
    std::vector<ShotSegment> shots;
    std::vector<SceneSegment> scenes;
    std::vector<SequenceSegment> sequences;
    std::vector<ActSegment> acts;

    // Face tracking info
    std::vector<FaceTrack> faceTracks;
    int dominantTrackId{-1};              // -1 = no reliable dominant face
    bool hasDominantFace{false};
};

// Analysis option flags (arcscope.md 5.6).
enum FilmAnalysisOptions : uint32_t {
    // Scene/Sequence segmentation runs change-point detection on curves reparameterized by emotion progress u(t) (arcscope.md 4.6.5).
    // Boundaries are mapped back to seconds for storage/UI.
    StructureUseEmotionProgressReparam = 1u << 0,
};

struct FilmAnalysisRequest {
    std::string filmId;                   // logical_id (e.g., "blade-runner-1982")
    std::string filmTitle;                // Display title
    std::string filePath;                 // Actual file path (optional for pure-data imports)
    std::string imdbId;                   // IMDb ID (optional)
    std::string source;                   // Data source: "local" / "imdb" / "omdb"
    int year{0};
    std::string director;
    double durationSeconds{0.0};
    std::vector<FeatureSample> samples;
    std::vector<FaceTrack> tracks;        // Face tracks from Bridge layer
    // Optional: exact cut timestamps in seconds (monotonic, unique, within [0,duration]).
    // When present, ShotSegments are built from these boundaries (doc 4.1.2 / 5.5).
    std::vector<double> cutTimes;

    uint32_t options{0};                  // FilmAnalysisOptions bitset
};

// ========== Internal Analysis Context ==========
// Used internally by FilmEngine during analysis
struct AnalysisContext {
    std::string filmPath;
    double runtime_sec{0.0};

    // unified 1Hz timeline centers: t_k = k + 0.5, k=0..N-1
    std::vector<double> t_center;

    std::vector<ShotSegment> shots;
    std::vector<AudioFeatures> audioFeatures;
    std::vector<ColorFeatures> colorFeatures;
    std::vector<InfoFeatures>  infoFeatures;
    std::vector<FaceFeatures>  faceFeatures;

    // Face analysis metadata (from FaceAnalyzer)
    int dominantTrackId{-1};          // -1 if no reliable dominant face
    bool faceEnabled{false};          // dominantTrackId >= 0

    // optional: if you have true 1Hz motion from optical flow aggregation, prefer this
    std::vector<double> motion_1hz; // length N, M(k)
    std::vector<double> camera_motion_1hz; // length N, camera motion magnitude (Bridge observation)
    std::vector<double> object_motion_1hz; // length N, object/residual motion magnitude (Bridge observation)
    std::vector<int32_t> camera_motion_type_1hz; // length N, discrete type (Bridge observation)
};

}  // namespace arcscope
