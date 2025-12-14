#import "arcscope_bridge.h"

#import "../VideoFeatureExtractor.h"
#import "../FaceTrackingEngine.h"

#import <Foundation/Foundation.h>

#include <algorithm>
#include <cctype>
#include <cstring>
#include <map>
#include <memory>
#include <stdexcept>
#include <string>
#include <vector>

#include <sqlite3.h>

#include "arcscope/FilmEngine.h"
#include "arcscope/Utility.h"

struct ArcScopeBridgeAnalyzer {
    std::string dbPath;
    std::unique_ptr<arcscope::FilmEngine> engine;
};

namespace {

char* copy_string(const std::string& value) {
    char* buffer = new char[value.size() + 1];
    std::memcpy(buffer, value.c_str(), value.size() + 1);
    return buffer;
}

std::string sanitize_identifier(const std::string& raw) {
    std::string result;
    result.reserve(raw.size());
    for (char ch : raw) {
        if (std::isalnum(static_cast<unsigned char>(ch))) {
            result.push_back(static_cast<char>(std::tolower(ch)));
        } else if (!result.empty() && result.back() != '-') {
            result.push_back('-');
        }
    }
    if (!result.empty() && result.back() == '-') {
        result.pop_back();
    }
    if (result.empty()) {
        return "film";
    }
    return result;
}

std::string column_text(sqlite3_stmt* stmt, int idx) {
    const unsigned char* ptr = sqlite3_column_text(stmt, idx);
    return ptr ? reinterpret_cast<const char*>(ptr) : std::string();
}

std::string derive_film_id(const char* explicit_id,
                           const char* title,
                           const char* input_path) {
    if (explicit_id && explicit_id[0] != '\0') {
        return sanitize_identifier(explicit_id);
    }
    if (title && title[0] != '\0') {
        return sanitize_identifier(title);
    }
    if (input_path && input_path[0] != '\0') {
        std::string path = input_path;
        const auto slash = path.find_last_of("/\\");
        if (slash != std::string::npos) {
            path = path.substr(slash + 1);
        }
        const auto dot = path.find_last_of('.');
        if (dot != std::string::npos) {
            path = path.substr(0, dot);
        }
        return sanitize_identifier(path);
    }
    return "film";
}

// Convert C bridge FaceTrack to C++ core FaceTrack
arcscope::FaceTrack ConvertTrack(const ArcScopeBridgeFaceTrack& bridgeTrack) {
    arcscope::FaceTrack track;
    track.trackId = bridgeTrack.track_id;
    track.appearanceDuration = bridgeTrack.appearance_duration;
    track.avgFaceArea = bridgeTrack.avg_face_area;
    track.avgConfidence = bridgeTrack.avg_confidence;
    track.score = bridgeTrack.score;

    if (bridgeTrack.point_count > 0) {
        track.timePoints.assign(bridgeTrack.time_points, bridgeTrack.time_points + bridgeTrack.point_count);
        track.faceAreas.assign(bridgeTrack.face_areas, bridgeTrack.face_areas + bridgeTrack.point_count);
        track.confidences.assign(bridgeTrack.confidences, bridgeTrack.confidences + bridgeTrack.point_count);
        track.valences.assign(bridgeTrack.valences, bridgeTrack.valences + bridgeTrack.point_count);
        track.arousals.assign(bridgeTrack.arousals, bridgeTrack.arousals + bridgeTrack.point_count);
        if (bridgeTrack.yaws) track.yaws.assign(bridgeTrack.yaws, bridgeTrack.yaws + bridgeTrack.point_count);
        if (bridgeTrack.rolls) track.rolls.assign(bridgeTrack.rolls, bridgeTrack.rolls + bridgeTrack.point_count);
        if (bridgeTrack.downcast_ratios) track.downcastRatios.assign(bridgeTrack.downcast_ratios, bridgeTrack.downcast_ratios + bridgeTrack.point_count);
    }

    return track;
}

// Convert Objective-C FaceTrack array to C bridge tracks
std::vector<ArcScopeBridgeFaceTrack> ConvertObjCTracksToC(NSArray* trackArray) {
    std::vector<ArcScopeBridgeFaceTrack> result;
    if (!trackArray) return result;

    result.reserve(trackArray.count);
    for (FaceTrack* objcTrack in trackArray) {
        ArcScopeBridgeFaceTrack bridgeTrack;
        bridgeTrack.track_id = objcTrack.trackId;
        bridgeTrack.appearance_duration = objcTrack.appearanceDuration;
        bridgeTrack.avg_face_area = objcTrack.avgFaceArea;
        bridgeTrack.avg_confidence = objcTrack.avgConfidence;
        bridgeTrack.score = objcTrack.score;
        bridgeTrack.point_count = objcTrack.timePoints.count;

        if (bridgeTrack.point_count > 0) {
            // Allocate arrays - these will need to be freed later
            bridgeTrack.time_points = new double[bridgeTrack.point_count];
            bridgeTrack.face_areas = new double[bridgeTrack.point_count];
            bridgeTrack.confidences = new double[bridgeTrack.point_count];
        bridgeTrack.valences = new double[bridgeTrack.point_count];
        bridgeTrack.arousals = new double[bridgeTrack.point_count];
        bridgeTrack.yaws = new double[bridgeTrack.point_count];
        bridgeTrack.rolls = new double[bridgeTrack.point_count];
        bridgeTrack.downcast_ratios = new double[bridgeTrack.point_count];

            for (size_t i = 0; i < bridgeTrack.point_count; ++i) {
                bridgeTrack.time_points[i] = [objcTrack.timePoints[i] doubleValue];
                bridgeTrack.face_areas[i] = [objcTrack.faceAreas[i] doubleValue];
                bridgeTrack.confidences[i] = [objcTrack.confidences[i] doubleValue];
                bridgeTrack.valences[i] = [objcTrack.valences[i] doubleValue];
                bridgeTrack.arousals[i] = [objcTrack.arousals[i] doubleValue];
                bridgeTrack.yaws[i] = (objcTrack.yaws.count > i) ? [objcTrack.yaws[i] doubleValue] : 0.0;
                bridgeTrack.rolls[i] = (objcTrack.rolls.count > i) ? [objcTrack.rolls[i] doubleValue] : 0.0;
                bridgeTrack.downcast_ratios[i] = (objcTrack.downcastRatios.count > i) ? [objcTrack.downcastRatios[i] doubleValue] : 0.5;
            }
        } else {
            bridgeTrack.time_points = nullptr;
            bridgeTrack.face_areas = nullptr;
            bridgeTrack.confidences = nullptr;
            bridgeTrack.valences = nullptr;
            bridgeTrack.arousals = nullptr;
            bridgeTrack.yaws = nullptr;
            bridgeTrack.rolls = nullptr;
            bridgeTrack.downcast_ratios = nullptr;
        }

        result.push_back(bridgeTrack);
    }
    return result;
}

arcscope::FilmAnalysisRequest BuildRequest(const std::string& filmId,
                                           const std::string& title,
                                           const std::string& filePath,
                                           const std::string& imdbId,
                                           const std::string& source,
                                           int year,
                                           const std::string& director,
                                           const ArcScopeBridgeSample* samples,
                                           size_t sample_count,
                                           const ArcScopeBridgeFaceTrack* tracks,
                                           size_t track_count,
                                           const std::vector<double>& cutTimesSec,
                                           uint32_t options) {
    arcscope::FilmAnalysisRequest request;
    request.filmId = filmId;
    request.filmTitle = title.empty() ? filmId : title;
    request.filePath = filePath;
    request.imdbId = imdbId;
    request.source = source.empty() ? "local" : source;
    request.year = year;
    request.director = director;
    request.options = options;

    // TIME INVARIANT CONTRACT ENFORCEMENT (§ I.R3):
    // - sample_count MUST equal duration_seconds (integer)
    // - samples[i].time_seconds MUST equal i + 0.5 (center-aligned 1Hz)
    // - Any deviation violates the "统一时间轴" contract
    for (size_t i = 0; i < sample_count; ++i) {
        const double expected_time = static_cast<double>(i) + 0.5;
        const double actual_time = samples[i].time_seconds;
        if (std::abs(actual_time - expected_time) > 1e-3) {
            throw std::runtime_error(
                "TIME INVARIANT VIOLATION: sample[" + std::to_string(i) +
                "].time_seconds = " + std::to_string(actual_time) +
                ", expected = " + std::to_string(expected_time) +
                " (Contract: times must be center-aligned at i+0.5 for 1Hz)"
            );
        }
    }
    request.durationSeconds = static_cast<double>(sample_count);
    request.samples.reserve(sample_count);
    request.cutTimes = cutTimesSec;

    // Bridge 只输出原始观测 SecondObs，不做融合
    // info_density, dialogue_density, arousal_proxy 这些融合后的结果必须在 C++ FilmEngine 中计算
    for (size_t i = 0; i < sample_count; ++i) {
        arcscope::FeatureSample sample;
        sample.timeSeconds = samples[i].time_seconds;

        // 原始观测字段（来自视频分析）
        sample.cutDensity = samples[i].cut_density;
        sample.motionAmplitude = samples[i].motion_amplitude;
        sample.cameraMotion = samples[i].camera_motion;
        sample.objectMotion = samples[i].object_motion;
        sample.cameraMotionType = samples[i].camera_motion_type;
        sample.audioRms = samples[i].audio_rms;
        sample.audioTransient = samples[i].audio_transient;
        sample.spectralBalance = samples[i].spectral_balance;

        // FFT-derived audio observations (Bridge AudioAnalysisEngine)
        sample.audioLoudness = samples[i].audio_loudness;
        sample.audioSpectralFlux = samples[i].audio_spectral_flux;
        sample.audioSpectralCentroid = samples[i].audio_spectral_centroid;
        sample.audioZeroCrossingRate = samples[i].audio_zero_crossing_rate;
        sample.audioBpm = samples[i].audio_bpm;
        sample.audioTempoConfidence = samples[i].audio_tempo_confidence;
        sample.audioRhythmEnergy = samples[i].audio_rhythm_energy;
        sample.audioEstimatedKey = samples[i].audio_estimated_key;
        sample.audioHarmonicTension = samples[i].audio_harmonic_tension;
        sample.audioDialogueClarity = samples[i].audio_dialogue_clarity;
        sample.audioSpeechProbability = samples[i].audio_speech_probability;
        sample.audioChroma.assign(samples[i].audio_chroma, samples[i].audio_chroma + 12);
        sample.audioStereoWidth = samples[i].audio_stereo_width;
        sample.audioReverbAmount = samples[i].audio_reverb_amount;
        sample.audioDialogueEnergy = samples[i].audio_dialogue_energy;
        sample.audioMusicEnergy = samples[i].audio_music_energy;
        sample.audioEffectsEnergy = samples[i].audio_effects_energy;
        sample.audioDialogueDominance = samples[i].audio_dialogue_dominance;
        sample.audioMusicDominance = samples[i].audio_music_dominance;
        sample.audioEffectsDominance = samples[i].audio_effects_dominance;

        // CAM16-UCS color observations (Bridge color science pipeline)
        sample.cam16J = samples[i].cam16_j;
        sample.cam16M = samples[i].cam16_m;
        sample.cam16H = samples[i].cam16_h;
        sample.cam16Jp = samples[i].cam16_jp;
        sample.cam16Ap = samples[i].cam16_ap;
        sample.cam16Bp = samples[i].cam16_bp;
        sample.cam16DeltaE = samples[i].cam16_delta_e;
        for (size_t k = 0; k < 12; ++k) {
            sample.cam16HueHist[k] = samples[i].cam16_hue_hist[k];
        }
        sample.cam16Harmony = samples[i].cam16_harmony;

        sample.brightness = samples[i].brightness;
        sample.saturation = samples[i].saturation;
        sample.warmth = samples[i].warmth;
        sample.grayness = samples[i].grayness;
        sample.avgHue = samples[i].avg_hue;

        // 字幕密度（原始 words/sec）
        sample.subtitleDensity = samples[i].subtitle_density;
        sample.expositionProbability = samples[i].exposition_probability;

        // NLP features (from SubtitleNLPProcessor)
        sample.nlpWordCount = samples[i].nlp_word_count;
        sample.nlpSentenceComplexity = samples[i].nlp_sentence_complexity;
        sample.nlpNewConceptRatio = samples[i].nlp_new_concept_ratio;
        sample.nlpEntityDensity = samples[i].nlp_entity_density;

        // ⚠️ WARNING: info_density is METADATA ONLY (see header § 67-84)
        // Core MUST NOT use this in PCA/InfoAnalyzer to avoid version drift
        // Core should derive InfoDensity from atomic features instead
        sample.infoDensity = samples[i].info_density;

        // 人脸观测（如果模型没接，这些字段应该标记为缺失或设为 0）
        sample.facePresence = samples[i].face_presence;
        sample.faceArousal = samples[i].face_arousal;
        sample.faceValence = samples[i].face_valence;
        sample.faceConfidence = samples[i].face_confidence;

        request.samples.push_back(sample);
    }

    // Convert face tracks
    request.tracks.reserve(track_count);
    for (size_t i = 0; i < track_count; ++i) {
        request.tracks.push_back(ConvertTrack(tracks[i]));
    }

    return request;
}

struct FilmRow {
    std::string id;           // Database primary key (INTEGER)
    std::string logical_id;   // Logical ID (TEXT, stable identifier)
    std::string title;
    std::string director;
    std::string file_path;
    std::string imdb_id;
    std::string source;
    int year{0};
    double runtime_sec{0.0};
};

bool fetch_latest_film(sqlite3* db, FilmRow& row) {
    const char* sql = "SELECT id, COALESCE(logical_id,''), COALESCE(title,''), COALESCE(director,''), "
                      "COALESCE(file_path,''), COALESCE(imdb_id,''), COALESCE(source,'local'), COALESCE(year,0), runtime_sec "
                      "FROM films ORDER BY updated_at DESC LIMIT 1";
    sqlite3_stmt* stmt = nullptr;
    if (sqlite3_prepare_v2(db, sql, -1, &stmt, nullptr) != SQLITE_OK) {
        return false;
    }
    const bool hasRow = sqlite3_step(stmt) == SQLITE_ROW;
    if (hasRow) {
        row.id = column_text(stmt, 0);
        row.logical_id = column_text(stmt, 1);
        row.title = column_text(stmt, 2);
        row.director = column_text(stmt, 3);
        row.file_path = column_text(stmt, 4);
        row.imdb_id = column_text(stmt, 5);
        row.source = column_text(stmt, 6);
        row.year = sqlite3_column_int(stmt, 7);
        row.runtime_sec = sqlite3_column_double(stmt, 8);
    }
    sqlite3_finalize(stmt);
    return hasRow;
}

bool fetch_film_by_id(sqlite3* db, const std::string& logical_id, FilmRow& row) {
    const char* sql = "SELECT id, COALESCE(logical_id,''), COALESCE(title,''), COALESCE(director,''), "
                      "COALESCE(file_path,''), COALESCE(imdb_id,''), COALESCE(source,'local'), COALESCE(year,0), runtime_sec "
                      "FROM films WHERE logical_id=?";
    sqlite3_stmt* stmt = nullptr;
    if (sqlite3_prepare_v2(db, sql, -1, &stmt, nullptr) != SQLITE_OK) {
        return false;
    }
    sqlite3_bind_text(stmt, 1, logical_id.c_str(), -1, SQLITE_TRANSIENT);
    const bool hasRow = sqlite3_step(stmt) == SQLITE_ROW;
    if (hasRow) {
        row.id = column_text(stmt, 0);
        row.logical_id = column_text(stmt, 1);
        row.title = column_text(stmt, 2);
        row.director = column_text(stmt, 3);
        row.file_path = column_text(stmt, 4);
        row.imdb_id = column_text(stmt, 5);
        row.source = column_text(stmt, 6);
        row.year = sqlite3_column_int(stmt, 7);
        row.runtime_sec = sqlite3_column_double(stmt, 8);
    }
    sqlite3_finalize(stmt);
    return hasRow;
}

void reset_detail(ArcScopeFilmDetail* detail) {
    if (!detail) {
        return;
    }
    detail->film_id = nullptr;
    detail->title = nullptr;
    detail->director = nullptr;
    detail->file_path = nullptr;
    detail->year = 0;
    detail->duration_seconds = 0.0;
    detail->curves = nullptr;
    detail->curve_count = 0;
    detail->overlays = nullptr;
    detail->overlay_count = 0;
    detail->issues = nullptr;
    detail->issue_count = 0;
    detail->shots = nullptr;
    detail->shot_count = 0;
    detail->scenes = nullptr;
    detail->scene_count = 0;
    detail->sequences = nullptr;
    detail->sequence_count = 0;
    detail->acts = nullptr;
    detail->act_count = 0;
}

ArcScopeFilmSummary make_summary(const FilmRow& row) {
    ArcScopeFilmSummary summary{};
    summary.film_id = copy_string(row.logical_id.empty() ? row.id : row.logical_id);  // Use logical_id
    summary.title = copy_string(row.title.empty() ? row.logical_id : row.title);
    summary.director = copy_string(row.director);
    summary.year = row.year;
    summary.duration_seconds = row.runtime_sec;
    return summary;
}

void fill_issue_data(sqlite3* db,
                     const std::string& filmId,
                     ArcScopeFilmDetail* detail) {
    // CONTRACT ENFORCEMENT: segments 表只存储诊断问题 (issues)
    // 结构分割 (scene/sequence/act) 必须存储在专用表:
    //   - scene_segments
    //   - sequence_segments
    //   - act_segments
    // 混用会导致 UI 把 scene 边界当 issue 绘制到 heatmap
    const char* sql = "SELECT type, start_sec, end_sec, severity, COALESCE(explanation,'') FROM segments WHERE film_id=? ORDER BY start_sec";
    sqlite3_stmt* stmt = nullptr;
    if (sqlite3_prepare_v2(db, sql, -1, &stmt, nullptr) != SQLITE_OK) {
        throw std::runtime_error("Failed to prepare segments query");
    }
    sqlite3_bind_text(stmt, 1, filmId.c_str(), -1, SQLITE_TRANSIENT);
    std::vector<ArcScopeBridgeIssue> issues;
    while (sqlite3_step(stmt) == SQLITE_ROW) {
        ArcScopeBridgeIssue issue;
        issue.type = copy_string(column_text(stmt, 0));
        issue.start_time = sqlite3_column_double(stmt, 1);
        issue.end_time = sqlite3_column_double(stmt, 2);
        issue.severity = sqlite3_column_double(stmt, 3);
        issue.explanation = copy_string(column_text(stmt, 4));
        issues.push_back(issue);
    }
    sqlite3_finalize(stmt);
    detail->issue_count = issues.size();
    detail->issues = new ArcScopeBridgeIssue[issues.size()];
    std::copy(issues.begin(), issues.end(), detail->issues);
}

void fill_scene_data(sqlite3* db,
                     const std::string& filmId,
                     ArcScopeFilmDetail* detail) {
    // CONTRACT: scene_segments 表存储结构分割（与 issues 分离）
    const char* sql = "SELECT start_sec, end_sec, COALESCE(scene_type,'') FROM scene_segments WHERE film_id=? ORDER BY start_sec";
    sqlite3_stmt* stmt = nullptr;
    if (sqlite3_prepare_v2(db, sql, -1, &stmt, nullptr) != SQLITE_OK) {
        throw std::runtime_error("Failed to prepare scene_segments query");
    }
    sqlite3_bind_text(stmt, 1, filmId.c_str(), -1, SQLITE_TRANSIENT);
    std::vector<ArcScopeSceneSegment> scenes;
    while (sqlite3_step(stmt) == SQLITE_ROW) {
        ArcScopeSceneSegment scene;
        scene.start_time = sqlite3_column_double(stmt, 0);
        scene.end_time = sqlite3_column_double(stmt, 1);
        scene.label = copy_string(column_text(stmt, 2));
        scenes.push_back(scene);
    }
    sqlite3_finalize(stmt);
    detail->scene_count = scenes.size();
    detail->scenes = scenes.empty() ? nullptr : new ArcScopeSceneSegment[scenes.size()];
    std::copy(scenes.begin(), scenes.end(), detail->scenes);
}

void fill_shot_data(sqlite3* db,
                    const std::string& filmId,
                    ArcScopeFilmDetail* detail) {
    const char* sql = "SELECT start_sec, end_sec, duration, COALESCE(avg_motion,0.0) FROM shot_segments WHERE film_id=? ORDER BY start_sec";
    sqlite3_stmt* stmt = nullptr;
    if (sqlite3_prepare_v2(db, sql, -1, &stmt, nullptr) != SQLITE_OK) {
        throw std::runtime_error("Failed to prepare shot_segments query");
    }
    sqlite3_bind_text(stmt, 1, filmId.c_str(), -1, SQLITE_TRANSIENT);
    std::vector<ArcScopeShotSegment> shots;
    while (sqlite3_step(stmt) == SQLITE_ROW) {
        ArcScopeShotSegment shot;
        shot.start_time = sqlite3_column_double(stmt, 0);
        shot.end_time = sqlite3_column_double(stmt, 1);
        shot.duration = sqlite3_column_double(stmt, 2);
        shot.avg_motion = sqlite3_column_double(stmt, 3);
        shots.push_back(shot);
    }
    sqlite3_finalize(stmt);
    detail->shot_count = shots.size();
    detail->shots = shots.empty() ? nullptr : new ArcScopeShotSegment[shots.size()];
    std::copy(shots.begin(), shots.end(), detail->shots);
}

void fill_sequence_data(sqlite3* db,
                        const std::string& filmId,
                        ArcScopeFilmDetail* detail) {
    const char* sql = "SELECT start_sec, end_sec FROM sequence_segments WHERE film_id=? ORDER BY start_sec";
    sqlite3_stmt* stmt = nullptr;
    if (sqlite3_prepare_v2(db, sql, -1, &stmt, nullptr) != SQLITE_OK) {
        throw std::runtime_error("Failed to prepare sequence_segments query");
    }
    sqlite3_bind_text(stmt, 1, filmId.c_str(), -1, SQLITE_TRANSIENT);
    std::vector<ArcScopeSequenceSegment> sequences;
    while (sqlite3_step(stmt) == SQLITE_ROW) {
        ArcScopeSequenceSegment seq;
        seq.start_time = sqlite3_column_double(stmt, 0);
        seq.end_time = sqlite3_column_double(stmt, 1);
        sequences.push_back(seq);
    }
    sqlite3_finalize(stmt);
    detail->sequence_count = sequences.size();
    detail->sequences = sequences.empty() ? nullptr : new ArcScopeSequenceSegment[sequences.size()];
    std::copy(sequences.begin(), sequences.end(), detail->sequences);
}

void fill_act_data(sqlite3* db,
                   const std::string& filmId,
                   ArcScopeFilmDetail* detail) {
    const char* sql = "SELECT act_number, start_sec, end_sec FROM act_segments WHERE film_id=? ORDER BY start_sec";
    sqlite3_stmt* stmt = nullptr;
    if (sqlite3_prepare_v2(db, sql, -1, &stmt, nullptr) != SQLITE_OK) {
        throw std::runtime_error("Failed to prepare act_segments query");
    }
    sqlite3_bind_text(stmt, 1, filmId.c_str(), -1, SQLITE_TRANSIENT);
    std::vector<ArcScopeActSegment> acts;
    while (sqlite3_step(stmt) == SQLITE_ROW) {
        ArcScopeActSegment act;
        act.act_number = sqlite3_column_int(stmt, 0);
        act.start_time = sqlite3_column_double(stmt, 1);
        act.end_time = sqlite3_column_double(stmt, 2);
        acts.push_back(act);
    }
    sqlite3_finalize(stmt);
    detail->act_count = acts.size();
    detail->acts = acts.empty() ? nullptr : new ArcScopeActSegment[acts.size()];
    std::copy(acts.begin(), acts.end(), detail->acts);
}

// decode_f32_blob: 从 SQLite blob 中解码 float32 数组
static std::vector<float> decode_f32_blob(sqlite3_stmt* stmt, int col_blob, int length) {
    const void* blob = sqlite3_column_blob(stmt, col_blob);
    const int bytes = sqlite3_column_bytes(stmt, col_blob);
    if (!blob || bytes < (int)(length * sizeof(float))) {
        return {};
    }
    std::vector<float> out(length);
    std::memcpy(out.data(), blob, length * sizeof(float));
    return out;
}

void fill_curve_data(sqlite3* db,
                     const std::string& filmId,
                     ArcScopeFilmDetail* detail) {
    const char* sql = "SELECT curve_type, fps, length, data_blob FROM curves WHERE film_id=? ORDER BY curve_type";
    sqlite3_stmt* stmt = nullptr;
    if (sqlite3_prepare_v2(db, sql, -1, &stmt, nullptr) != SQLITE_OK) {
        throw std::runtime_error("Failed to prepare curves query");
    }
    sqlite3_bind_text(stmt, 1, filmId.c_str(), -1, SQLITE_TRANSIENT);

    std::vector<ArcScopeBridgeCurve> curves;
    while (sqlite3_step(stmt) == SQLITE_ROW) {
        std::string kind = column_text(stmt, 0);
        const float fps = static_cast<float>(sqlite3_column_double(stmt, 1));
        const int length = sqlite3_column_int(stmt, 2);
        auto values = decode_f32_blob(stmt, 3, length);

        // TIME INVARIANT CONTRACT ENFORCEMENT (§ I.R3):
        // - fps MUST be exactly 1.0 (统一时间轴)
        // - times[i] MUST be i + 0.5 (center-aligned 1Hz)
        // - Non-1Hz curves are FORBIDDEN (would break UI timeline consistency)
        if (std::abs(fps - 1.0f) > 1e-6f) {
            // Skip curves with fps != 1.0 (violates contract)
            continue;
        }

        // Validate basic data integrity
        if ((int)values.size() != length || length <= 0) {
            continue;
        }

        ArcScopeBridgeCurve curve{};
        curve.kind = copy_string(kind);
        curve.count = static_cast<size_t>(length);
        curve.times = new float[curve.count];
        curve.values = new float[curve.count];

        // STRICT CONTRACT: times[i] = i + 0.5 (center-aligned, 1Hz)
        // This MUST match the time semantics in ArcScopeBridgeSample
        for (int i = 0; i < length; ++i) {
            curve.times[i] = static_cast<float>(i) + 0.5f;
            curve.values[i] = values[i];
        }
        curves.push_back(curve);
    }
    sqlite3_finalize(stmt);

    detail->curve_count = curves.size();
    detail->curves = curves.empty() ? nullptr : new ArcScopeBridgeCurve[curves.size()];
    for (size_t i = 0; i < curves.size(); ++i) {
        detail->curves[i] = curves[i];
    }
}

std::vector<int32_t> decode_i32_blob(sqlite3_stmt* stmt, int col, int count) {
    std::vector<int32_t> out;
    if (count <= 0) return out;
    const void* blob = sqlite3_column_blob(stmt, col);
    const int bytes = sqlite3_column_bytes(stmt, col);
    if (!blob || bytes < (int)(count * (int)sizeof(int32_t))) return out;
    out.resize((size_t)count);
    std::memcpy(out.data(), blob, (size_t)count * sizeof(int32_t));
    return out;
}

void fill_overlay_data(sqlite3* db,
                       const std::string& filmId,
                       ArcScopeFilmDetail* detail) {
    const char* sql =
        "SELECT overlay_type, fps, length, channels, value_type, data_blob "
        "FROM metadata_overlays WHERE film_id=? ORDER BY overlay_type";
    sqlite3_stmt* stmt = nullptr;
    if (sqlite3_prepare_v2(db, sql, -1, &stmt, nullptr) != SQLITE_OK) {
        return;
    }
    sqlite3_bind_text(stmt, 1, filmId.c_str(), -1, SQLITE_TRANSIENT);

    std::vector<ArcScopeBridgeOverlay> overlays;
    while (sqlite3_step(stmt) == SQLITE_ROW) {
        std::string kind = column_text(stmt, 0);
        const float fps = static_cast<float>(sqlite3_column_double(stmt, 1));
        const int length = sqlite3_column_int(stmt, 2);
        const int channels = sqlite3_column_int(stmt, 3);
        const std::string valueType = column_text(stmt, 4);

        if (std::abs(fps - 1.0f) > 1e-6f) continue;
        if (length <= 0 || channels <= 0) continue;

        ArcScopeBridgeOverlay ov{};
        ov.kind = copy_string(kind);
        ov.count = (size_t)length;
        ov.channels = channels;
        ov.times = new float[ov.count];
        for (int i = 0; i < length; ++i) ov.times[i] = static_cast<float>(i) + 0.5f;

        const int total = length * channels;
        if (valueType == "i32") {
            ov.value_type = ARCSCOPE_OVERLAY_I32;
            auto values = decode_i32_blob(stmt, 5, total);
            if ((int)values.size() != total) {
                delete[] ov.kind;
                delete[] ov.times;
                continue;
            }
            ov.values_i32 = new int32_t[total];
            ov.values_f32 = nullptr;
            std::memcpy(ov.values_i32, values.data(), (size_t)total * sizeof(int32_t));
        } else {
            ov.value_type = ARCSCOPE_OVERLAY_F32;
            auto values = decode_f32_blob(stmt, 5, total);
            if ((int)values.size() != total) {
                delete[] ov.kind;
                delete[] ov.times;
                continue;
            }
            ov.values_f32 = new float[total];
            ov.values_i32 = nullptr;
            std::memcpy(ov.values_f32, values.data(), (size_t)total * sizeof(float));
        }

        overlays.push_back(ov);
    }
    sqlite3_finalize(stmt);

    detail->overlay_count = overlays.size();
    detail->overlays = overlays.empty() ? nullptr : new ArcScopeBridgeOverlay[overlays.size()];
    for (size_t i = 0; i < overlays.size(); ++i) {
        detail->overlays[i] = overlays[i];
    }
}

}  // namespace

int arcscope_list_films(const char* db_path, ArcScopeFilmList* out_list) {
    if (!db_path || !out_list) {
        return -1;
    }
    sqlite3* db = nullptr;
    if (sqlite3_open(db_path, &db) != SQLITE_OK) {
        return -2;
    }
    const char* sql = "SELECT id, COALESCE(logical_id,''), COALESCE(title,''), COALESCE(director,''), "
                      "COALESCE(imdb_id,''), COALESCE(source,'local'), COALESCE(year,0), runtime_sec "
                      "FROM films ORDER BY updated_at DESC";
    sqlite3_stmt* stmt = nullptr;
    if (sqlite3_prepare_v2(db, sql, -1, &stmt, nullptr) != SQLITE_OK) {
        sqlite3_close(db);
        return -3;
    }
    std::vector<ArcScopeFilmSummary> summaries;
    while (sqlite3_step(stmt) == SQLITE_ROW) {
        FilmRow row;
        row.id = column_text(stmt, 0);
        row.logical_id = column_text(stmt, 1);
        row.title = column_text(stmt, 2);
        row.director = column_text(stmt, 3);
        row.imdb_id = column_text(stmt, 4);
        row.source = column_text(stmt, 5);
        row.year = sqlite3_column_int(stmt, 6);
        row.runtime_sec = sqlite3_column_double(stmt, 7);
        summaries.push_back(make_summary(row));
    }
    sqlite3_finalize(stmt);
    sqlite3_close(db);
    out_list->count = summaries.size();
    out_list->items = new ArcScopeFilmSummary[summaries.size()];
    for (size_t i = 0; i < summaries.size(); ++i) {
        out_list->items[i] = summaries[i];
    }
    return 0;
}

int arcscope_list_reference_films(const char* db_path, ArcScopeFilmList* out_list) {
    if (!db_path || !out_list) {
        return -1;
    }
    sqlite3* db = nullptr;
    if (sqlite3_open(db_path, &db) != SQLITE_OK) {
        return -2;
    }
    const char* sql =
        "SELECT f.id, COALESCE(f.logical_id,''), COALESCE(f.title,''), COALESCE(f.director,''), "
        "COALESCE(f.imdb_id,''), COALESCE(f.source,'local'), COALESCE(f.year,0), f.runtime_sec "
        "FROM films f "
        "JOIN references_library r ON r.film_id = f.id "
        "ORDER BY f.updated_at DESC";
    sqlite3_stmt* stmt = nullptr;
    if (sqlite3_prepare_v2(db, sql, -1, &stmt, nullptr) != SQLITE_OK) {
        sqlite3_close(db);
        return -3;
    }
    std::vector<ArcScopeFilmSummary> summaries;
    while (sqlite3_step(stmt) == SQLITE_ROW) {
        FilmRow row;
        row.id = column_text(stmt, 0);
        row.logical_id = column_text(stmt, 1);
        row.title = column_text(stmt, 2);
        row.director = column_text(stmt, 3);
        row.imdb_id = column_text(stmt, 4);
        row.source = column_text(stmt, 5);
        row.year = sqlite3_column_int(stmt, 6);
        row.runtime_sec = sqlite3_column_double(stmt, 7);
        summaries.push_back(make_summary(row));
    }
    sqlite3_finalize(stmt);
    sqlite3_close(db);
    out_list->count = summaries.size();
    out_list->items = new ArcScopeFilmSummary[summaries.size()];
    for (size_t i = 0; i < summaries.size(); ++i) {
        out_list->items[i] = summaries[i];
    }
    return 0;
}

int arcscope_set_reference_film(const char* db_path,
                                const char* film_id,
                                int is_reference,
                                const char* note) {
    if (!db_path || !film_id) {
        return -1;
    }
    sqlite3* db = nullptr;
    if (sqlite3_open(db_path, &db) != SQLITE_OK) {
        return -2;
    }

    // Resolve logical_id -> films.id
    int64_t filmPk = -1;
    {
        const char* sql = "SELECT id FROM films WHERE logical_id = ? LIMIT 1";
        sqlite3_stmt* stmt = nullptr;
        if (sqlite3_prepare_v2(db, sql, -1, &stmt, nullptr) != SQLITE_OK) {
            sqlite3_close(db);
            return -3;
        }
        sqlite3_bind_text(stmt, 1, film_id, -1, SQLITE_TRANSIENT);
        if (sqlite3_step(stmt) == SQLITE_ROW) {
            filmPk = sqlite3_column_int64(stmt, 0);
        }
        sqlite3_finalize(stmt);
    }
    if (filmPk < 0) {
        sqlite3_close(db);
        return -4;
    }

    int rc = SQLITE_OK;
    if (is_reference) {
        const char* sql =
            "INSERT INTO references_library (film_id, note) VALUES (?, ?) "
            "ON CONFLICT(film_id) DO UPDATE SET note = excluded.note";
        sqlite3_stmt* stmt = nullptr;
        if (sqlite3_prepare_v2(db, sql, -1, &stmt, nullptr) != SQLITE_OK) {
            sqlite3_close(db);
            return -5;
        }
        sqlite3_bind_int64(stmt, 1, filmPk);
        if (note) sqlite3_bind_text(stmt, 2, note, -1, SQLITE_TRANSIENT);
        else sqlite3_bind_null(stmt, 2);
        rc = sqlite3_step(stmt);
        sqlite3_finalize(stmt);
    } else {
        const char* sql = "DELETE FROM references_library WHERE film_id = ?";
        sqlite3_stmt* stmt = nullptr;
        if (sqlite3_prepare_v2(db, sql, -1, &stmt, nullptr) != SQLITE_OK) {
            sqlite3_close(db);
            return -5;
        }
        sqlite3_bind_int64(stmt, 1, filmPk);
        rc = sqlite3_step(stmt);
        sqlite3_finalize(stmt);
    }

    sqlite3_close(db);
    return (rc == SQLITE_DONE) ? 0 : -6;
}

void arcscope_bridge_free_film_list(ArcScopeFilmList* list) {
    if (!list) {
        return;
    }
    for (size_t i = 0; i < list->count; ++i) {
        delete[] list->items[i].film_id;
        delete[] list->items[i].title;
        delete[] list->items[i].director;
    }
    delete[] list->items;
    list->items = nullptr;
    list->count = 0;
}

ArcScopeBridgeAnalyzer* arcscope_bridge_create(const char* database_path) {
    if (!database_path) {
        return nullptr;
    }
    try {
        auto analyzer = new ArcScopeBridgeAnalyzer;
        analyzer->dbPath = database_path;
        analyzer->engine = std::make_unique<arcscope::FilmEngine>(analyzer->dbPath);
        return analyzer;
    } catch (...) {
        return nullptr;
    }
}

void arcscope_bridge_destroy(ArcScopeBridgeAnalyzer* analyzer) {
    delete analyzer;
}

int arcscope_bridge_analyze(ArcScopeBridgeAnalyzer* analyzer,
                            const char* film_id,
                            const char* film_title,
                            const ArcScopeBridgeSample* samples,
                            size_t sample_count,
                            const ArcScopeBridgeFaceTrack* tracks,
                            size_t track_count) {
    return arcscope_bridge_analyze_with_options(analyzer,
                                               film_id,
                                               film_title,
                                               samples,
                                               sample_count,
                                               tracks,
                                               track_count,
                                               0);
}

int arcscope_bridge_analyze_with_options(ArcScopeBridgeAnalyzer* analyzer,
                                        const char* film_id,
                                        const char* film_title,
                                        const ArcScopeBridgeSample* samples,
                                        size_t sample_count,
                                        const ArcScopeBridgeFaceTrack* tracks,
                                        size_t track_count,
                                        uint32_t options) {
    if (!analyzer || !samples || sample_count == 0) {
        return -1;
    }
    const std::string resolvedId = derive_film_id(film_id, film_title, nullptr);
    try {
        const std::vector<double> cutTimesSec; // external callers: no exact cut timestamps
        auto request = BuildRequest(resolvedId,
                                    film_title ? film_title : resolvedId,
                                    "",        // filePath (empty for pure-data imports)
                                    "",        // imdbId (empty)
                                    "local",   // source
                                    0,         // year
                                    "",        // director
                                    samples,
                                    sample_count,
                                    tracks,
                                    track_count,
                                    cutTimesSec,
                                    options);
        analyzer->engine->analyze_and_store(request);
        return 0;
    } catch (...) {
        return -2;
    }
}

namespace arcscope {

// FilmAnalyzer (arcscope.md 7.2/9): Bridge-layer orchestrator for
// input_path -> raw observations -> FilmEngine -> SQLite.
class FilmAnalyzer {
public:
    FilmAnalyzer(const char* input_path,
                 const char* db_path,
                 const char* film_id,
                 const char* title,
                 int year,
                 const char* director,
                 uint32_t options)
        : inputPath_(input_path ? input_path : ""),
          dbPath_(db_path ? db_path : ""),
          explicitFilmId_(film_id ? film_id : ""),
          title_(title ? title : ""),
          year_(year),
          director_(director ? director : ""),
          options_(options) {}

    int run(ArcScopeProgressCallback progress_cb, void* user_data) {
        if (inputPath_.empty() || dbPath_.empty()) {
            return -1;
        }
        auto notify = [&](float value) {
            if (progress_cb) progress_cb(value, user_data);
        };

        notify(0.0f);
        @autoreleasepool {
            NSString* path = [NSString stringWithUTF8String:inputPath_.c_str()];
            NSURL* url = [NSURL fileURLWithPath:path];
            if (!url) {
                return -2;
            }

            notify(0.05f);

            std::vector<ArcScopeBridgeSample> samples;
            NSArray* trackArray = nil;
            std::vector<double> cutTimesSec;
            if (!ExtractSamplesFromVideo(url, samples, &trackArray, &cutTimesSec) || samples.empty()) {
                return -3;
            }

            notify(0.70f);

            std::vector<ArcScopeBridgeFaceTrack> bridgeTracks;
            if (trackArray && trackArray.count > 0) {
                bridgeTracks = ConvertObjCTracksToC(trackArray);
            }

            notify(0.75f);

            const std::string logicalId = derive_film_id(explicitFilmId_.empty() ? nullptr : explicitFilmId_.c_str(),
                                                        title_.empty() ? nullptr : title_.c_str(),
                                                        inputPath_.c_str());
            arcscope::FilmAnalysisRequest request = BuildRequest(
                logicalId,
                title_.empty() ? logicalId : title_,
                inputPath_,
                "",
                "local",
                year_,
                director_,
                samples.data(),
                samples.size(),
                bridgeTracks.data(),
                bridgeTracks.size(),
                cutTimesSec,
                options_);

            for (auto& track : bridgeTracks) {
                delete[] track.time_points;
                delete[] track.face_areas;
                delete[] track.confidences;
                delete[] track.valences;
                delete[] track.arousals;
                delete[] track.yaws;
                delete[] track.rolls;
                delete[] track.downcast_ratios;
            }

            notify(0.80f);

            arcscope::FilmEngine engine(dbPath_);
            engine.analyze_and_store(request);

            notify(0.95f);
            notify(1.0f);
            return 0;
        }
    }

private:
    std::string inputPath_;
    std::string dbPath_;
    std::string explicitFilmId_;
    std::string title_;
    int year_{0};
    std::string director_;
    uint32_t options_{0};
};

} // namespace arcscope

int arcscope_analyze_film(const char* input_path,
                          const char* db_path,
                          const char* film_id,
                          const char* title,
                          int year,
                          const char* director,
                          ArcScopeProgressCallback progress_cb,
                          void* user_data) {
    return arcscope_analyze_film_with_options(input_path,
                                             db_path,
                                             film_id,
                                             title,
                                             year,
                                             director,
                                             0,
                                             progress_cb,
                                             user_data);
}

int arcscope_analyze_film_with_options(const char* input_path,
                                      const char* db_path,
                                      const char* film_id,
                                      const char* title,
                                      int year,
                                      const char* director,
                                      uint32_t options,
                                      ArcScopeProgressCallback progress_cb,
                                      void* user_data) {
    arcscope::FilmAnalyzer analyzer(input_path, db_path, film_id, title, year, director, options);
    return analyzer.run(progress_cb, user_data);
}

int arcscope_load_film_detail(const char* db_path,
                              const char* film_id,
                              ArcScopeFilmDetail* out_detail) {
    if (!db_path || !out_detail) {
        return -1;
    }
    reset_detail(out_detail);
    sqlite3* db = nullptr;
    if (sqlite3_open(db_path, &db) != SQLITE_OK) {
        return -2;
    }
    FilmRow row;
    bool ok = false;
    if (film_id && film_id[0] != '\0') {
        ok = fetch_film_by_id(db, film_id, row);
    } else {
        ok = fetch_latest_film(db, row);
    }
    if (!ok) {
        sqlite3_close(db);
        return -3;
    }
    try {
        out_detail->film_id = copy_string(row.logical_id.empty() ? row.id : row.logical_id);
        out_detail->title = copy_string(row.title.empty() ? row.logical_id : row.title);
        out_detail->director = copy_string(row.director);
        out_detail->file_path = row.file_path.empty() ? nullptr : copy_string(row.file_path);
        out_detail->year = row.year;
        out_detail->duration_seconds = row.runtime_sec;

        // Bridge 只负责从数据库读取数据，不做结构推断
        // CONTRACT: issues 和 structure segments 已经在数据库层分离
        // Use database primary key (row.id) for querying child tables
        fill_curve_data(db, row.id, out_detail);
        fill_overlay_data(db, row.id, out_detail);
        fill_issue_data(db, row.id, out_detail);
        fill_shot_data(db, row.id, out_detail);
        fill_scene_data(db, row.id, out_detail);
        fill_sequence_data(db, row.id, out_detail);
        fill_act_data(db, row.id, out_detail);
    } catch (...) {
        sqlite3_close(db);
        arcscope_bridge_free_detail(out_detail);
        return -4;
    }
    sqlite3_close(db);
    return 0;
}

void arcscope_bridge_free_detail(ArcScopeFilmDetail* detail) {
    if (!detail) {
        return;
    }
    delete[] detail->film_id;
    delete[] detail->title;
    delete[] detail->director;
    delete[] detail->file_path;
    for (size_t i = 0; i < detail->curve_count; ++i) {
        delete[] detail->curves[i].kind;
        delete[] detail->curves[i].times;
        delete[] detail->curves[i].values;
    }
    delete[] detail->curves;
    for (size_t i = 0; i < detail->overlay_count; ++i) {
        delete[] detail->overlays[i].kind;
        delete[] detail->overlays[i].times;
        delete[] detail->overlays[i].values_f32;
        delete[] detail->overlays[i].values_i32;
    }
    delete[] detail->overlays;
    for (size_t i = 0; i < detail->issue_count; ++i) {
        delete[] detail->issues[i].type;
        delete[] detail->issues[i].explanation;
    }
    delete[] detail->issues;
    // Free shots (no string fields to delete)
    delete[] detail->shots;
    // Free scenes (has label string field)
    for (size_t i = 0; i < detail->scene_count; ++i) {
        delete[] detail->scenes[i].label;
    }
    delete[] detail->scenes;
    // Free sequences (no string fields to delete)
    delete[] detail->sequences;
    // Free acts (no string fields to delete)
    delete[] detail->acts;
    reset_detail(detail);
}

void arcscope_bridge_free_tracks(ArcScopeBridgeFaceTrack* tracks, size_t track_count) {
    if (!tracks) {
        return;
    }
    // Free each track's dynamically allocated arrays
    for (size_t i = 0; i < track_count; ++i) {
        delete[] tracks[i].time_points;
        delete[] tracks[i].face_areas;
        delete[] tracks[i].confidences;
        delete[] tracks[i].valences;
        delete[] tracks[i].arousals;
        delete[] tracks[i].yaws;
        delete[] tracks[i].rolls;
        delete[] tracks[i].downcast_ratios;
    }
    // Note: We do NOT delete[] tracks itself, as the caller owns that array
    // This function only frees the per-track arrays allocated by ConvertObjCTracksToC
}
