#include "arcscope/SQLiteStore.h"

#include <stdexcept>
#include <string>

#include "arcscope/Utility.h"

namespace arcscope {

namespace {

void exec_or_throw(sqlite3* db, const std::string& sql) {
    char* err = nullptr;
    if (sqlite3_exec(db, sql.c_str(), nullptr, nullptr, &err) != SQLITE_OK) {
        std::string message = err ? err : "Unknown sqlite error";
        sqlite3_free(err);
        throw std::runtime_error(message);
    }
}

sqlite3_stmt* prepare_or_throw(sqlite3* db, const std::string& sql) {
    sqlite3_stmt* stmt = nullptr;
    if (sqlite3_prepare_v2(db, sql.c_str(), -1, &stmt, nullptr) != SQLITE_OK) {
        throw std::runtime_error(sqlite3_errmsg(db));
    }
    return stmt;
}

bool column_exists(sqlite3* db, const std::string& table, const std::string& column) {
    std::string pragma = "PRAGMA table_info(" + table + ");";
    sqlite3_stmt* stmt = nullptr;
    if (sqlite3_prepare_v2(db, pragma.c_str(), -1, &stmt, nullptr) != SQLITE_OK) {
        return false;
    }
    bool exists = false;
    while (sqlite3_step(stmt) == SQLITE_ROW) {
        const unsigned char* name = sqlite3_column_text(stmt, 1);
        if (name && column == reinterpret_cast<const char*>(name)) {
            exists = true;
            break;
        }
    }
    sqlite3_finalize(stmt);
    return exists;
}

void ensure_film_columns(sqlite3* db) {
    if (!column_exists(db, "films", "title")) {
        exec_or_throw(db, "ALTER TABLE films ADD COLUMN title TEXT;");
    }
    if (!column_exists(db, "films", "director")) {
        exec_or_throw(db, "ALTER TABLE films ADD COLUMN director TEXT;");
    }
    if (!column_exists(db, "films", "year")) {
        exec_or_throw(db, "ALTER TABLE films ADD COLUMN year INTEGER;");
    }
}

}  // namespace

SQLiteStore::SQLiteStore(const std::string& path) {
    if (sqlite3_open(path.c_str(), &db_) != SQLITE_OK) {
        throw std::runtime_error("Failed to open SQLite database at " + path);
    }
}

SQLiteStore::~SQLiteStore() {
    if (db_) {
        sqlite3_close(db_);
    }
}

void SQLiteStore::initialize() {
    const char* schema = R"SQL(
        PRAGMA journal_mode=WAL;

        CREATE TABLE IF NOT EXISTS films (
            id INTEGER PRIMARY KEY AUTOINCREMENT,
            logical_id TEXT NOT NULL UNIQUE,
            title TEXT NOT NULL,
            year INTEGER,
            director TEXT,
            country TEXT,
            runtime_sec REAL,
            file_path TEXT,
            imdb_id TEXT,
            tmdb_id INTEGER,
            source TEXT,
            tags TEXT,
            my_rating REAL,
            ext_rating REAL,
            created_at TEXT DEFAULT CURRENT_TIMESTAMP,
            updated_at TEXT DEFAULT CURRENT_TIMESTAMP
        );

        CREATE INDEX IF NOT EXISTS idx_films_logical ON films(logical_id);
        CREATE INDEX IF NOT EXISTS idx_films_imdb ON films(imdb_id);
        CREATE INDEX IF NOT EXISTS idx_films_updated ON films(updated_at);

        -- Trigger to auto-update updated_at on any film modification
        CREATE TRIGGER IF NOT EXISTS update_films_timestamp
            AFTER UPDATE ON films FOR EACH ROW
        BEGIN
            UPDATE films SET updated_at = CURRENT_TIMESTAMP WHERE id = NEW.id;
        END;

	        CREATE TABLE IF NOT EXISTS curves (
	            id INTEGER PRIMARY KEY AUTOINCREMENT,
	            film_id INTEGER NOT NULL,
	            curve_type TEXT NOT NULL,
	            fps REAL NOT NULL,
	            length INTEGER NOT NULL,
	            data_blob BLOB NOT NULL,
	            FOREIGN KEY(film_id) REFERENCES films(id) ON DELETE CASCADE
	        );

	        CREATE INDEX IF NOT EXISTS idx_curves_film ON curves(film_id);

	        -- Metadata overlays: display-only tracks (never enter PCA/Diagnostics)
	        CREATE TABLE IF NOT EXISTS metadata_overlays (
	            id INTEGER PRIMARY KEY AUTOINCREMENT,
	            film_id INTEGER NOT NULL,
	            overlay_type TEXT NOT NULL,
	            fps REAL NOT NULL,
	            length INTEGER NOT NULL,
	            channels INTEGER NOT NULL DEFAULT 1,
	            value_type TEXT NOT NULL, -- "f32" or "i32"
	            data_blob BLOB NOT NULL,
	            FOREIGN KEY(film_id) REFERENCES films(id) ON DELETE CASCADE
	        );

	        CREATE INDEX IF NOT EXISTS idx_overlays_film ON metadata_overlays(film_id);

	        CREATE TABLE IF NOT EXISTS metrics (
	            id INTEGER PRIMARY KEY AUTOINCREMENT,
	            film_id INTEGER NOT NULL,
            name TEXT NOT NULL,
            value REAL NOT NULL,
            FOREIGN KEY(film_id) REFERENCES films(id) ON DELETE CASCADE
        );

        CREATE INDEX IF NOT EXISTS idx_metrics_film ON metrics(film_id);

        CREATE TABLE IF NOT EXISTS segments (
            id INTEGER PRIMARY KEY AUTOINCREMENT,
            film_id INTEGER NOT NULL,
            type TEXT NOT NULL,
            start_sec REAL NOT NULL,
            end_sec REAL NOT NULL,
            severity REAL NOT NULL,
            explanation TEXT,
            FOREIGN KEY(film_id) REFERENCES films(id) ON DELETE CASCADE
        );

        CREATE INDEX IF NOT EXISTS idx_segments_film ON segments(film_id);

        CREATE TABLE IF NOT EXISTS references_library (
            id INTEGER PRIMARY KEY AUTOINCREMENT,
            film_id INTEGER NOT NULL UNIQUE,
            note TEXT,
            FOREIGN KEY(film_id) REFERENCES films(id) ON DELETE CASCADE
        );

        CREATE TABLE IF NOT EXISTS shot_segments (
            id INTEGER PRIMARY KEY AUTOINCREMENT,
            film_id INTEGER NOT NULL,
            shot_number INTEGER NOT NULL,
            start_sec REAL NOT NULL,
            end_sec REAL NOT NULL,
            duration REAL NOT NULL,
            avg_motion REAL,
            FOREIGN KEY(film_id) REFERENCES films(id) ON DELETE CASCADE
        );

        CREATE INDEX IF NOT EXISTS idx_shot_segments_film ON shot_segments(film_id);

        CREATE TABLE IF NOT EXISTS scene_segments (
            id INTEGER PRIMARY KEY AUTOINCREMENT,
            film_id INTEGER NOT NULL,
            scene_number INTEGER NOT NULL,
            start_sec REAL NOT NULL,
            end_sec REAL NOT NULL,
            scene_type TEXT,
            avg_pace REAL,
            avg_sound REAL,
            avg_arousal REAL,
            avg_arousal_slope REAL,
            avg_pace_sound_divergence REAL,
            color_state INTEGER,
            FOREIGN KEY(film_id) REFERENCES films(id) ON DELETE CASCADE
        );

        CREATE TABLE IF NOT EXISTS sequence_segments (
            id INTEGER PRIMARY KEY AUTOINCREMENT,
            film_id INTEGER NOT NULL,
            sequence_number INTEGER NOT NULL,
            start_sec REAL NOT NULL,
            end_sec REAL NOT NULL,
            scene_indices TEXT,
            FOREIGN KEY(film_id) REFERENCES films(id) ON DELETE CASCADE
        );

        CREATE TABLE IF NOT EXISTS act_segments (
            id INTEGER PRIMARY KEY AUTOINCREMENT,
            film_id INTEGER NOT NULL,
            act_number INTEGER NOT NULL,
            start_sec REAL NOT NULL,
            end_sec REAL NOT NULL,
            sequence_indices TEXT,
            FOREIGN KEY(film_id) REFERENCES films(id) ON DELETE CASCADE
        );
    )SQL";

    exec_or_throw(db_, schema);
}

void SQLiteStore::save_analysis(const FilmAnalysisResult& result) {
    exec_or_throw(db_, "BEGIN TRANSACTION;");
    try {
        // Insert film record with proper UPSERT on logical_id
        int film_id = -1;
        {
            auto stmt = prepare_or_throw(db_,
                "INSERT INTO films (logical_id, title, year, director, runtime_sec, file_path, imdb_id, source, updated_at) "
                "VALUES (?,?,?,?,?,?,?,?, CURRENT_TIMESTAMP) "
                "ON CONFLICT(logical_id) DO UPDATE SET "
                "title=excluded.title, year=excluded.year, director=excluded.director, "
                "runtime_sec=excluded.runtime_sec, file_path=excluded.file_path, "
                "imdb_id=excluded.imdb_id, source=excluded.source, updated_at=CURRENT_TIMESTAMP "
                "RETURNING id;");

            sqlite3_bind_text(stmt, 1, result.filmId.c_str(), -1, SQLITE_TRANSIENT);        // logical_id
            sqlite3_bind_text(stmt, 2, result.filmTitle.c_str(), -1, SQLITE_TRANSIENT);     // title
            sqlite3_bind_int(stmt, 3, result.year);                                          // year
            sqlite3_bind_text(stmt, 4, result.director.c_str(), -1, SQLITE_TRANSIENT);      // director
            sqlite3_bind_double(stmt, 5, result.durationSeconds);                            // runtime_sec
            sqlite3_bind_text(stmt, 6, result.filePath.c_str(), -1, SQLITE_TRANSIENT);      // file_path

            // Handle optional imdb_id (bind NULL if empty)
            if (result.imdbId.empty()) {
                sqlite3_bind_null(stmt, 7);
            } else {
                sqlite3_bind_text(stmt, 7, result.imdbId.c_str(), -1, SQLITE_TRANSIENT);
            }

            // Handle optional source (default to "local" if empty)
            std::string source = result.source.empty() ? "local" : result.source;
            sqlite3_bind_text(stmt, 8, source.c_str(), -1, SQLITE_TRANSIENT);               // source

            if (sqlite3_step(stmt) == SQLITE_ROW) {
                film_id = sqlite3_column_int(stmt, 0);
            }
            sqlite3_finalize(stmt);
        }

        if (film_id < 0) {
            throw std::runtime_error("Failed to insert/get film_id");
        }

	        // Delete existing curves/metrics/segments for this film
	        {
	            auto stmt = prepare_or_throw(db_, "DELETE FROM curves WHERE film_id=?;");
	            sqlite3_bind_int(stmt, 1, film_id);
	            sqlite3_step(stmt);
	            sqlite3_finalize(stmt);
	        }
	        {
	            auto stmt = prepare_or_throw(db_, "DELETE FROM metadata_overlays WHERE film_id=?;");
	            sqlite3_bind_int(stmt, 1, film_id);
	            sqlite3_step(stmt);
	            sqlite3_finalize(stmt);
	        }
	        {
	            auto stmt = prepare_or_throw(db_, "DELETE FROM metrics WHERE film_id=?;");
	            sqlite3_bind_int(stmt, 1, film_id);
	            sqlite3_step(stmt);
            sqlite3_finalize(stmt);
        }
        {
            auto stmt = prepare_or_throw(db_, "DELETE FROM segments WHERE film_id=?;");
            sqlite3_bind_int(stmt, 1, film_id);
            sqlite3_step(stmt);
            sqlite3_finalize(stmt);
        }

	        // Insert curves as BLOBs
	        {
	            auto stmt = prepare_or_throw(db_,
	                "INSERT INTO curves (film_id, curve_type, fps, length, data_blob) VALUES (?,?,?,?,?);");

            for (const auto& curve : result.curves) {
                const std::string kindString = curve_kind_to_string(curve.kind);

                // Serialize values to BLOB (float array)
                std::vector<float> float_values(curve.values.begin(), curve.values.end());
                const void* blob_data = float_values.data();
                int blob_size = (int)(float_values.size() * sizeof(float));

                sqlite3_bind_int(stmt, 1, film_id);
                sqlite3_bind_text(stmt, 2, kindString.c_str(), -1, SQLITE_TRANSIENT);
                sqlite3_bind_double(stmt, 3, 1.0);  // fps = 1 Hz
                sqlite3_bind_int(stmt, 4, (int)curve.values.size());
                sqlite3_bind_blob(stmt, 5, blob_data, blob_size, SQLITE_TRANSIENT);

                if (sqlite3_step(stmt) != SQLITE_DONE) {
                    throw std::runtime_error(sqlite3_errmsg(db_));
                }
                sqlite3_reset(stmt);
            }
	            sqlite3_finalize(stmt);
	        }

	        // Insert metadata overlays as BLOBs
	        {
	            auto stmt = prepare_or_throw(db_,
	                "INSERT INTO metadata_overlays (film_id, overlay_type, fps, length, channels, value_type, data_blob) "
	                "VALUES (?,?,?,?,?,?,?);");

	            for (const auto& ov : result.overlays) {
	                sqlite3_bind_int(stmt, 1, film_id);
	                sqlite3_bind_text(stmt, 2, ov.type.c_str(), -1, SQLITE_TRANSIENT);
	                sqlite3_bind_double(stmt, 3, ov.fps);

	                const int length = (ov.channels > 0) ? (int)((ov.valueType == OverlayValueType::Float32 ? ov.valuesF32.size() : ov.valuesI32.size()) / (size_t)ov.channels)
	                                                   : 0;
	                sqlite3_bind_int(stmt, 4, length);
	                sqlite3_bind_int(stmt, 5, ov.channels);

	                if (ov.valueType == OverlayValueType::Float32) {
	                    sqlite3_bind_text(stmt, 6, "f32", -1, SQLITE_TRANSIENT);
	                    const void* blob_data = ov.valuesF32.data();
	                    int blob_size = (int)(ov.valuesF32.size() * sizeof(float));
	                    sqlite3_bind_blob(stmt, 7, blob_data, blob_size, SQLITE_TRANSIENT);
	                } else {
	                    sqlite3_bind_text(stmt, 6, "i32", -1, SQLITE_TRANSIENT);
	                    const void* blob_data = ov.valuesI32.data();
	                    int blob_size = (int)(ov.valuesI32.size() * sizeof(int32_t));
	                    sqlite3_bind_blob(stmt, 7, blob_data, blob_size, SQLITE_TRANSIENT);
	                }

	                if (sqlite3_step(stmt) != SQLITE_DONE) {
	                    throw std::runtime_error(sqlite3_errmsg(db_));
	                }
	                sqlite3_reset(stmt);
	            }
	            sqlite3_finalize(stmt);
	        }

	        // Insert metrics
	        {
	            auto stmt = prepare_or_throw(db_, "INSERT INTO metrics (film_id, name, value) VALUES (?,?,?);");
            for (const auto& [name, value] : result.metrics.scalars) {
                sqlite3_bind_int(stmt, 1, film_id);
                sqlite3_bind_text(stmt, 2, name.c_str(), -1, SQLITE_TRANSIENT);
                sqlite3_bind_double(stmt, 3, value);

                if (sqlite3_step(stmt) != SQLITE_DONE) {
                    throw std::runtime_error(sqlite3_errmsg(db_));
                }
                sqlite3_reset(stmt);
            }
            sqlite3_finalize(stmt);
        }

        // Insert issue segments
        {
            auto stmt = prepare_or_throw(db_,
                "INSERT INTO segments (film_id, type, start_sec, end_sec, severity, explanation) VALUES (?,?,?,?,?,?);");

            for (const auto& issue : result.issues) {
                sqlite3_bind_int(stmt, 1, film_id);
                sqlite3_bind_text(stmt, 2, issue.type.c_str(), -1, SQLITE_TRANSIENT);
                sqlite3_bind_double(stmt, 3, issue.startTime);
                sqlite3_bind_double(stmt, 4, issue.endTime);
                sqlite3_bind_double(stmt, 5, issue.severity);
                sqlite3_bind_text(stmt, 6, issue.explanation.c_str(), -1, SQLITE_TRANSIENT);

                if (sqlite3_step(stmt) != SQLITE_DONE) {
                    throw std::runtime_error(sqlite3_errmsg(db_));
                }
                sqlite3_reset(stmt);
            }
            sqlite3_finalize(stmt);
        }

        // Delete and insert shot segments
        {
            auto del_stmt = prepare_or_throw(db_, "DELETE FROM shot_segments WHERE film_id=?;");
            sqlite3_bind_int(del_stmt, 1, film_id);
            sqlite3_step(del_stmt);
            sqlite3_finalize(del_stmt);

            auto stmt = prepare_or_throw(db_,
                "INSERT INTO shot_segments (film_id, shot_number, start_sec, end_sec, duration, avg_motion) "
                "VALUES (?,?,?,?,?,?);");

            int shot_num = 0;
            for (const auto& shot : result.shots) {
                sqlite3_bind_int(stmt, 1, film_id);
                sqlite3_bind_int(stmt, 2, shot_num++);
                sqlite3_bind_double(stmt, 3, shot.startTime);
                sqlite3_bind_double(stmt, 4, shot.endTime);
                sqlite3_bind_double(stmt, 5, shot.duration);
                sqlite3_bind_double(stmt, 6, shot.avgMotion);

                if (sqlite3_step(stmt) != SQLITE_DONE) {
                    throw std::runtime_error(sqlite3_errmsg(db_));
                }
                sqlite3_reset(stmt);
            }
            sqlite3_finalize(stmt);
        }

        // Delete and insert scene segments
        {
            auto del_stmt = prepare_or_throw(db_, "DELETE FROM scene_segments WHERE film_id=?;");
            sqlite3_bind_int(del_stmt, 1, film_id);
            sqlite3_step(del_stmt);
            sqlite3_finalize(del_stmt);

            auto stmt = prepare_or_throw(db_,
                "INSERT INTO scene_segments (film_id, scene_number, start_sec, end_sec, scene_type, "
                "avg_pace, avg_sound, avg_arousal, avg_arousal_slope, avg_pace_sound_divergence, color_state) "
                "VALUES (?,?,?,?,?,?,?,?,?,?,?);");

            int scene_num = 0;
            for (const auto& scene : result.scenes) {
                sqlite3_bind_int(stmt, 1, film_id);
                sqlite3_bind_int(stmt, 2, scene_num++);
                sqlite3_bind_double(stmt, 3, scene.startTime);
                sqlite3_bind_double(stmt, 4, scene.endTime);

                // Convert SceneType enum to string
                const char* scene_type_str = "Unknown";
                switch(scene.type) {
                    case SceneType::LowEnergyBasin: scene_type_str = "LowEnergyBasin"; break;
                    case SceneType::BuildUp: scene_type_str = "BuildUp"; break;
                    case SceneType::HighEnergyPeak: scene_type_str = "HighEnergyPeak"; break;
                    case SceneType::EmotionalReversal: scene_type_str = "EmotionalReversal"; break;
                    case SceneType::MisalignedTension: scene_type_str = "MisalignedTension"; break;
                    case SceneType::Release: scene_type_str = "Release"; break;
                    default: break;
                }
                sqlite3_bind_text(stmt, 5, scene_type_str, -1, SQLITE_TRANSIENT);

                sqlite3_bind_double(stmt, 6, scene.avgPace);
                sqlite3_bind_double(stmt, 7, scene.avgSound);
                sqlite3_bind_double(stmt, 8, scene.avgArousal);
                sqlite3_bind_double(stmt, 9, scene.avgArousalSlope);
                sqlite3_bind_double(stmt, 10, scene.avgPaceSoundDivergence);
                sqlite3_bind_int(stmt, 11, scene.colorState);

                if (sqlite3_step(stmt) != SQLITE_DONE) {
                    throw std::runtime_error(sqlite3_errmsg(db_));
                }
                sqlite3_reset(stmt);
            }
            sqlite3_finalize(stmt);
        }

        // Delete and insert sequence segments
        {
            auto del_stmt = prepare_or_throw(db_, "DELETE FROM sequence_segments WHERE film_id=?;");
            sqlite3_bind_int(del_stmt, 1, film_id);
            sqlite3_step(del_stmt);
            sqlite3_finalize(del_stmt);

            auto stmt = prepare_or_throw(db_,
                "INSERT INTO sequence_segments (film_id, sequence_number, start_sec, end_sec, scene_indices) "
                "VALUES (?,?,?,?,?);");

            int seq_num = 0;
            for (const auto& seq : result.sequences) {
                sqlite3_bind_int(stmt, 1, film_id);
                sqlite3_bind_int(stmt, 2, seq_num++);
                sqlite3_bind_double(stmt, 3, seq.startTime);
                sqlite3_bind_double(stmt, 4, seq.endTime);

                // Serialize scene indices as comma-separated string
                std::string indices;
                for (size_t i = 0; i < seq.sceneIndices.size(); ++i) {
                    if (i > 0) indices += ",";
                    indices += std::to_string(seq.sceneIndices[i]);
                }
                sqlite3_bind_text(stmt, 5, indices.c_str(), -1, SQLITE_TRANSIENT);

                if (sqlite3_step(stmt) != SQLITE_DONE) {
                    throw std::runtime_error(sqlite3_errmsg(db_));
                }
                sqlite3_reset(stmt);
            }
            sqlite3_finalize(stmt);
        }

        // Delete and insert act segments
        {
            auto del_stmt = prepare_or_throw(db_, "DELETE FROM act_segments WHERE film_id=?;");
            sqlite3_bind_int(del_stmt, 1, film_id);
            sqlite3_step(del_stmt);
            sqlite3_finalize(del_stmt);

            auto stmt = prepare_or_throw(db_,
                "INSERT INTO act_segments (film_id, act_number, start_sec, end_sec, sequence_indices) "
                "VALUES (?,?,?,?,?);");

            for (const auto& act : result.acts) {
                sqlite3_bind_int(stmt, 1, film_id);
                sqlite3_bind_int(stmt, 2, act.actNumber);
                sqlite3_bind_double(stmt, 3, act.startTime);
                sqlite3_bind_double(stmt, 4, act.endTime);

                // Serialize sequence indices as comma-separated string
                std::string indices;
                for (size_t i = 0; i < act.sequenceIndices.size(); ++i) {
                    if (i > 0) indices += ",";
                    indices += std::to_string(act.sequenceIndices[i]);
                }
                sqlite3_bind_text(stmt, 5, indices.c_str(), -1, SQLITE_TRANSIENT);

                if (sqlite3_step(stmt) != SQLITE_DONE) {
                    throw std::runtime_error(sqlite3_errmsg(db_));
                }
                sqlite3_reset(stmt);
            }
            sqlite3_finalize(stmt);
        }

        exec_or_throw(db_, "COMMIT;");
    } catch (...) {
        exec_or_throw(db_, "ROLLBACK;");
        throw;
    }
}

}  // namespace arcscope
