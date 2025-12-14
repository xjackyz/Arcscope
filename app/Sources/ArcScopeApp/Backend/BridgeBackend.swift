import CAArcScopeBridge
import Foundation

struct BridgeBackend: ArcScopeBackend {
    func loadFilm(databaseURL: URL?, filmId: String?) async throws -> FilmAnalysis {
        let dbURL = try resolveDatabaseURL(from: databaseURL)
        var detail = ArcScopeFilmDetail()
        let status = SecurityScopedAccess.withAccess(dbURL) {
            dbURL.path.withCString { dbCString in
                filmId?.withCString { filmCString in
                    arcscope_load_film_detail(dbCString, filmCString, &detail)
                } ?? arcscope_load_film_detail(dbCString, nil, &detail)
            }
        }
        guard status == 0 else {
            arcscope_bridge_free_detail(&detail)
            throw ArcScopeBackendError.decodeFailed
        }
        defer { arcscope_bridge_free_detail(&detail) }
        return try FilmAnalysis(detail: detail)
    }

    // MARK: - Helpers

    private func resolveDatabaseURL(from input: URL?) throws -> URL {
        if let input {
            return input
        }
        let defaultPath = AppPaths.defaultDatabaseURL()
        guard FileManager.default.fileExists(atPath: defaultPath.path) else {
            throw ArcScopeBackendError.fileMissing
        }
        return defaultPath
    }
}

private extension FilmAnalysis {
    init(detail: ArcScopeFilmDetail) throws {
        guard let filmIdC = detail.film_id,
              let titleC = detail.title else {
            throw ArcScopeBackendError.decodeFailed
        }
        let filePath = detail.file_path.map { String(cString: $0) }
        let director = detail.director.map { String(cString: $0) } ?? ""
        let curvesBuffer = detail.curves
        let curveModels: [FilmAnalysis.Curve] = (0..<detail.curve_count).compactMap { curveIndex in
            guard let buffer = curvesBuffer else { return nil }
            let curve = buffer[curveIndex]
            guard let kindCString = curve.kind,
                  let kind = CurveKind(rawValue: String(cString: kindCString)),
                  let times = curve.times,
                  let values = curve.values else {
                return nil
            }
            let samples: [FilmAnalysis.Sample] = (0..<curve.count).map { idx in
                FilmAnalysis.Sample(time: Double(times[idx]), value: Double(values[idx]))
            }
            return FilmAnalysis.Curve(kind: kind, samples: samples)
        }

        let overlaysBuffer = detail.overlays
        let overlayModels: [FilmAnalysis.Overlay] = (0..<detail.overlay_count).compactMap { overlayIndex in
            guard let buffer = overlaysBuffer else { return nil }
            let ov = buffer[overlayIndex]
            guard let kindCString = ov.kind,
                  let times = ov.times else {
                return nil
            }

            let kind = String(cString: kindCString)
            let count = Int(ov.count)
            let channels = Int(ov.channels)
            guard count > 0, channels > 0 else { return nil }

            let tvals: [Double] = (0..<count).map { Double(times[$0]) }

            if ov.value_type == ARCSCOPE_OVERLAY_I32 {
                guard let values = ov.values_i32 else { return nil }
                if channels == 1 {
                    let samples: [FilmAnalysis.Sample] = (0..<count).map { idx in
                        FilmAnalysis.Sample(time: tvals[idx], value: Double(values[idx]))
                    }
                    return FilmAnalysis.Overlay(kind: kind, valueType: .i32, channels: 1, samples: samples, multiValues: [])
                }
                var multi = Array(repeating: Array(repeating: 0.0, count: count), count: channels)
                for c in 0..<channels {
                    for i in 0..<count {
                        multi[c][i] = Double(values[i * channels + c])
                    }
                }
                return FilmAnalysis.Overlay(kind: kind, valueType: .i32, channels: channels, samples: [], multiValues: multi)
            } else {
                guard let values = ov.values_f32 else { return nil }
                if channels == 1 {
                    let samples: [FilmAnalysis.Sample] = (0..<count).map { idx in
                        FilmAnalysis.Sample(time: tvals[idx], value: Double(values[idx]))
                    }
                    return FilmAnalysis.Overlay(kind: kind, valueType: .f32, channels: 1, samples: samples, multiValues: [])
                }
                var multi = Array(repeating: Array(repeating: 0.0, count: count), count: channels)
                for c in 0..<channels {
                    for i in 0..<count {
                        multi[c][i] = Double(values[i * channels + c])
                    }
                }
                return FilmAnalysis.Overlay(kind: kind, valueType: .f32, channels: channels, samples: [], multiValues: multi)
            }
        }

        let issuesBuffer = detail.issues
        let issueModels: [FilmAnalysis.Issue] = (0..<detail.issue_count).compactMap { index in
            guard let buffer = issuesBuffer else { return nil }
            let issue = buffer[index]
            return FilmAnalysis.Issue(
                type: issue.type.map { String(cString: $0) } ?? "",
                start: issue.start_time,
                end: issue.end_time,
                severity: issue.severity,
                explanation: issue.explanation.map { String(cString: $0) } ?? ""
            )
        }
        let scenesBuffer = detail.scenes
        let sceneModels: [FilmAnalysis.SceneSegment] = (0..<detail.scene_count).compactMap { index in
            guard let buffer = scenesBuffer else { return nil }
            let scene = buffer[index]
            let label = scene.label.map { String(cString: $0) } ?? "Scene \(index + 1)"
            return FilmAnalysis.SceneSegment(label: label, start: scene.start_time, end: scene.end_time)
        }

        let shotsBuffer = detail.shots
        let shotModels: [FilmAnalysis.ShotSegment] = (0..<detail.shot_count).compactMap { index in
            guard let buffer = shotsBuffer else { return nil }
            let shot = buffer[index]
            return FilmAnalysis.ShotSegment(start: shot.start_time,
                                            end: shot.end_time,
                                            duration: shot.duration,
                                            avgMotion: shot.avg_motion)
        }

        let sequencesBuffer = detail.sequences
        let sequenceModels: [FilmAnalysis.SequenceSegment] = (0..<detail.sequence_count).compactMap { index in
            guard let buffer = sequencesBuffer else { return nil }
            let seg = buffer[index]
            return FilmAnalysis.SequenceSegment(start: seg.start_time, end: seg.end_time)
        }

        let actsBuffer = detail.acts
        let actModels: [FilmAnalysis.ActSegment] = (0..<detail.act_count).compactMap { index in
            guard let buffer = actsBuffer else { return nil }
            let act = buffer[index]
            return FilmAnalysis.ActSegment(number: Int(act.act_number),
                                           start: act.start_time,
                                           end: act.end_time)
        }
        self.init(
            filmId: String(cString: filmIdC),
            title: String(cString: titleC),
            year: Int(detail.year),
            director: director,
            filePath: filePath,
            duration: detail.duration_seconds,
            curves: curveModels,
            overlays: overlayModels,
            issues: issueModels,
            shots: shotModels,
            scenes: sceneModels,
            sequences: sequenceModels,
            acts: actModels
        )
    }
}
