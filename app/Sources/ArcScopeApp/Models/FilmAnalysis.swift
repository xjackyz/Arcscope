import Foundation
import SwiftUI

struct FilmAnalysis: Identifiable, Codable {
    var id: String { filmId }
    let filmId: String
    let title: String
    let year: Int
    let director: String
    let filePath: String?
    let duration: Double
    let curves: [Curve]
    let overlays: [Overlay]
    let issues: [Issue]
    let shots: [ShotSegment]
    let scenes: [SceneSegment]
    let sequences: [SequenceSegment]
    let acts: [ActSegment]

    struct Curve: Identifiable, Codable {
        let kind: CurveKind
        let samples: [Sample]

        var id: String { kind.rawValue }
        var color: Color {
            kind.color
        }
    }

    struct Sample: Codable {
        let time: Double
        let value: Double
    }

    enum OverlayValueType: String, Codable {
        case f32
        case i32
    }

    struct Overlay: Identifiable, Codable {
        let kind: String
        let valueType: OverlayValueType
        let channels: Int
        let samples: [Sample]             // used when channels==1
        let multiValues: [[Double]]       // used when channels>1, shape: [channels][count]

        var id: String { kind }
    }

    struct Issue: Identifiable, Codable {
        let type: String
        let start: Double
        let end: Double
        let severity: Double
        let explanation: String

        var id: String {
            "\(type)-\(start)-\(end)"
        }
    }

    struct SceneSegment: Identifiable, Codable {
        let label: String
        let start: Double
        let end: Double

        var id: String { "\(label)-\(start)" }
        var duration: Double { max(0, end - start) }
    }

    struct ShotSegment: Identifiable, Codable {
        let start: Double
        let end: Double
        let duration: Double
        let avgMotion: Double

        var id: String { "shot-\(start)" }
    }

    struct SequenceSegment: Identifiable, Codable {
        let start: Double
        let end: Double

        var id: String { "seq-\(start)" }
        var duration: Double { max(0, end - start) }
    }

    struct ActSegment: Identifiable, Codable {
        let number: Int
        let start: Double
        let end: Double

        var id: String { "act-\(number)-\(start)" }
        var duration: Double { max(0, end - start) }
    }
}

enum CurveKind: String, Codable, CaseIterable, Identifiable {
    case pace
    case sound
    case color
    case info
    case arousal
    case faceAffect = "face_affect"
    case facePresence = "face_presence"

    var id: String { rawValue }

    var label: String {
        switch self {
        case .pace: return "Pace"
        case .sound: return "Sound"
        case .color: return "Color"
        case .info: return "Info"
        case .arousal: return "Arousal"
        case .faceAffect: return "FaceAffect"
        case .facePresence: return "FacePresence"
        }
    }

    var color: Color {
        switch self {
        case .pace:
            return Color.orange
        case .sound:
            return Color.cyan
        case .color:
            return Color.pink
        case .info:
            return Color.purple
        case .arousal:
            return Color.red
        case .faceAffect:
            return Color.green
        case .facePresence:
            return Color.white.opacity(0.55)
        }
    }
}

extension FilmAnalysis {
    private enum CodingKeys: String, CodingKey {
        case filmId, title, year, director, filePath, duration, curves, overlays, issues, shots, scenes, sequences, acts
    }

    init(from decoder: Decoder) throws {
        let container = try decoder.container(keyedBy: CodingKeys.self)
        filmId = try container.decode(String.self, forKey: .filmId)
        title = try container.decode(String.self, forKey: .title)
        year = try container.decodeIfPresent(Int.self, forKey: .year) ?? 0
        director = try container.decodeIfPresent(String.self, forKey: .director) ?? ""
        filePath = try container.decodeIfPresent(String.self, forKey: .filePath)
        duration = try container.decode(Double.self, forKey: .duration)
        curves = try container.decode([Curve].self, forKey: .curves)
        overlays = try container.decodeIfPresent([Overlay].self, forKey: .overlays) ?? []
        issues = try container.decode([Issue].self, forKey: .issues)
        shots = try container.decodeIfPresent([ShotSegment].self, forKey: .shots) ?? []
        scenes = try container.decodeIfPresent([SceneSegment].self, forKey: .scenes) ?? []
        sequences = try container.decodeIfPresent([SequenceSegment].self, forKey: .sequences) ?? []
        acts = try container.decodeIfPresent([ActSegment].self, forKey: .acts) ?? []
    }

    func encode(to encoder: Encoder) throws {
        var container = encoder.container(keyedBy: CodingKeys.self)
        try container.encode(filmId, forKey: .filmId)
        try container.encode(title, forKey: .title)
        try container.encode(year, forKey: .year)
        try container.encode(director, forKey: .director)
        try container.encodeIfPresent(filePath, forKey: .filePath)
        try container.encode(duration, forKey: .duration)
        try container.encode(curves, forKey: .curves)
        try container.encode(overlays, forKey: .overlays)
        try container.encode(issues, forKey: .issues)
        try container.encode(shots, forKey: .shots)
        try container.encode(scenes, forKey: .scenes)
        try container.encode(sequences, forKey: .sequences)
        try container.encode(acts, forKey: .acts)
    }
}

extension FilmAnalysis {
    func overlayValue(kind: String, at time: Double) -> Double? {
        guard let ov = overlays.first(where: { $0.kind == kind && $0.channels == 1 }),
              !ov.samples.isEmpty else { return nil }
        if let sample = ov.samples.last(where: { $0.time <= time }) { return sample.value }
        return ov.samples.first?.value
    }

    func emotionProgress(at time: Double) -> Double? {
        overlayValue(kind: "emotion_progress", at: time)
    }
}
