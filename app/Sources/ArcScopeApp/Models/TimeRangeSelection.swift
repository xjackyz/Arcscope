import Foundation

struct TimeRangeSelection: Equatable {
    let kind: Kind
    let start: Double
    let end: Double
    let label: String

    enum Kind: String, Codable {
        case shot
        case scene
        case sequence
        case act
        case issue
    }

    var duration: Double { max(0, end - start) }
}

extension TimeRangeSelection.Kind {
    var uiIcon: String {
        switch self {
        case .shot: return "🎞"
        case .scene: return "🧩"
        case .sequence: return "🧱"
        case .act: return "Ⅰ"
        case .issue: return ""
        }
    }

    var uiTitle: String {
        switch self {
        case .shot: return "Shot"
        case .scene: return "Scene"
        case .sequence: return "Sequence"
        case .act: return "Act"
        case .issue: return "Issue"
        }
    }

    func uiLabel(detail: String? = nil) -> String {
        let d = (detail ?? "").trimmingCharacters(in: .whitespacesAndNewlines)
        let base = uiIcon.isEmpty ? uiTitle : "\(uiIcon) \(uiTitle)"
        if d.isEmpty { return base }
        return "\(base): \(d)"
    }
}
