import SwiftUI

// arcscope.md §6.4: Diagnostics heatmap hover/click should use a consistent
// icon + text + color vocabulary across the UI (no drift between views).
struct IssuePresentation: Sendable {
    let type: String
    let icon: String
    let color: Color

    var displayName: String {
        let key = "issue.type.\(type)"
        let localized = NSLocalizedString(key, bundle: .module, comment: "")
        return (localized == key) ? type : localized
    }

    var shortLabel: String {
        if icon.isEmpty { return displayName }
        return "\(icon) \(displayName)"
    }

    static func forType(_ type: String) -> IssuePresentation {
        switch type {
        // LowActivity: blue family
        case "LowActivity":
            return IssuePresentation(type: type, icon: "⏳", color: .blue)

        // Overcut / high tension-low pace: orange family
        case "OvercutFlat":
            return IssuePresentation(type: type, icon: "✂︎", color: .orange)
        case "HighTensionLowPace":
            return IssuePresentation(type: type, icon: "⚡︎", color: .orange)

        // AV / sound dynamics: purple family
        case "AudioVisualMisalign":
            return IssuePresentation(type: type, icon: "🔀", color: .purple)
        case "SoundOverdriveFlatEmotion":
            return IssuePresentation(type: type, icon: "🔊", color: .purple)
        case "SilentHighTension":
            return IssuePresentation(type: type, icon: "🔇", color: .purple)

        // Info issues: red family; use the doc icon "…" for info.
        case "InfoKillsEmotion":
            return IssuePresentation(type: type, icon: "…", color: .red)
        case "CognitiveOverload":
            return IssuePresentation(type: type, icon: "…", color: .red)

        // Color / multimodal counterpoint: cyan-teal-indigo family
        case "SoundColorConflict":
            return IssuePresentation(type: type, icon: "🔊⇄🎨", color: .cyan)
        case "ColorEmotionMismatchWarm", "ColorEmotionMismatchCold":
            return IssuePresentation(type: type, icon: "🎨", color: .teal)
        case "ModeCounterpoint":
            return IssuePresentation(type: type, icon: "↯", color: .indigo)

        // Face-related: green family
        case "FacialAffectSuppressed", "NonCharacterDrivenPeak":
            return IssuePresentation(type: type, icon: "🙂", color: .green)

        default:
            return IssuePresentation(type: type, icon: "", color: .gray)
        }
    }
}
