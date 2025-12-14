import Foundation

enum TimelineFormatting {
    static func timePositional(_ seconds: Double) -> String {
        let formatter = DateComponentsFormatter()
        formatter.unitsStyle = .positional
        formatter.zeroFormattingBehavior = [.pad]
        let s = max(0, Int(seconds.rounded()))
        formatter.allowedUnits = (s >= 3600) ? [.hour, .minute, .second] : [.minute, .second]
        return formatter.string(from: TimeInterval(s)) ?? "0:00"
    }

    static func percent01(_ u: Double) -> String {
        if !u.isFinite { return "—" }
        let number = NSNumber(value: max(0.0, min(1.0, u)))
        let formatter = NumberFormatter()
        formatter.numberStyle = .percent
        formatter.maximumFractionDigits = 0
        formatter.minimumFractionDigits = 0
        return formatter.string(from: number) ?? "0%"
    }

    static func decimal(_ value: Double, maxFractionDigits: Int) -> String {
        let formatter = NumberFormatter()
        formatter.numberStyle = .decimal
        formatter.maximumFractionDigits = max(0, maxFractionDigits)
        formatter.minimumFractionDigits = 0
        return formatter.string(from: NSNumber(value: value)) ?? "\(value)"
    }

    static func significant(_ value: Double, maxSignificantDigits: Int) -> String {
        let formatter = NumberFormatter()
        formatter.numberStyle = .decimal
        formatter.usesSignificantDigits = true
        formatter.maximumSignificantDigits = max(1, maxSignificantDigits)
        formatter.minimumSignificantDigits = 1
        return formatter.string(from: NSNumber(value: value)) ?? "\(value)"
    }

    static func durationAbbreviated(_ seconds: Double) -> String {
        let formatter = DateComponentsFormatter()
        formatter.unitsStyle = .abbreviated
        formatter.allowedUnits = [.hour, .minute]
        formatter.maximumUnitCount = 2
        let s = max(0, Int(seconds.rounded()))
        return formatter.string(from: TimeInterval(s)) ?? "0m"
    }

    static func secondsShort(_ seconds: Double) -> String {
        let s = max(0, Int(seconds.rounded()))
        return L10n.format("term.seconds_short", s)
    }

    static func secondsDecimalShort(_ seconds: Double, maxFractionDigits: Int) -> String {
        let number = decimal(seconds, maxFractionDigits: maxFractionDigits)
        return L10n.format("term.seconds_value", number)
    }

    static func filmMetaLine(year: Int, director: String, durationSeconds: Double) -> String {
        let yearPart = (year > 0) ? String(year) : "—"
        let directorPart = director.trimmingCharacters(in: .whitespacesAndNewlines).isEmpty ? "—" : director
        let durationPart = durationAbbreviated(durationSeconds)
        return "\(yearPart) · \(directorPart) · \(durationPart)"
    }

    static func structureCountsLine(shots: Int, scenes: Int, sequences: Int, acts: Int) -> String {
        L10n.format("export.structure.counts", shots, scenes, sequences, acts)
    }
}
