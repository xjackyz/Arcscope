import Foundation
import SwiftUI

struct DiagnosticsHeatmap: View {
    let analysis: FilmAnalysis
    @Binding var hoveredSecond: Double?
    @Binding var selectedRange: TimeRangeSelection?

    var body: some View {
        GeometryReader { geo in
            ZStack(alignment: .leading) {
                Color.black.opacity(0.15)
                if let selectedRange {
                    let startX = CGFloat(selectedRange.start / analysis.duration) * geo.size.width
                    let width = CGFloat((selectedRange.end - selectedRange.start) / analysis.duration) * geo.size.width
                    RoundedRectangle(cornerRadius: 4, style: .continuous)
                        .fill(Color.white.opacity(0.10))
                        .frame(width: max(width, 2), height: geo.size.height)
                        .offset(x: startX)
                        .overlay(
                            RoundedRectangle(cornerRadius: 4, style: .continuous)
                                .stroke(Color.white.opacity(0.18), lineWidth: 1)
                                .frame(width: max(width, 2), height: geo.size.height)
                                .offset(x: startX)
                        )
                }
                ForEach(analysis.issues) { issue in
                    let presentation = IssuePresentation.forType(issue.type)
                    let startX = CGFloat(issue.start / analysis.duration) * geo.size.width
                    let width = CGFloat((issue.end - issue.start) / analysis.duration) * geo.size.width
                    RoundedRectangle(cornerRadius: 4, style: .continuous)
                        .fill(presentation.color.opacity(issue.severity))
                        .frame(width: max(width, 4), height: geo.size.height * 0.6)
                        .offset(x: startX, y: geo.size.height * 0.2)
                        .overlay(
                            Text(presentation.shortLabel)
                                .font(.caption2)
                                .padding(4)
                                .foregroundStyle(.white)
                                .background(Color.black.opacity(0.3), in: RoundedRectangle(cornerRadius: 6, style: .continuous))
                                .opacity(width > 120 ? 1 : 0)
                                .offset(y: -geo.size.height * 0.4)
                            , alignment: .topLeading
                        )
                        .contentShape(Rectangle())
                        .help(tooltipText(for: issue))
                        .onTapGesture {
                            selectedRange = TimeRangeSelection(kind: .issue,
                                                               start: issue.start,
                                                               end: issue.end,
                                                               label: presentation.shortLabel)
                        }
                }
                if let hoveredSecond {
                    let x = CGFloat(hoveredSecond / analysis.duration) * geo.size.width
                    Rectangle()
                        .fill(Color.white.opacity(0.5))
                        .frame(width: 1)
                        .offset(x: x)
                }
            }
            .clipShape(RoundedRectangle(cornerRadius: 10, style: .continuous))
        }
    }

    private func color(for type: String) -> Color {
        IssuePresentation.forType(type).color
    }

    private func shortLabel(for type: String) -> String {
        IssuePresentation.forType(type).shortLabel
    }

    private func icon(for type: String) -> String {
        IssuePresentation.forType(type).icon
    }

    private func tooltipText(for issue: FilmAnalysis.Issue) -> String {
        let dur = max(0, issue.end - issue.start)
        let p = IssuePresentation.forType(issue.type)
        let timeRange = "\(TimelineFormatting.timePositional(issue.start))–\(TimelineFormatting.timePositional(issue.end))"
        let duration = TimelineFormatting.secondsShort(dur)
        let severity = TimelineFormatting.decimal(issue.severity, maxFractionDigits: 2)
        let explanation = issue.explanation.trimmingCharacters(in: .whitespacesAndNewlines)
        if explanation.isEmpty {
            return "\(p.shortLabel)\n\(timeRange) (\(duration))\n\(L10n.string("diagnostics.severity")) \(severity)"
        }
        return "\(p.shortLabel)\n\(timeRange) (\(duration))\n\(L10n.string("diagnostics.severity")) \(severity)\n\(explanation)"
    }
}
