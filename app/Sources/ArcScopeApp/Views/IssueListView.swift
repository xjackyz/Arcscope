import SwiftUI

struct IssueListView: View {
    let issues: [FilmAnalysis.Issue]

    var body: some View {
        VStack(alignment: .leading, spacing: 8) {
            Text("diagnostics.title")
                .font(.title2.bold())
            if issues.isEmpty {
                Text("diagnostics.none")
                    .foregroundStyle(.secondary)
            } else {
                ScrollView(.horizontal, showsIndicators: false) {
                    HStack(spacing: 12) {
                        ForEach(issues) { issue in
                            let presentation = IssuePresentation.forType(issue.type)
                            VStack(alignment: .leading, spacing: 4) {
                                Text(presentation.shortLabel)
                                    .font(.headline)
                                Text(issue.explanation)
                                    .font(.footnote)
                                    .foregroundStyle(.secondary)
                                Text(verbatim: "\(TimelineFormatting.secondsShort(issue.start)) – \(TimelineFormatting.secondsShort(issue.end))")
                                    .font(.caption)
                                    .foregroundStyle(.secondary)
                            }
                            .padding(12)
                            .frame(width: 220, alignment: .leading)
                            .background(presentation.color.opacity(0.15), in: RoundedRectangle(cornerRadius: 12, style: .continuous))
                            .overlay(
                                RoundedRectangle(cornerRadius: 12, style: .continuous)
                                    .stroke(presentation.color.opacity(0.6 + issue.severity * 0.4), lineWidth: 1)
                            )
                        }
                    }
                }
            }
        }
    }

    private func color(for type: String) -> Color {
        IssuePresentation.forType(type).color
    }
}
