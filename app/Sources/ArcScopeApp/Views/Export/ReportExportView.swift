import SwiftUI
#if canImport(UniformTypeIdentifiers)
import UniformTypeIdentifiers
#endif
#if canImport(AppKit)
import AppKit
#endif

struct ReportExportView: View {
    let analysis: FilmAnalysis

    @Environment(\.dismiss) private var dismiss
    @State private var markdown: String = ""

    var body: some View {
        VStack(spacing: 12) {
            header
            HStack(spacing: 12) {
                reportPreview
                markdownPreview
            }
        }
        .padding(18)
        .frame(minWidth: 1040, minHeight: 720)
        .task {
            markdown = ReportBuilder(analysis: analysis).markdown()
        }
    }

    private var header: some View {
        HStack(spacing: 10) {
            VStack(alignment: .leading, spacing: 3) {
                Text("export.report.title")
                    .font(.title3.weight(.semibold))
                Text("export.report.subtitle")
                    .font(.callout)
                    .foregroundStyle(.secondary)
            }
            Spacer()
            Button("export.report.export_markdown") { exportMarkdown() }
                .buttonStyle(.borderedProminent)
            Button("export.report.export_pdf") { exportPDF() }
                .buttonStyle(.bordered)
            Button("export.close") { dismiss() }
                .buttonStyle(.bordered)
        }
    }

    private var reportPreview: some View {
        ReportDocumentView(analysis: analysis)
            .background(Color.black.opacity(0.92))
            .clipShape(RoundedRectangle(cornerRadius: 16, style: .continuous))
            .overlay(
                RoundedRectangle(cornerRadius: 16, style: .continuous)
                    .stroke(Color.white.opacity(0.12), lineWidth: 1)
            )
    }

    private var markdownPreview: some View {
        VStack(alignment: .leading, spacing: 8) {
            Text("export.markdown_preview")
                .font(.headline)
                .foregroundStyle(.secondary)
            ScrollView {
                Text(markdown)
                    .font(.caption.monospaced())
                    .foregroundStyle(.white.opacity(0.85))
                    .frame(maxWidth: .infinity, alignment: .leading)
                    .padding(12)
            }
            .background(Color.white.opacity(0.06), in: RoundedRectangle(cornerRadius: 12, style: .continuous))
        }
        .frame(width: 420)
    }

    private func exportMarkdown() {
        #if canImport(AppKit)
        let panel = NSSavePanel()
        panel.canCreateDirectories = true
        panel.nameFieldStringValue = "\(safeFileName(analysis.title))-report.md"
        panel.allowedContentTypes = [UTType(filenameExtension: "md") ?? .plainText]
        panel.begin { response in
            guard response == .OK, let url = panel.url else { return }
            try? markdown.data(using: .utf8)?.write(to: url)
        }
        #endif
    }

    private func exportPDF() {
        #if canImport(AppKit)
        let panel = NSSavePanel()
        panel.canCreateDirectories = true
        panel.nameFieldStringValue = "\(safeFileName(analysis.title))-report.pdf"
        panel.allowedContentTypes = [.pdf]
        panel.begin { response in
            guard response == .OK, let url = panel.url else { return }
            Task { @MainActor in
                let view = ReportDocumentView(analysis: analysis)
                    .frame(width: 980, height: 1380)
                    .background(Color.black.opacity(0.92))
                let host = NSHostingView(rootView: view)
                host.frame = CGRect(x: 0, y: 0, width: 980, height: 1380)
                let data = host.dataWithPDF(inside: host.bounds)
                try? data.write(to: url)
            }
        }
        #endif
    }

    private func safeFileName(_ s: String) -> String {
        let cleaned = s.replacingOccurrences(of: "/", with: "_")
        return cleaned.isEmpty ? "arcscope" : cleaned
    }
}

private struct ReportDocumentView: View {
    let analysis: FilmAnalysis

    var body: some View {
        VStack(alignment: .leading, spacing: 14) {
            header
            curveStats
            structureStats
            issueStats
            Spacer(minLength: 8)
            footer
        }
        .padding(28)
        .frame(width: 980, height: 1380, alignment: .topLeading)
    }

    private var header: some View {
        VStack(alignment: .leading, spacing: 6) {
            Text(analysis.title)
                .font(.system(size: 34, weight: .semibold))
            Text(verbatim: TimelineFormatting.filmMetaLine(year: analysis.year,
                                                          director: analysis.director,
                                                          durationSeconds: analysis.duration))
                .font(.callout)
                .foregroundStyle(.secondary)
            if let path = analysis.filePath, !path.isEmpty {
                Text(path)
                    .font(.caption.monospaced())
                    .foregroundStyle(.secondary)
                    .lineLimit(1)
            }
        }
    }

    private var curveStats: some View {
        VStack(alignment: .leading, spacing: 8) {
            Text("export.section.curves")
                .font(.headline)
            ForEach(analysis.curves) { c in
                let stats = Stats(samples: c.samples)
                HStack {
                    Text(c.kind.label)
                        .font(.caption.weight(.semibold))
                        .frame(width: 110, alignment: .leading)
                    StatPill(title: "mean", value: stats.mean)
                    StatPill(title: "min", value: stats.min)
                    StatPill(title: "max", value: stats.max)
                    Spacer()
                }
            }
        }
        .padding(14)
        .background(Color.white.opacity(0.05), in: RoundedRectangle(cornerRadius: 14, style: .continuous))
    }

    private var structureStats: some View {
        VStack(alignment: .leading, spacing: 8) {
            Text("export.section.structure")
                .font(.headline)
            Text(verbatim: TimelineFormatting.structureCountsLine(shots: analysis.shots.count,
                                                                 scenes: analysis.scenes.count,
                                                                 sequences: analysis.sequences.count,
                                                                 acts: analysis.acts.count))
            .font(.callout)
            .foregroundStyle(.secondary)
        }
        .padding(14)
        .background(Color.white.opacity(0.05), in: RoundedRectangle(cornerRadius: 14, style: .continuous))
    }

    private var issueStats: some View {
        VStack(alignment: .leading, spacing: 10) {
            Text("export.section.diagnostics")
                .font(.headline)
            if analysis.issues.isEmpty {
                Text("export.no_issues")
                    .font(.callout)
                    .foregroundStyle(.secondary)
            } else {
                ForEach(analysis.issues.prefix(18)) { issue in
                    VStack(alignment: .leading, spacing: 2) {
                        HStack {
                            Text(issue.type)
                                .font(.caption.weight(.semibold))
                            Spacer()
                            Text(verbatim: "\(TimelineFormatting.timePositional(issue.start))–\(TimelineFormatting.timePositional(issue.end))  s=\(TimelineFormatting.decimal(issue.severity, maxFractionDigits: 2))")
                                .font(.caption.monospacedDigit())
                                .foregroundStyle(.secondary)
                        }
                        if !issue.explanation.isEmpty {
                            Text(issue.explanation)
                                .font(.caption)
                                .foregroundStyle(.secondary)
                                .lineLimit(2)
                        }
                    }
                    Divider().overlay(Color.white.opacity(0.10))
                }
                if analysis.issues.count > 18 {
                    Text(verbatim: L10n.format("export.more_issues", analysis.issues.count - 18))
                        .font(.caption)
                        .foregroundStyle(.secondary)
                }
            }
        }
        .padding(14)
        .background(Color.white.opacity(0.05), in: RoundedRectangle(cornerRadius: 14, style: .continuous))
    }

    private var footer: some View {
        VStack(alignment: .leading, spacing: 6) {
            Divider().overlay(Color.white.opacity(0.12))
            Text("export.report.footer")
                .font(.caption.weight(.semibold))
                .foregroundStyle(.secondary)
        }
    }

    private struct Stats {
        let mean: Double
        let min: Double
        let max: Double

        init(samples: [FilmAnalysis.Sample]) {
            let v = samples.map(\.value)
            if v.isEmpty {
                mean = 0
                min = 0
                max = 0
                return
            }
            let sum = v.reduce(0.0, +)
            mean = sum / Double(v.count)
            min = v.min() ?? 0
            max = v.max() ?? 0
        }
    }
}

private struct StatPill: View {
    let title: String
    let value: Double

    var body: some View {
        HStack(spacing: 6) {
            Text(title)
                .font(.caption2.weight(.semibold))
                .foregroundStyle(.secondary)
            Text(verbatim: TimelineFormatting.decimal(value, maxFractionDigits: 2))
                .font(.caption2.monospacedDigit().weight(.semibold))
        }
        .padding(.horizontal, 10)
        .padding(.vertical, 6)
        .background(.ultraThinMaterial, in: Capsule(style: .continuous))
    }
}

private struct ReportBuilder {
    let analysis: FilmAnalysis

    func markdown() -> String {
        var lines: [String] = []
        lines.append("# \(analysis.title)")
        lines.append("")
        lines.append("- \(L10n.string("report.md.year")): \(analysis.year > 0 ? String(analysis.year) : "—")")
        lines.append("- \(L10n.string("report.md.director")): \(analysis.director.isEmpty ? "—" : analysis.director)")
        lines.append("- \(L10n.string("report.md.duration")): \(TimelineFormatting.secondsShort(analysis.duration))")
        if let path = analysis.filePath, !path.isEmpty {
            lines.append("- \(L10n.string("report.md.file")): `\(path)`")
        }
        lines.append("")
        lines.append("## \(L10n.string("report.md.curves"))")
        for c in analysis.curves {
            let stats = stats(samples: c.samples)
            lines.append("- \(c.kind.rawValue): mean=\(fmt(stats.mean)) min=\(fmt(stats.min)) max=\(fmt(stats.max))")
        }
        lines.append("")
        lines.append("## \(L10n.string("report.md.structure"))")
        lines.append("- \(TimelineFormatting.structureCountsLine(shots: analysis.shots.count, scenes: analysis.scenes.count, sequences: analysis.sequences.count, acts: analysis.acts.count))")
        lines.append("")
        lines.append("## \(L10n.string("report.md.diagnostics"))")
        if analysis.issues.isEmpty {
            lines.append("- (\(L10n.string("report.md.none")))")
        } else {
            for issue in analysis.issues {
                lines.append("- \(issue.type) \(fmtTime(issue.start))–\(fmtTime(issue.end)) severity=\(fmt(issue.severity))")
                if !issue.explanation.isEmpty {
                    lines.append("  - \(issue.explanation)")
                }
            }
        }
        lines.append("")
        lines.append("_\(L10n.string("report.md.generated_by"))_")
        return lines.joined(separator: "\n")
    }

    private func stats(samples: [FilmAnalysis.Sample]) -> (mean: Double, min: Double, max: Double) {
        let v = samples.map(\.value)
        guard !v.isEmpty else { return (0, 0, 0) }
        let sum = v.reduce(0.0, +)
        return (sum / Double(v.count), v.min() ?? 0, v.max() ?? 0)
    }

    private func fmt(_ x: Double) -> String { TimelineFormatting.decimal(x, maxFractionDigits: 3) }

    private func fmtTime(_ sec: Double) -> String { TimelineFormatting.timePositional(sec) }
}
