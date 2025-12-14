import AVFoundation
import SwiftUI
#if canImport(UniformTypeIdentifiers)
import UniformTypeIdentifiers
#endif
#if canImport(AppKit)
import AppKit
#endif

struct PosterExportView: View {
    let analysis: FilmAnalysis

    @Environment(\.dismiss) private var dismiss
    @StateObject private var model = PosterExportViewModel()

    private let posterSize = CGSize(width: 980, height: 1380)

    var body: some View {
        VStack(spacing: 12) {
            header
            ScrollView {
                PosterDocumentView(analysis: analysis,
                                   posterSize: posterSize,
                                   thumbnails: model.thumbnailsBySecond)
                    .frame(width: posterSize.width, height: posterSize.height)
                    .background(Color.black.opacity(0.92))
                    .clipShape(RoundedRectangle(cornerRadius: 18, style: .continuous))
                    .overlay(
                        RoundedRectangle(cornerRadius: 18, style: .continuous)
                            .stroke(Color.white.opacity(0.12), lineWidth: 1)
                    )
            }
        }
        .padding(18)
        .frame(minWidth: 1040, minHeight: 860)
        .task {
            await model.prepare(for: analysis)
        }
    }

    private var header: some View {
        HStack(spacing: 10) {
            VStack(alignment: .leading, spacing: 3) {
                Text("export.poster.title")
                    .font(.title3.weight(.semibold))
                Text("export.poster.subtitle")
                    .font(.callout)
                    .foregroundStyle(.secondary)
            }
            Spacer()
            if model.isLoading {
                ProgressView()
                    .progressViewStyle(.circular)
            }
            Button("export.poster.export_png") {
                model.exportPNG(analysis: analysis, posterSize: posterSize)
            }
            .buttonStyle(.borderedProminent)
            Button("export.poster.export_pdf") {
                model.exportPDF(analysis: analysis, posterSize: posterSize)
            }
            .buttonStyle(.bordered)
            Button("export.close") { dismiss() }
                .buttonStyle(.bordered)
        }
    }
}

private struct PosterDocumentView: View {
    let analysis: FilmAnalysis
    let posterSize: CGSize
    let thumbnails: [Int: CGImage]

    var body: some View {
        VStack(alignment: .leading, spacing: 16) {
            header
            curveStack
            keyframes
            footer
        }
        .padding(28)
        .frame(width: posterSize.width, height: posterSize.height, alignment: .topLeading)
    }

    private var header: some View {
        VStack(alignment: .leading, spacing: 6) {
            Text(analysis.title)
                .font(.system(size: 34, weight: .semibold, design: .default))
            HStack(spacing: 10) {
                Text(verbatim: TimelineFormatting.filmMetaLine(year: analysis.year,
                                                              director: analysis.director,
                                                              durationSeconds: analysis.duration))
            }
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

    private var curveStack: some View {
        VStack(alignment: .leading, spacing: 12) {
            Text("export.poster.core_curves")
                .font(.headline)
            VStack(spacing: 10) {
                curveRow(kind: .pace)
                curveRow(kind: .sound)
                curveRow(kind: .color)
                curveRow(kind: .info)
                curveRow(kind: .arousal)
                if analysis.curves.contains(where: { $0.kind == .faceAffect }) {
                    curveRow(kind: .faceAffect)
                }
            }
        }
    }

    private func curveRow(kind: CurveKind) -> some View {
        let curve = analysis.curves.first(where: { $0.kind == kind })
        return HStack(spacing: 10) {
            Text(kind.label)
                .font(.caption.weight(.semibold))
                .foregroundStyle(.secondary)
                .frame(width: 88, alignment: .leading)
            CurveMiniChart(samples: curve?.samples ?? [], color: kind.color)
                .frame(height: 40)
                .background(Color.white.opacity(0.05), in: RoundedRectangle(cornerRadius: 10, style: .continuous))
        }
    }

    private var keyframes: some View {
        VStack(alignment: .leading, spacing: 10) {
            Text("export.poster.keyframes")
                .font(.headline)
            let secs = keyframeSeconds()
            if secs.isEmpty {
                Text("export.poster.keyframes.unavailable")
                    .font(.callout)
                    .foregroundStyle(.secondary)
            } else {
                LazyVGrid(columns: Array(repeating: GridItem(.flexible(), spacing: 10), count: 3), spacing: 10) {
                    ForEach(secs, id: \.self) { sec in
                        let img = thumbnails[sec]
                        ZStack {
                            if let img {
                                Image(decorative: img, scale: 1.0, orientation: .up)
                                    .resizable()
                                    .scaledToFill()
                                    .frame(height: 160)
                                    .clipped()
                            } else {
                                Rectangle()
                                    .fill(Color.white.opacity(0.06))
                                    .frame(height: 160)
                            }
                            VStack {
                                Spacer()
                                HStack {
                                    Text(verbatim: TimelineFormatting.timePositional(Double(sec)))
                                        .font(.caption2.monospacedDigit().weight(.semibold))
                                        .padding(.horizontal, 8)
                                        .padding(.vertical, 4)
                                        .background(.ultraThinMaterial, in: Capsule(style: .continuous))
                                    Spacer()
                                }
                                .padding(8)
                            }
                        }
                        .clipShape(RoundedRectangle(cornerRadius: 12, style: .continuous))
                        .overlay(
                            RoundedRectangle(cornerRadius: 12, style: .continuous)
                                .stroke(Color.white.opacity(0.12), lineWidth: 1)
                        )
                    }
                }
            }
        }
    }

    private var footer: some View {
        VStack(alignment: .leading, spacing: 6) {
            Divider().overlay(Color.white.opacity(0.12))
            Text("export.poster.footer")
                .font(.caption.weight(.semibold))
                .foregroundStyle(.secondary)
        }
        .padding(.top, 4)
    }

    private func keyframeSeconds() -> [Int] {
        guard let path = analysis.filePath, !path.isEmpty else { return [] }
        let pace = analysis.curves.first(where: { $0.kind == .pace })?.samples ?? []
        let arousal = analysis.curves.first(where: { $0.kind == .arousal })?.samples ?? []
        let peaks = peaksUnionSeconds(pace, arousal)
        return peaks
    }

    private func peaksUnionSeconds(_ a: [FilmAnalysis.Sample], _ b: [FilmAnalysis.Sample]) -> [Int] {
        var out = Set<Int>()
        for sec in peakSeconds(from: a) { out.insert(sec) }
        for sec in peakSeconds(from: b) { out.insert(sec) }
        return out.sorted()
    }

    private func peakSeconds(from samples: [FilmAnalysis.Sample]) -> [Int] {
        guard samples.count >= 3 else { return [] }
        let values = samples.map(\.value).sorted()
        let thr = quantile(values, 0.95)

        var peaks: [Int] = []
        var i = 1
        while i + 1 < samples.count {
            let v0 = samples[i - 1].value
            let v1 = samples[i].value
            let v2 = samples[i + 1].value
            if v1 >= thr, v1 >= v0, v1 >= v2 {
                peaks.append(Int(floor(samples[i].time)))
            }
            i += 1
        }
        // Dedup near neighbors (no fixed seconds; collapse contiguous seconds).
        peaks.sort()
        var collapsed: [Int] = []
        for p in peaks {
            if let last = collapsed.last, p <= last + 1 { continue }
            collapsed.append(p)
        }
        return collapsed
    }

    private func quantile(_ sorted: [Double], _ q: Double) -> Double {
        guard !sorted.isEmpty else { return 0.0 }
        let qq = min(max(q, 0.0), 1.0)
        let pos = qq * Double(sorted.count - 1)
        let k = Int(floor(pos))
        let k2 = min(k + 1, sorted.count - 1)
        let frac = pos - Double(k)
        return sorted[k] + frac * (sorted[k2] - sorted[k])
    }

    // Time formatting is handled by TimelineFormatting.
}

private struct CurveMiniChart: View {
    let samples: [FilmAnalysis.Sample]
    let color: Color

    var body: some View {
        GeometryReader { geo in
            Canvas { ctx, size in
                guard samples.count >= 2 else { return }
                let maxT = samples.last?.time ?? 1.0
                let denT = max(maxT, 1e-6)
                var path = Path()
                for (idx, s) in samples.enumerated() {
                    let x = CGFloat(s.time / denT) * size.width
                    let y = (1.0 - CGFloat(min(max(s.value, 0.0), 1.0))) * size.height
                    if idx == 0 { path.move(to: CGPoint(x: x, y: y)) }
                    else { path.addLine(to: CGPoint(x: x, y: y)) }
                }
                ctx.stroke(path, with: .color(color.opacity(0.85)), lineWidth: 2.0)
            }
        }
        .padding(.horizontal, 10)
        .padding(.vertical, 8)
    }
}

@MainActor
private final class PosterExportViewModel: ObservableObject {
    @Published var thumbnailsBySecond: [Int: CGImage] = [:]
    @Published var isLoading: Bool = false
    private let bookmarkStore = SecurityScopedBookmarkStore.shared

    func prepare(for analysis: FilmAnalysis) async {
        thumbnailsBySecond = [:]
        guard let path = analysis.filePath, !path.isEmpty else { return }
        isLoading = true
        defer { isLoading = false }

        let seconds = keyframeSeconds(for: analysis)
        guard !seconds.isEmpty else { return }
        let url = bookmarkStore.restoreVideoFile(forPath: path) ?? URL(fileURLWithPath: path)
        thumbnailsBySecond = await SecurityScopedAccess.withAccessAsync(url) {
            await generateThumbnails(url: url, seconds: seconds)
        }
    }

    func exportPNG(analysis: FilmAnalysis, posterSize: CGSize) {
        export(kind: .png, analysis: analysis, posterSize: posterSize)
    }

    func exportPDF(analysis: FilmAnalysis, posterSize: CGSize) {
        export(kind: .pdf, analysis: analysis, posterSize: posterSize)
    }

    private enum ExportKind { case png, pdf }

    private func export(kind: ExportKind, analysis: FilmAnalysis, posterSize: CGSize) {
        #if canImport(AppKit)
        let panel = NSSavePanel()
        panel.canCreateDirectories = true
        panel.nameFieldStringValue = "\(analysis.title.replacingOccurrences(of: "/", with: "_"))-poster.\(kind == .png ? "png" : "pdf")"
        panel.allowedContentTypes = kind == .png ? [.png] : [.pdf]
        panel.begin { response in
            guard response == .OK, let url = panel.url else { return }
            Task { @MainActor in
                let poster = PosterDocumentView(analysis: analysis, posterSize: posterSize, thumbnails: self.thumbnailsBySecond)
                    .frame(width: posterSize.width, height: posterSize.height)
                    .background(Color.black.opacity(0.92))
                switch kind {
                case .png:
                    let renderer = ImageRenderer(content: poster)
                    renderer.scale = 2
                    guard let nsImage = renderer.nsImage,
                          let tiff = nsImage.tiffRepresentation,
                          let rep = NSBitmapImageRep(data: tiff),
                          let png = rep.representation(using: .png, properties: [:]) else { return }
                    try? png.write(to: url)
                case .pdf:
                    let host = NSHostingView(rootView: poster)
                    host.frame = CGRect(origin: .zero, size: posterSize)
                    let data = host.dataWithPDF(inside: host.bounds)
                    try? data.write(to: url)
                }
            }
        }
        #endif
    }

    private func keyframeSeconds(for analysis: FilmAnalysis) -> [Int] {
        let pace = analysis.curves.first(where: { $0.kind == .pace })?.samples ?? []
        let arousal = analysis.curves.first(where: { $0.kind == .arousal })?.samples ?? []
        var set = Set<Int>()
        for sec in peakSeconds(from: pace) { set.insert(sec) }
        for sec in peakSeconds(from: arousal) { set.insert(sec) }
        return set.sorted()
    }

    private func peakSeconds(from samples: [FilmAnalysis.Sample]) -> [Int] {
        guard samples.count >= 3 else { return [] }
        let sorted = samples.map(\.value).sorted()
        let thr = quantile(sorted, 0.95)
        var peaks: [Int] = []
        for i in 1..<(samples.count - 1) {
            let v0 = samples[i - 1].value
            let v1 = samples[i].value
            let v2 = samples[i + 1].value
            if v1 >= thr, v1 >= v0, v1 >= v2 {
                peaks.append(Int(floor(samples[i].time)))
            }
        }
        peaks.sort()
        var out: [Int] = []
        for p in peaks {
            if let last = out.last, p <= last + 1 { continue }
            out.append(p)
        }
        return out
    }

    private func quantile(_ sorted: [Double], _ q: Double) -> Double {
        guard !sorted.isEmpty else { return 0.0 }
        let qq = min(max(q, 0.0), 1.0)
        let pos = qq * Double(sorted.count - 1)
        let k = Int(floor(pos))
        let k2 = min(k + 1, sorted.count - 1)
        let frac = pos - Double(k)
        return sorted[k] + frac * (sorted[k2] - sorted[k])
    }

    private func generateThumbnails(url: URL, seconds: [Int]) async -> [Int: CGImage] {
        await Task.detached(priority: .utility) { () -> [Int: CGImage] in
            let asset = AVURLAsset(url: url)
            let gen = AVAssetImageGenerator(asset: asset)
            gen.appliesPreferredTrackTransform = true
            gen.requestedTimeToleranceAfter = .zero
            gen.requestedTimeToleranceBefore = .zero

            var out: [Int: CGImage] = [:]
            for sec in seconds {
                let t = CMTime(seconds: Double(sec), preferredTimescale: 600)
                if let cg = try? gen.copyCGImage(at: t, actualTime: nil) {
                    out[sec] = cg
                }
            }
            return out
        }.value
    }
}
