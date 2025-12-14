import SwiftUI

struct StructureBandsView: View {
    let analysis: FilmAnalysis
    @Binding var hoveredSecond: Double?
    @Binding var selectedRange: TimeRangeSelection?

    var body: some View {
        VStack(spacing: 8) {
            if !analysis.shots.isEmpty {
                ShotGridBand(shots: analysis.shots,
                             duration: analysis.duration,
                             hoveredSecond: $hoveredSecond,
                             selectedRange: $selectedRange)
                    .frame(height: 14)
            }
            if !analysis.scenes.isEmpty {
                SceneBand(scenes: analysis.scenes,
                          duration: analysis.duration,
                          hoveredSecond: $hoveredSecond,
                          selectedRange: $selectedRange)
                    .frame(height: 26)
            }
            if !analysis.sequences.isEmpty {
                SpanBand(spans: analysis.sequences.enumerated().map { (idx, s) in (s.start, s.end, "\(idx + 1)") },
                         duration: analysis.duration,
                         fill: Color.white.opacity(0.10),
                         stroke: Color.white.opacity(0.18),
                         label: "Sequence",
                         hoveredSecond: $hoveredSecond,
                         selectedRange: $selectedRange,
                         selectionKind: .sequence)
                    .frame(height: 18)
            }
            if !analysis.acts.isEmpty {
                SpanBand(spans: analysis.acts.map { ($0.start, $0.end, "\($0.number)") },
                         duration: analysis.duration,
                         fill: Color.white.opacity(0.16),
                         stroke: Color.white.opacity(0.26),
                         label: "Act",
                         hoveredSecond: $hoveredSecond,
                         selectedRange: $selectedRange,
                         selectionKind: .act)
                    .frame(height: 22)
            }
        }
    }
}

private struct ShotGridBand: View {
    let shots: [FilmAnalysis.ShotSegment]
    let duration: Double
    @Binding var hoveredSecond: Double?
    @Binding var selectedRange: TimeRangeSelection?

    var body: some View {
        GeometryReader { geo in
            ZStack(alignment: .leading) {
                Color.black.opacity(0.18)
                Canvas { ctx, size in
                    guard duration > 0 else { return }

                    let durations = shots.map(\.duration).sorted()
                    let shortThr = durations.isEmpty ? 0.0 : durations[max(0, min(durations.count - 1, Int(round(0.10 * Double(durations.count - 1)))))]

                    // Draw cut boundaries as thin ticks at each shot start.
                    for (idx, shot) in shots.enumerated() {
                        guard shot.start > 0.0 else { continue }
                        let x = CGFloat(shot.start / duration) * size.width
                        var p = Path()
                        p.move(to: CGPoint(x: x, y: 0))
                        p.addLine(to: CGPoint(x: x, y: size.height))
                        let isShort = (shot.duration > 0.0 && shortThr > 0.0) ? (shot.duration <= shortThr) : false
                        let alpha = isShort ? 0.42 : 0.22
                        ctx.stroke(p, with: .color(Color.white.opacity(alpha)), lineWidth: isShort ? 1.4 : 1.0)

                        // If selected shot, draw a stronger marker at its start.
                        if let sel = selectedRange, sel.kind == .shot,
                           abs(sel.start - shot.start) < 1e-6, abs(sel.end - shot.end) < 1e-6 {
                            var p2 = Path()
                            p2.move(to: CGPoint(x: x, y: 0))
                            p2.addLine(to: CGPoint(x: x, y: size.height))
                            ctx.stroke(p2, with: .color(Color.white.opacity(0.65)), lineWidth: 2.0)
                        }

                        // Small implicit numbering cadence: every 20th shot gets a tiny tick cap.
                        if idx % 20 == 0 {
                            let cap = CGRect(x: x - 1, y: 0, width: 2, height: 3)
                            ctx.fill(Path(cap), with: .color(Color.white.opacity(0.28)))
                        }
                    }

                    // Subtle density cue: higher avgMotion -> brighter background stripe.
                    for shot in shots {
                        let x0 = CGFloat(shot.start / duration) * size.width
                        let x1 = CGFloat(shot.end / duration) * size.width
                        let w = max(1, x1 - x0)
                        let v = max(0.0, min(1.0, shot.avgMotion))
                        let r = CGRect(x: x0, y: 0, width: w, height: size.height)
                        ctx.fill(Path(r), with: .color(Color.white.opacity(0.03 + 0.10 * v)))
                    }
                }

                if let hoveredSecond {
                    let x = CGFloat(hoveredSecond / max(duration, 1)) * geo.size.width
                    Rectangle()
                        .fill(Color.white.opacity(0.5))
                        .frame(width: 1)
                        .offset(x: x)
                }
            }
            .clipShape(RoundedRectangle(cornerRadius: 8, style: .continuous))
            .contentShape(Rectangle())
            .gesture(
                DragGesture(minimumDistance: 0)
                    .onChanged { value in
                        let ratio = max(0, min(1, value.location.x / geo.size.width))
                        hoveredSecond = ratio * duration
                    }
                    .onEnded { _ in hoveredSecond = nil }
            )
            .onTapGesture {
                guard let t = hoveredSecond else { return }
                if let shotIndex = shots.firstIndex(where: { t >= $0.start && t < $0.end }) {
                    let shot = shots[shotIndex]
                    selectedRange = TimeRangeSelection(kind: .shot,
                                                       start: shot.start,
                                                       end: shot.end,
                                                       label: TimeRangeSelection.Kind.shot.uiLabel(detail: "\(shotIndex + 1)"))
                }
            }
            .overlay(alignment: .topTrailing) {
                if let t = hoveredSecond,
                   let shotIndex = shots.firstIndex(where: { t >= $0.start && t < $0.end }) {
                    let shot = shots[shotIndex]
                    Text(verbatim: L10n.format("structure.shot_hover",
                                               shotIndex + 1,
                                               TimelineFormatting.secondsDecimalShort(shot.duration, maxFractionDigits: 1)))
                        .font(.caption2.monospacedDigit())
                        .foregroundStyle(.secondary)
                        .padding(.horizontal, 8)
                        .padding(.vertical, 4)
                        .background(.ultraThinMaterial, in: Capsule(style: .continuous))
                        .padding(6)
                }
            }
        }
    }
}

private struct SceneBand: View {
    let scenes: [FilmAnalysis.SceneSegment]
    let duration: Double
    @Binding var hoveredSecond: Double?
    @Binding var selectedRange: TimeRangeSelection?

    var body: some View {
        GeometryReader { geo in
            ZStack(alignment: .leading) {
                ForEach(Array(scenes.enumerated()), id: \.offset) { index, scene in
                    let startX = CGFloat(scene.start / max(duration, 1)) * geo.size.width
                    let width = CGFloat(scene.duration / max(duration, 1)) * geo.size.width
                    RoundedRectangle(cornerRadius: 6, style: .continuous)
                        .fill(sceneColor(label: scene.label, index: index))
                        .frame(width: max(width, 2), height: geo.size.height)
                        .overlay(
                            RoundedRectangle(cornerRadius: 6, style: .continuous)
                                .stroke(Color.white.opacity(isSelected(scene) ? 0.55 : 0.0),
                                        lineWidth: isSelected(scene) ? 2 : 0)
                        )
                        .overlay(
                            Text(scene.label)
                                .font(.caption2)
                                .padding(4)
                                .foregroundStyle(.white.opacity(0.82))
                                .opacity(width > 120 ? 1 : 0),
                            alignment: .topLeading
                        )
                        .offset(x: startX)
                        .contentShape(Rectangle())
                        .onTapGesture {
                            selectedRange = TimeRangeSelection(kind: .scene,
                                                               start: scene.start,
                                                               end: scene.end,
                                                               label: TimeRangeSelection.Kind.scene.uiLabel(detail: scene.label))
                        }
                }
                if let hoveredSecond {
                    let x = CGFloat(hoveredSecond / max(duration, 1)) * geo.size.width
                    Rectangle()
                        .fill(Color.white.opacity(0.35))
                        .frame(width: 1)
                        .offset(x: x)
                }
            }
            .clipShape(RoundedRectangle(cornerRadius: 10, style: .continuous))
            .contentShape(Rectangle())
            .gesture(
                DragGesture(minimumDistance: 0)
                    .onChanged { value in
                        let ratio = max(0, min(1, value.location.x / geo.size.width))
                        hoveredSecond = ratio * duration
                    }
                    .onEnded { _ in hoveredSecond = nil }
            )
        }
    }

    private func isSelected(_ scene: FilmAnalysis.SceneSegment) -> Bool {
        guard let sel = selectedRange, sel.kind == .scene else { return false }
        return abs(sel.start - scene.start) < 1e-6 && abs(sel.end - scene.end) < 1e-6
    }

    private func sceneColor(label: String, index: Int) -> Color {
        switch label {
        case "LowEnergyBasin": return .blue.opacity(0.22)
        case "BuildUp": return .orange.opacity(0.20)
        case "HighEnergyPeak": return .red.opacity(0.22)
        case "EmotionalReversal": return .purple.opacity(0.22)
        case "MisalignedTension": return .mint.opacity(0.20)
        case "Release": return .teal.opacity(0.20)
        default:
            let base: [Color] = [.white.opacity(0.10), .white.opacity(0.14), .white.opacity(0.12)]
            return base[index % base.count]
        }
    }
}

private struct SpanBand: View {
    let spans: [(start: Double, end: Double, label: String)]
    let duration: Double
    let fill: Color
    let stroke: Color
    let label: String
    @Binding var hoveredSecond: Double?
    @Binding var selectedRange: TimeRangeSelection?
    let selectionKind: TimeRangeSelection.Kind

    var body: some View {
        GeometryReader { geo in
            ZStack(alignment: .leading) {
                ForEach(Array(spans.enumerated()), id: \.offset) { _, span in
                    let startX = CGFloat(span.start / max(duration, 1)) * geo.size.width
                    let width = CGFloat(max(0, span.end - span.start) / max(duration, 1)) * geo.size.width
                    RoundedRectangle(cornerRadius: 6, style: .continuous)
                        .fill(fill)
                        .overlay(
                            RoundedRectangle(cornerRadius: 6, style: .continuous)
                                .stroke(stroke, lineWidth: 1)
                        )
                        .overlay(
                            RoundedRectangle(cornerRadius: 6, style: .continuous)
                                .stroke(Color.white.opacity(isSelected(span) ? 0.55 : 0.0),
                                        lineWidth: isSelected(span) ? 2 : 0)
                        )
                        .frame(width: max(width, 2), height: geo.size.height)
                        .overlay(
                            Text(span.label.isEmpty ? label : span.label)
                                .font(.caption2)
                                .padding(4)
                                .foregroundStyle(.white.opacity(0.78))
                                .opacity(width > 120 ? 1 : 0),
                            alignment: .topLeading
                        )
                        .offset(x: startX)
                        .contentShape(Rectangle())
                        .onTapGesture {
                            let detail = span.label.isEmpty ? "" : span.label
                            selectedRange = TimeRangeSelection(kind: selectionKind,
                                                               start: span.start,
                                                               end: span.end,
                                                               label: selectionKind.uiLabel(detail: detail))
                        }
                }
                if let hoveredSecond {
                    let x = CGFloat(hoveredSecond / max(duration, 1)) * geo.size.width
                    Rectangle()
                        .fill(Color.white.opacity(0.28))
                        .frame(width: 1)
                        .offset(x: x)
                }
            }
            .clipShape(RoundedRectangle(cornerRadius: 10, style: .continuous))
            .contentShape(Rectangle())
            .gesture(
                DragGesture(minimumDistance: 0)
                    .onChanged { value in
                        let ratio = max(0, min(1, value.location.x / geo.size.width))
                        hoveredSecond = ratio * duration
                    }
                    .onEnded { _ in hoveredSecond = nil }
            )
        }
    }

    private func isSelected(_ span: (start: Double, end: Double, label: String)) -> Bool {
        guard let sel = selectedRange, sel.kind == selectionKind else { return false }
        return abs(sel.start - span.start) < 1e-6 && abs(sel.end - span.end) < 1e-6
    }
}
