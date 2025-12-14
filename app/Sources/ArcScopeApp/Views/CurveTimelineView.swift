import Charts
import Foundation
import SwiftUI

struct CurveTimelineView: View {
    let analysis: FilmAnalysis
    let referenceAnalysis: FilmAnalysis?
    @Binding var hoveredSecond: Double?
    @Binding var selectedRange: TimeRangeSelection?
    @State private var visibleKinds: Set<CurveKind> = [.pace, .sound, .color, .info, .arousal, .faceAffect]

    var body: some View {
        let visibleCurves = analysis.curves.filter { visibleKinds.contains($0.kind) }
        let ghost = ghostReference(for: analysis)
        let ghostCurves = ghost?.curves.filter { visibleKinds.contains($0.kind) }
        return Chart {
            ForEach(visibleCurves) { curve in
                curveMarks(curve)
            }
            if let ghostCurves {
                ForEach(ghostCurves) { curve in
                    ghostCurveMarks(curve)
                }
            }
        }
        .chartXAxis {
            AxisMarks(position: .bottom) { value in
                AxisGridLine().foregroundStyle(Color.white.opacity(0.1))
                AxisTick().foregroundStyle(.white.opacity(0.6))
                AxisValueLabel(format: FloatingPointFormatStyle<Double>.number.precision(.fractionLength(0)))
                    .foregroundStyle(.white.opacity(0.8))
            }
        }
        .chartYAxis(.hidden)
        .chartXScale(domain: 0...analysis.duration)
        .padding(.horizontal, -8)
        .background(SwiftUI.Color(.sRGB, red: 0, green: 0, blue: 0, opacity: 0.2))
        .clipShape(RoundedRectangle(cornerRadius: 16, style: .continuous))
        .overlay(alignment: .topLeading) {
            legend
                .padding(12)
        }
        .chartOverlay { proxy in
            GeometryReader { geo in
                ZStack(alignment: .topLeading) {
                    if let selectedRange,
                       let x0 = proxy.position(forX: selectedRange.start),
                       let x1 = proxy.position(forX: selectedRange.end) {
                        let left = min(x0, x1)
                        let width = max(2, abs(x1 - x0))
                        Rectangle()
                            .fill(Color.white.opacity(0.08))
                            .frame(width: width)
                            .offset(x: left)
                            .overlay(
                                Rectangle().stroke(Color.white.opacity(0.18), lineWidth: 1)
                                    .frame(width: width)
                                    .offset(x: left)
                            )
                    }

                    Rectangle()
                        .fill(Color.clear)
                        .contentShape(Rectangle())
                        .gesture(
                            DragGesture(minimumDistance: 0)
                                .onChanged { value in
                                    let ratio = max(0, min(1, value.location.x / geo.size.width))
                                    hoveredSecond = ratio * analysis.duration
                                }
                                .onEnded { _ in hoveredSecond = nil }
                        )
                        .simultaneousGesture(
                            TapGesture(count: 2).onEnded {
                                selectedRange = nil
                            }
                        )
                }
            }
        }
        .overlay(alignment: .bottomLeading) {
            if let hoveredSecond {
                HUDView(analysis: analysis,
                        referenceTitle: ghost?.referenceTitle,
                        referenceTimeSeconds: ghost?.referenceTimeSeconds(hoveredSecond),
                        ghostCurves: ghostCurves,
                        hoveredSecond: hoveredSecond,
                        visibleKinds: visibleKinds)
                    .padding(12)
            }
        }
        .overlay(alignment: .topTrailing) {
            if let selectedRange {
                VStack(alignment: .trailing, spacing: 4) {
                    Text("\(selectedRange.label)")
                        .font(.caption.weight(.semibold))
                    Text(verbatim: "\(TimelineFormatting.timePositional(selectedRange.start))–\(TimelineFormatting.timePositional(selectedRange.end))  (\(TimelineFormatting.timePositional(selectedRange.duration)))")
                        .font(.caption2.monospacedDigit())
                        .foregroundStyle(.secondary)
                    if let u0 = analysis.emotionProgress(at: selectedRange.start),
                       let u1 = analysis.emotionProgress(at: selectedRange.end) {
                        Text(verbatim: L10n.format("timeline.selected.u_range",
                                                   TimelineFormatting.percent01(u0),
                                                   TimelineFormatting.percent01(u1)))
                            .font(.caption2.monospacedDigit())
                            .foregroundStyle(.secondary)
                    }
                }
                .padding(10)
                .background(.ultraThinMaterial, in: RoundedRectangle(cornerRadius: 12, style: .continuous))
                .padding(12)
            }
        }
    }

    private var legend: some View {
        HStack(spacing: 12) {
            ForEach(CurveKind.allCases) { kind in
                Button {
                    if visibleKinds.contains(kind) {
                        visibleKinds.remove(kind)
                    } else {
                        visibleKinds.insert(kind)
                    }
                } label: {
                    HStack(spacing: 4) {
                        Circle()
                            .fill(kind.color.opacity(visibleKinds.contains(kind) ? 1.0 : 0.25))
                            .frame(width: 8, height: 8)
                        Text(kind.label)
                        .font(.caption2)
                        .textCase(.uppercase)
                        .foregroundStyle(.secondary)
                        .opacity(visibleKinds.contains(kind) ? 1.0 : 0.55)
                    }
                }
                .buttonStyle(.plain)
            }
        }
    }

    @ChartContentBuilder
    private func curveMarks(_ curve: FilmAnalysis.Curve) -> some ChartContent {
        ForEach(curve.samples, id: \.time) { sample in
            LineMark(
                x: .value("Time", sample.time),
                y: .value(curve.kind.label, sample.value)
            )
            .interpolationMethod(.catmullRom)
            .foregroundStyle(curve.color.gradient)
            .lineStyle(StrokeStyle(lineWidth: curve.kind == .arousal ? 3 : 1.8))
            .opacity(curve.kind == .arousal ? 0.95 : 0.75)
        }
    }

    @ChartContentBuilder
    private func ghostCurveMarks(_ curve: FilmAnalysis.Curve) -> some ChartContent {
        ForEach(curve.samples, id: \.time) { sample in
            LineMark(
                x: .value("Time", sample.time),
                y: .value(curve.kind.label, sample.value)
            )
            .interpolationMethod(.catmullRom)
            .foregroundStyle(curve.kind.color.opacity(0.35))
            .lineStyle(StrokeStyle(lineWidth: 1.2, dash: [5, 5]))
        }
    }

    private struct GhostReference {
        let referenceTitle: String
        let curves: [FilmAnalysis.Curve]
        let referenceTimeSeconds: (Double) -> Double?
    }

    private func ghostReference(for analysis: FilmAnalysis) -> GhostReference? {
        guard let referenceAnalysis else { return nil }

        // Doc 10.6: align by emotional progress u(t) (doc 4.6.5)
        if let uCur = emotionProgress(from: analysis),
           let uRef = emotionProgress(from: referenceAnalysis),
           let mapping = buildUInverseMapping(uCur: uCur, uRef: uRef) {
            let curves = referenceAnalysis.curves.map { refCurve in
                let ghostSamples = resampleCurve(refCurve.samples, mapping: mapping)
                return FilmAnalysis.Curve(kind: refCurve.kind, samples: ghostSamples)
            }
            let refTimeAt: (Double) -> Double? = { t in
                let k = Int(floor(t))
                if k < 0 || k >= mapping.count { return nil }
                return mapping[k].refTime
            }
            return GhostReference(referenceTitle: referenceAnalysis.title,
                                  curves: curves,
                                  referenceTimeSeconds: refTimeAt)
        }

        // Fallback: linear time scaling (only if u(t) is missing)
        guard referenceAnalysis.duration > 0 else { return nil }
        let scale = analysis.duration / referenceAnalysis.duration
        let curves = referenceAnalysis.curves.map { curve in
            let scaledSamples = curve.samples.map { sample in
                FilmAnalysis.Sample(time: sample.time * scale, value: sample.value)
            }
            return FilmAnalysis.Curve(kind: curve.kind, samples: scaledSamples)
        }
        let refTimeAt: (Double) -> Double? = { t in
            if !(scale > 1e-9) { return nil }
            return t / scale
        }
        return GhostReference(referenceTitle: referenceAnalysis.title,
                              curves: curves,
                              referenceTimeSeconds: refTimeAt)
    }

    private func emotionProgress(from analysis: FilmAnalysis) -> [Double]? {
        guard let ov = analysis.overlays.first(where: { $0.kind == "emotion_progress" }),
              ov.channels == 1 else { return nil }
        let values = ov.samples.map(\.value)
        if values.isEmpty { return nil }
        return values
    }

    private struct UMapPoint {
        let j0: Int
        let j1: Int
        let w: Double
        let time: Double
        let refTime: Double
    }

    private func buildUInverseMapping(uCur: [Double], uRef: [Double]) -> [UMapPoint]? {
        let N = uCur.count
        guard N > 0, uRef.count > 0 else { return nil }
        var map: [UMapPoint] = []
        map.reserveCapacity(N)

        var j = 0
        for k in 0..<N {
            let u = max(0.0, min(1.0, uCur[k]))
            while j + 1 < uRef.count && uRef[j + 1] < u { j += 1 }

            let j0 = j
            let j1 = min(j + 1, uRef.count - 1)
            let u0 = uRef[j0]
            let u1 = uRef[j1]
            let w = (j1 != j0 && u1 > u0) ? ((u - u0) / (u1 - u0)) : 0.0
            let ww = max(0.0, min(1.0, w))
            let refIndex = Double(j0) + ww * Double(j1 - j0)
            map.append(UMapPoint(j0: j0, j1: j1, w: ww, time: Double(k) + 0.5, refTime: refIndex + 0.5))
        }
        if let first = map.first, let last = map.last, !(last.time > first.time) {
            return nil
        }
        return map
    }

    private func resampleCurve(_ refSamples: [FilmAnalysis.Sample], mapping: [UMapPoint]) -> [FilmAnalysis.Sample] {
        guard !refSamples.isEmpty else { return mapping.map { FilmAnalysis.Sample(time: $0.time, value: 0) } }
        let refValues = refSamples.map(\.value)
        let n = refValues.count

        return mapping.map { m in
            let a = refValues[min(max(m.j0, 0), n - 1)]
            let b = refValues[min(max(m.j1, 0), n - 1)]
            let v = a + m.w * (b - a)
            return FilmAnalysis.Sample(time: m.time, value: v)
        }
    }
}

private struct HUDView: View {
    let analysis: FilmAnalysis
    let referenceTitle: String?
    let referenceTimeSeconds: Double?
    let ghostCurves: [FilmAnalysis.Curve]?
    let hoveredSecond: Double
    let visibleKinds: Set<CurveKind>

    var body: some View {
        let rows = analysis.curves.filter { visibleKinds.contains($0.kind) }
        return VStack(alignment: .leading, spacing: 6) {
            header
            ForEach(rows) { curve in
                curveRow(curve)
            }
        }
        .padding(12)
        .background(.ultraThinMaterial, in: RoundedRectangle(cornerRadius: 12, style: .continuous))
    }

    @ViewBuilder
    private var header: some View {
        VStack(alignment: .leading, spacing: 2) {
            let tLabel = TimelineFormatting.timePositional(hoveredSecond)
            if let u = analysis.emotionProgress(at: hoveredSecond) {
                Text(verbatim: L10n.format("timeline.hud.t_u",
                                           tLabel,
                                           TimelineFormatting.percent01(u)))
                    .font(.headline)
            } else {
                Text(verbatim: L10n.format("timeline.hud.t", tLabel))
                    .font(.headline)
            }
            HStack(spacing: 8) {
                Text(verbatim: TimelineFormatting.secondsShort(hoveredSecond))
                    .font(.caption2.monospacedDigit())
                    .foregroundStyle(.secondary)
                if let referenceTitle,
                   let referenceTimeSeconds {
                    Text(verbatim: L10n.format("timeline.hud.ref",
                                               referenceTitle,
                                               TimelineFormatting.timePositional(referenceTimeSeconds)))
                        .font(.caption2)
                        .foregroundStyle(.secondary)
                }
            }
        }
    }

    @ViewBuilder
    private func curveRow(_ curve: FilmAnalysis.Curve) -> some View {
        let curValue = value(at: hoveredSecond, in: curve)
        let ghost = ghostCurves?.first(where: { $0.kind == curve.kind })
        let refValue = ghost.map { value(at: hoveredSecond, in: $0) }

        HStack {
            Text(curve.kind.label)
                .foregroundStyle(curve.color)
            Spacer()
            if let refValue {
                Text(verbatim: "\(TimelineFormatting.decimal(curValue, maxFractionDigits: 2))  |  \(TimelineFormatting.decimal(refValue, maxFractionDigits: 2))")
                    .foregroundStyle(.secondary)
            } else {
                Text(verbatim: TimelineFormatting.decimal(curValue, maxFractionDigits: 2))
                    .foregroundStyle(.secondary)
            }
        }
        .font(.caption)
    }

    private func value(at second: Double, in curve: FilmAnalysis.Curve) -> Double {
        guard let sample = curve.samples.last(where: { $0.time <= second }) else {
            return curve.samples.first?.value ?? 0
        }
        return sample.value
    }
}
