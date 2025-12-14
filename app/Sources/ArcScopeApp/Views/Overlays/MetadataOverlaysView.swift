import Charts
import SwiftUI

struct MetadataOverlaysView: View {
    let overlays: [FilmAnalysis.Overlay]
    let duration: Double
    @Binding var hoveredSecond: Double?
    @Binding var selectedRange: TimeRangeSelection?

    @State private var visible: Set<String> = []

    var body: some View {
        if overlays.isEmpty {
            EmptyView()
        } else {
            VStack(alignment: .leading, spacing: 10) {
                header
                ForEach(visibleOverlays, id: \.kind) { overlay in
                    overlayTrack(overlay)
                }
            }
        }
    }

    private var visibleOverlays: [FilmAnalysis.Overlay] {
        overlays.filter { visible.contains($0.kind) }
    }

    private var header: some View {
        ScrollView(.horizontal, showsIndicators: false) {
            HStack(spacing: 10) {
                ForEach(overlays, id: \.kind) { overlay in
                    Button {
                        if visible.contains(overlay.kind) {
                            visible.remove(overlay.kind)
                        } else {
                            visible.insert(overlay.kind)
                        }
                    } label: {
                        Text(label(for: overlay.kind))
                            .font(.caption2.weight(.semibold))
                            .padding(.horizontal, 10)
                            .padding(.vertical, 6)
                            .background(.ultraThinMaterial, in: Capsule(style: .continuous))
                            .overlay(
                                Capsule(style: .continuous)
                                    .strokeBorder(Color.white.opacity(visible.contains(overlay.kind) ? 0.35 : 0.12),
                                                  lineWidth: visible.contains(overlay.kind) ? 1 : 0.8)
                            )
                            .opacity(visible.contains(overlay.kind) ? 1.0 : 0.65)
                    }
                    .buttonStyle(.plain)
                }
            }
        }
    }

    @ViewBuilder
    private func overlayTrack(_ overlay: FilmAnalysis.Overlay) -> some View {
        let counterpoint = overlays.first(where: { $0.kind == "counterpoint" && $0.channels == 1 })
        VStack(alignment: .leading, spacing: 6) {
            Text(label(for: overlay.kind))
                .font(.caption)
                .foregroundStyle(.secondary)

            if overlay.kind == "hue_hist12", overlay.channels == 12, !overlay.multiValues.isEmpty {
                HueHistogramStrip(valuesByChannel: overlay.multiValues, duration: duration, hoveredSecond: $hoveredSecond)
                    .frame(height: 14)
            } else if overlay.kind == "chroma12", overlay.channels == 12, !overlay.multiValues.isEmpty {
                ChromaTrack(valuesByChannel: overlay.multiValues, duration: duration, hoveredSecond: $hoveredSecond)
                    .frame(height: 40)
            } else if overlay.kind == "align_norm", overlay.channels == 1 {
                AlignNormTrack(samples: overlay.samples,
                               counterpoint: counterpoint?.samples ?? [],
                               duration: duration,
                               hoveredSecond: $hoveredSecond)
                    .frame(height: 44)
            } else if overlay.kind == "key_pc", overlay.valueType == .i32, overlay.channels == 1 {
                KeyTrack(samples: overlay.samples, duration: duration, hoveredSecond: $hoveredSecond)
                    .frame(height: 18)
            } else if overlay.kind == "harmonic_tension", overlay.channels == 1 {
                QuantileBandLineTrack(samples: overlay.samples,
                                      duration: duration,
                                      q: 0.7,
                                      bandColor: .purple.opacity(0.18),
                                      strokeColor: .purple.opacity(0.75),
                                      hoveredSecond: $hoveredSecond)
                    .frame(height: 44)
            } else if overlay.kind == "tempo_confidence", overlay.channels == 1 {
                ConfidenceBandTrack(samples: overlay.samples, duration: duration, hoveredSecond: $hoveredSecond)
                    .frame(height: 18)
            } else if overlay.kind == "speech_probability", overlay.channels == 1 {
                SpeechBandTrack(samples: overlay.samples, duration: duration, hoveredSecond: $hoveredSecond)
                    .frame(height: 30)
            } else if overlay.kind == "dialogue_clarity", overlay.channels == 1 {
                SpeechBandTrack(samples: overlay.samples, duration: duration, hoveredSecond: $hoveredSecond)
                    .frame(height: 30)
            } else if overlay.kind == "dialogue_energy", overlay.channels == 1 {
                AutoScaledLineTrack(samples: overlay.samples, duration: duration, hoveredSecond: $hoveredSecond)
                    .frame(height: 38)
            } else if overlay.kind == "music_energy", overlay.channels == 1 {
                AutoScaledLineTrack(samples: overlay.samples, duration: duration, hoveredSecond: $hoveredSecond)
                    .frame(height: 38)
            } else if overlay.kind == "effects_energy", overlay.channels == 1 {
                AutoScaledLineTrack(samples: overlay.samples, duration: duration, hoveredSecond: $hoveredSecond)
                    .frame(height: 38)
            } else if overlay.kind == "camera_motion_type", overlay.valueType == .i32, overlay.channels == 1 {
                CameraMotionTypeTrack(samples: overlay.samples, duration: duration)
                    .frame(height: 14)
            } else if overlay.kind == "shot_scale", overlay.valueType == .i32, overlay.channels == 1 {
                ShotScaleTrack(samples: overlay.samples, duration: duration)
                    .frame(height: 14)
            } else if overlay.kind == "state_shift", overlay.valueType == .i32, overlay.channels == 1 {
                EventTickTrack(samples: overlay.samples, duration: duration)
                    .frame(height: 14)
            } else if overlay.kind == "counterpoint", overlay.valueType == .i32, overlay.channels == 1 {
                CounterpointBand(samples: overlay.samples, duration: duration)
                    .frame(height: 14)
            } else if overlay.kind == "frontal_score", overlay.channels == 1 {
                Raw01LineTrack(samples: overlay.samples, duration: duration, hoveredSecond: $hoveredSecond)
                    .frame(height: 38)
            } else if overlay.kind == "downcast_score", overlay.channels == 1 {
                Raw01LineTrack(samples: overlay.samples, duration: duration, hoveredSecond: $hoveredSecond)
                    .frame(height: 38)
            } else if overlay.kind == "stereo_width", overlay.channels == 1 {
                Raw01LineTrack(samples: overlay.samples, duration: duration, hoveredSecond: $hoveredSecond)
                    .frame(height: 38)
            } else if overlay.kind == "reverb_amount", overlay.channels == 1 {
                Raw01LineTrack(samples: overlay.samples, duration: duration, hoveredSecond: $hoveredSecond)
                    .frame(height: 38)
            } else if overlay.kind == "silence_score", overlay.channels == 1 {
                Raw01LineTrack(samples: overlay.samples, duration: duration, hoveredSecond: $hoveredSecond)
                    .frame(height: 38)
            } else if overlay.valueType == .i32, overlay.channels == 1 {
                EventBandTrack(samples: overlay.samples, duration: duration)
                    .frame(height: 14)
            } else if overlay.channels == 1 {
                NormalizedLineTrack(samples: overlay.samples, duration: duration, hoveredSecond: $hoveredSecond)
                    .frame(height: 38)
            } else {
                Text(verbatim: L10n.format("overlay.multichannel", overlay.channels))
                    .font(.caption2)
                    .foregroundStyle(.secondary)
            }
        }
        .padding(10)
        .background(Color.black.opacity(0.18), in: RoundedRectangle(cornerRadius: 12, style: .continuous))
        .overlay {
            if let selectedRange {
                SelectionSpanOverlay(duration: duration,
                                     start: selectedRange.start,
                                     end: selectedRange.end)
                    .clipShape(RoundedRectangle(cornerRadius: 12, style: .continuous))
            }
        }
    }

    private func label(for kind: String) -> String {
        let fallback: String = {
            switch kind {
        case "color_brightness": return "ColorBrightness"
        case "color_saturation": return "ColorSaturation"
        case "color_warmth": return "ColorWarmth"
        case "color_grayness": return "ColorGrayness"
        case "color_energy": return "ColorEnergy"
        case "color_harmony": return "ColorHarmony"
        case "color_contrast": return "ColorContrast (ΔE)"
        case "hue_deg": return "Hue"
        case "hue_hist12": return "Hue Track"
        case "color_state": return "ColorState"
        case "state_shift": return "StateShift"
        case "align_norm": return "AlignNorm"
        case "counterpoint": return "Counterpoint"
        case "loudness": return "Loudness"
        case "rhythm_drive": return "RhythmDrive"
        case "spectral_brightness": return "SpectralBrightness"
        case "spectral_centroid_hz": return "SpectralCentroid (Hz)"
        case "spectral_hf_ratio": return "HF Energy Ratio"
        case "dialogue_clarity": return "DialogueClarity"
        case "speech_probability": return "SpeechProbability"
        case "bpm": return "BPM"
        case "tempo_confidence": return "TempoConfidence"
        case "harmonic_tension": return "HarmonicTension"
        case "stereo_width": return "StereoWidth"
        case "reverb_amount": return "ReverbAmount"
        case "silence_score": return "SilenceScore"
        case "key_change_density": return "KeyChangeDensity"
        case "is_speaking": return "IsSpeaking"
        case "look_at_camera": return "LookAtCamera"
        case "frontal_score": return "FrontalScore"
        case "downcast_score": return "DowncastScore"
        case "dialogue_dominance": return "DialogueDominance"
        case "music_dominance": return "MusicDominance"
        case "effects_dominance": return "EffectsDominance"
        case "dialogue_energy": return "DialogueEnergy"
        case "music_energy": return "MusicEnergy"
        case "effects_energy": return "EffectsEnergy"
        case "chroma12": return "Chroma"
        case "key_pc": return "Key"
        case "face_presence": return "FacePresence"
        case "camera_motion": return "CameraMotion"
        case "object_motion": return "ObjectMotion"
        case "camera_vs_object_delta": return "Camera−Object"
        case "camera_motion_type": return "CameraMotionType"
        case "shot_scale": return "ShotScale (ELS→ECU)"
        case "exposition_curve": return "Exposition"
        case "verbal_load": return "VerbalLoad"
        case "visual_load": return "VisualLoad"
        case "event_load": return "EventLoad"
        case "visual_entropy": return "VisualEntropy"
        case "shot_scale_numeric": return "ShotScaleNumeric"
        case "camera_motion_complexity": return "CameraMotionComplexity"
        case "face_change": return "FaceChange"
        case "emotion_progress": return "EmotionProgress u(t)"
        default: return kind
            }
        }()

        let key = "overlay.kind.\(kind)"
        let localized = NSLocalizedString(key, bundle: .module, comment: "")
        return (localized == key) ? fallback : localized
    }
}

private struct Raw01LineTrack: View {
    let samples: [FilmAnalysis.Sample]
    let duration: Double
    @Binding var hoveredSecond: Double?

    var body: some View {
        return Chart {
            ForEach(samples, id: \.time) { sample in
                LineMark(
                    x: .value("Time", sample.time),
                    y: .value("Value", max(0.0, min(1.0, sample.value)))
                )
                .interpolationMethod(.catmullRom)
                .foregroundStyle(Color.white.opacity(0.65))
                .lineStyle(StrokeStyle(lineWidth: 1.2))
            }
        }
        .chartXScale(domain: 0...duration)
        .chartYScale(domain: 0...1)
        .chartXAxis(.hidden)
        .chartYAxis(.hidden)
        .overlay(alignment: .topTrailing) {
            if let hoveredSecond,
               let i = samples.lastIndex(where: { $0.time <= hoveredSecond }) {
                Text(verbatim: TimelineFormatting.decimal(max(0.0, min(1.0, samples[i].value)), maxFractionDigits: 2))
                    .font(.caption2.monospacedDigit())
                    .foregroundStyle(.secondary)
            }
        }
    }
}

private struct AutoScaledLineTrack: View {
    let samples: [FilmAnalysis.Sample]
    let duration: Double
    @Binding var hoveredSecond: Double?

    var body: some View {
        let values = samples.map(\.value).filter { $0.isFinite }
        let q05 = quantile(values, 0.05)
        let q95 = max(quantile(values, 0.95), q05 + 1e-12)
        let ymin = max(0.0, q05)
        let ymax = q95

        return Chart {
            ForEach(samples, id: \.time) { sample in
                LineMark(
                    x: .value("Time", sample.time),
                    y: .value("Value", sample.value)
                )
                .interpolationMethod(.catmullRom)
                .foregroundStyle(Color.white.opacity(0.65))
                .lineStyle(StrokeStyle(lineWidth: 1.2))
            }
        }
        .chartXScale(domain: 0...duration)
        .chartYScale(domain: ymin...ymax)
        .chartXAxis(.hidden)
        .chartYAxis(.hidden)
        .overlay(alignment: .topTrailing) {
            if let hoveredSecond,
               let i = samples.lastIndex(where: { $0.time <= hoveredSecond }) {
                Text(verbatim: TimelineFormatting.significant(samples[i].value, maxSignificantDigits: 3))
                    .font(.caption2.monospacedDigit())
                    .foregroundStyle(.secondary)
            }
        }
        .background(Color.white.opacity(0.06), in: RoundedRectangle(cornerRadius: 6, style: .continuous))
    }

    private func quantile(_ values: [Double], _ q: Double) -> Double {
        guard !values.isEmpty else { return 0.0 }
        let qq = min(1.0, max(0.0, q))
        let sorted = values.sorted()
        if sorted.count == 1 { return sorted[0] }
        let pos = qq * Double(sorted.count - 1)
        let i0 = Int(floor(pos))
        let i1 = Int(ceil(pos))
        if i0 == i1 { return sorted[i0] }
        let a0 = sorted[i0]
        let a1 = sorted[i1]
        return a0 + (pos - Double(i0)) * (a1 - a0)
    }
}

private struct EventBandTrack: View {
    let samples: [FilmAnalysis.Sample]
    let duration: Double

    var body: some View {
        GeometryReader { geo in
            Canvas { ctx, size in
                guard duration > 0 else { return }
                let count = Int(ceil(duration))
                let valuesBySecond = valuesPerSecond(duration: duration)
                for sec in 0..<count {
                    let v = valuesBySecond[min(sec, valuesBySecond.count - 1)]
                    guard v > 0.5 else { continue }
                    let x0 = CGFloat(Double(sec) / duration) * size.width
                    let x1 = CGFloat(Double(sec + 1) / duration) * size.width
                    let rect = CGRect(x: x0, y: 0, width: max(1, x1 - x0), height: size.height)
                    ctx.fill(Path(roundedRect: rect, cornerRadius: 2), with: .color(.white.opacity(0.65)))
                }
            }
        }
        .background(Color.white.opacity(0.06), in: RoundedRectangle(cornerRadius: 6, style: .continuous))
    }

    private func valuesPerSecond(duration: Double) -> [Double] {
        let count = Int(ceil(duration))
        guard count > 0 else { return [] }
        var out = Array(repeating: 0.0, count: count)
        for s in samples {
            let sec = max(0, min(count - 1, Int(floor(s.time))))
            out[sec] = max(out[sec], s.value)
        }
        return out
    }
}

private struct CameraMotionTypeTrack: View {
    let samples: [FilmAnalysis.Sample]
    let duration: Double

    var body: some View {
        GeometryReader { geo in
            Canvas { ctx, size in
                guard duration > 0 else { return }
                let count = Int(ceil(duration))
                let valuesBySecond = valuesPerSecond(duration: duration)
                for sec in 0..<count {
                    let v = Int(valuesBySecond[min(sec, valuesBySecond.count - 1)].rounded())
                    guard v != 0 else { continue }
                    let x0 = CGFloat(Double(sec) / duration) * size.width
                    let x1 = CGFloat(Double(sec + 1) / duration) * size.width
                    let rect = CGRect(x: x0, y: 0, width: max(1, x1 - x0), height: size.height)
                    ctx.fill(Path(roundedRect: rect, cornerRadius: 2), with: .color(color(for: v).opacity(0.65)))
                }
            }
        }
        .background(Color.white.opacity(0.06), in: RoundedRectangle(cornerRadius: 6, style: .continuous))
    }

    private func valuesPerSecond(duration: Double) -> [Double] {
        let count = Int(ceil(duration))
        guard count > 0 else { return [] }
        var out = Array(repeating: 0.0, count: count)
        for s in samples {
            let sec = max(0, min(count - 1, Int(floor(s.time))))
            out[sec] = s.value
        }
        return out
    }

    private func color(for v: Int) -> Color {
        switch v {
        case 1: return .cyan   // pan
        case 2: return .blue   // tilt
        case 3: return .orange // push
        case 4: return .purple // pull
        case 5: return .green  // track
        default: return .clear
        }
    }
}

private struct ShotScaleTrack: View {
    let samples: [FilmAnalysis.Sample]
    let duration: Double

    var body: some View {
        GeometryReader { _ in
            Canvas { ctx, size in
                guard duration > 0 else { return }
                let count = Int(ceil(duration))
                let valuesBySecond = valuesPerSecond(duration: duration)
                for sec in 0..<count {
                    let v = Int(valuesBySecond[min(sec, valuesBySecond.count - 1)].rounded())
                    guard v != 0 else { continue }
                    let x0 = CGFloat(Double(sec) / duration) * size.width
                    let x1 = CGFloat(Double(sec + 1) / duration) * size.width
                    let rect = CGRect(x: x0, y: 0, width: max(1, x1 - x0), height: size.height)
                    ctx.fill(Path(roundedRect: rect, cornerRadius: 2), with: .color(color(for: v).opacity(0.75)))
                }
            }
        }
        .background(Color.white.opacity(0.06), in: RoundedRectangle(cornerRadius: 6, style: .continuous))
    }

    private func valuesPerSecond(duration: Double) -> [Double] {
        let count = Int(ceil(duration))
        guard count > 0 else { return [] }
        var out = Array(repeating: 0.0, count: count)
        for s in samples {
            let sec = max(0, min(count - 1, Int(floor(s.time))))
            out[sec] = s.value
        }
        return out
    }

    private func color(for v: Int) -> Color {
        // 1..6: ELS → ECU
        switch v {
        case 1: return .indigo  // ELS
        case 2: return .blue    // LS
        case 3: return .cyan    // MS
        case 4: return .green   // MCU
        case 5: return .yellow  // CU
        case 6: return .orange  // ECU
        default: return .clear
        }
    }
}

private struct SelectionSpanOverlay: View {
    let duration: Double
    let start: Double
    let end: Double

    var body: some View {
        GeometryReader { geo in
            if duration > 0, end > start {
                let x0 = CGFloat(start / duration) * geo.size.width
                let x1 = CGFloat(end / duration) * geo.size.width
                let w = max(2, x1 - x0)
                RoundedRectangle(cornerRadius: 10, style: .continuous)
                    .fill(Color.white.opacity(0.06))
                    .frame(width: w, height: geo.size.height)
                    .offset(x: x0)
                    .overlay(
                        RoundedRectangle(cornerRadius: 10, style: .continuous)
                            .stroke(Color.white.opacity(0.12), lineWidth: 1)
                            .frame(width: w, height: geo.size.height)
                            .offset(x: x0)
                    )
            }
        }
        .allowsHitTesting(false)
    }
}

private struct NormalizedLineTrack: View {
    let samples: [FilmAnalysis.Sample]
    let duration: Double
    @Binding var hoveredSecond: Double?

    var body: some View {
        let xs = samples.map(\.value)
        let normalized = robustSigmoid01(xs)
        return Chart {
            ForEach(Array(samples.enumerated()), id: \.offset) { idx, sample in
                LineMark(
                    x: .value("Time", sample.time),
                    y: .value("Value", normalized[idx])
                )
                .interpolationMethod(.catmullRom)
                .foregroundStyle(Color.white.opacity(0.65))
                .lineStyle(StrokeStyle(lineWidth: 1.2))
            }
        }
        .chartXScale(domain: 0...duration)
        .chartYScale(domain: 0...1)
        .chartXAxis(.hidden)
        .chartYAxis(.hidden)
        .overlay(alignment: .topTrailing) {
            if let hoveredSecond,
               let i = samples.lastIndex(where: { $0.time <= hoveredSecond }) {
                Text(verbatim: TimelineFormatting.decimal(normalized[i], maxFractionDigits: 2))
                    .font(.caption2.monospacedDigit())
                    .foregroundStyle(.secondary)
            }
        }
    }

    private func robustSigmoid01(_ v: [Double]) -> [Double] {
        guard !v.isEmpty else { return [] }
        let med = median(v)
        let mad = median(v.map { abs($0 - med) })
        let denom = mad + 1e-6
        return v.map { x in
            let z = (x - med) / denom
            let clipped = max(-40.0, min(40.0, z))
            return 1.0 / (1.0 + exp(-clipped))
        }
    }

    private func median(_ v: [Double]) -> Double {
        guard !v.isEmpty else { return 0 }
        let s = v.sorted()
        if s.count % 2 == 1 { return s[s.count / 2] }
        return 0.5 * (s[s.count / 2 - 1] + s[s.count / 2])
    }
}

private struct DiscreteBandTrack: View {
    let samples: [FilmAnalysis.Sample]
    let duration: Double
    @Binding var hoveredSecond: Double?

    var body: some View {
        GeometryReader { geo in
            Canvas { ctx, size in
                guard !samples.isEmpty else { return }
                let seconds = Int(ceil(duration))
                for sec in 0..<seconds {
                    let t0 = Double(sec)
                    let x0 = CGFloat(t0 / max(duration, 1)) * size.width
                    let w = CGFloat(1.0 / max(duration, 1)) * size.width
                    let v = value(at: t0)
                    let color = colorForDiscrete(Int(v))
                    ctx.fill(Path(CGRect(x: x0, y: 0, width: w + 1, height: size.height)),
                             with: .color(color))
                }
            }
        }
    }

    private func value(at time: Double) -> Double {
        samples.last(where: { $0.time <= time })?.value ?? samples.first?.value ?? 0
    }

    private func colorForDiscrete(_ v: Int) -> Color {
        // deterministic palette
        let colors: [Color] = [.red, .orange, .yellow, .green, .mint, .cyan, .blue, .indigo, .purple, .pink, .brown, .teal]
        if v < 0 { return Color.white.opacity(0.08) }
        return colors[v % colors.count].opacity(0.55)
    }
}

private struct HueHistogramStrip: View {
    let valuesByChannel: [[Double]]   // [12][N]
    let duration: Double
    @Binding var hoveredSecond: Double?

    var body: some View {
        GeometryReader { geo in
            Canvas { ctx, size in
                let channels = valuesByChannel.count
                guard channels == 12, let count = valuesByChannel.first?.count, count > 0 else { return }
                for i in 0..<count {
                    var hist = [Double](repeating: 0.0, count: 12)
                    for c in 0..<12 {
                        hist[c] = max(0.0, valuesByChannel[c][i])
                    }
                    let hue = circularMeanDeg(hist) / 360.0
                    let sum = hist.reduce(0.0, +)
                    let purity = (sum > 1e-12) ? ((hist.max() ?? 0.0) / sum) : 0.0
                    // Strict semantics: saturation encodes "hue concentration" (purity), not a fixed threshold.
                    let sat = max(0.0, min(1.0, purity))
                    let color = Color(hue: hue, saturation: sat, brightness: 0.9).opacity(0.8)

                    let t0 = Double(i)
                    let x0 = CGFloat(t0 / max(duration, 1)) * size.width
                    let w = CGFloat(1.0 / max(duration, 1)) * size.width
                    ctx.fill(Path(CGRect(x: x0, y: 0, width: w + 1, height: size.height)), with: .color(color))
                }
            }
        }
    }

    private func circularMeanDeg(_ hist: [Double]) -> Double {
        var sx = 0.0, sy = 0.0, wsum = 0.0
        for b in 0..<12 {
            let w = hist[b]
            let ang = Double(b) * 30.0 + 15.0
            let rad = ang * Double.pi / 180.0
            sx += w * cos(rad)
            sy += w * sin(rad)
            wsum += w
        }
        if wsum <= 1e-12 { return 0 }
        var deg = atan2(sy, sx) * 180.0 / Double.pi
        if deg < 0 { deg += 360 }
        return deg
    }
}

private struct ChromaTrack: View {
    let valuesByChannel: [[Double]]   // [12][N]
    let duration: Double
    @Binding var hoveredSecond: Double?

    var body: some View {
        VStack(spacing: 4) {
            ChromaHueStrip(valuesByChannel: valuesByChannel, duration: duration, hoveredSecond: $hoveredSecond)
                .frame(height: 14)
            ChromaStackStrip(valuesByChannel: valuesByChannel, duration: duration)
                .frame(height: 22)
                .clipShape(RoundedRectangle(cornerRadius: 6, style: .continuous))
        }
    }
}

private struct ChromaHueStrip: View {
    let valuesByChannel: [[Double]]   // [12][N]
    let duration: Double
    @Binding var hoveredSecond: Double?

    var body: some View {
        GeometryReader { geo in
            ZStack(alignment: .topTrailing) {
                Canvas { ctx, size in
                    let channels = valuesByChannel.count
                    guard channels == 12, let count = valuesByChannel.first?.count, count > 0 else { return }
                    for i in 0..<count {
                        var hist = [Double](repeating: 0.0, count: 12)
                        for c in 0..<12 {
                            hist[c] = max(0.0, valuesByChannel[c][i])
                        }
                        let sum = hist.reduce(0.0, +)
                        let maxBin = hist.enumerated().max(by: { $0.element < $1.element })?.offset ?? 0
                        let purity = (sum > 1e-12) ? ((hist[maxBin]) / sum) : 0.0

                        // Strict strip semantics:
                        // - hue encodes pitch-class center (blend between circular mean and dominant bin)
                        // - saturation encodes concentration (purity=max/sum)
                        let meanHue = circularMeanHue01(hist)
                        let maxHue = Double(maxBin) / 12.0
                        let hue = meanHue * (1.0 - purity) + maxHue * purity
                        let sat = max(0.0, min(1.0, purity))
                        let color = Color(hue: hue, saturation: sat, brightness: 0.92).opacity(0.85)

                        let t0 = Double(i)
                        let x0 = CGFloat(t0 / max(duration, 1)) * size.width
                        let w = CGFloat(1.0 / max(duration, 1)) * size.width
                        ctx.fill(Path(CGRect(x: x0, y: 0, width: w + 1, height: size.height)), with: .color(color))
                    }
                }

                if let hoveredSecond {
                    hoveredBadge(at: hoveredSecond)
                        .padding(6)
                }
            }
        }
    }

    private func hoveredBadge(at t: Double) -> some View {
        let s = hoveredSummary(at: t)
        return Text(verbatim: L10n.format("overlay.badge.pitch_purity",
                                          pitchName(s.dominantPc),
                                          TimelineFormatting.decimal(s.purity, maxFractionDigits: 2)))
            .font(.caption2.monospacedDigit())
            .foregroundStyle(.secondary)
            .padding(.horizontal, 8)
            .padding(.vertical, 4)
            .background(.ultraThinMaterial, in: Capsule(style: .continuous))
    }

    private func hoveredSummary(at t: Double) -> (dominantPc: Int, purity: Double) {
        let count = valuesByChannel.first?.count ?? 0
        guard count > 0 else { return (0, 0.0) }
        let i = max(0, min(count - 1, Int(floor(t))))
        var hist = [Double](repeating: 0.0, count: 12)
        for c in 0..<12 {
            if c < valuesByChannel.count, i < valuesByChannel[c].count {
                hist[c] = max(0.0, valuesByChannel[c][i])
            }
        }
        let sum = hist.reduce(0.0, +)
        let maxBin = hist.enumerated().max(by: { $0.element < $1.element })?.offset ?? 0
        let purity = (sum > 1e-12) ? ((hist[maxBin]) / sum) : 0.0
        return (maxBin, purity)
    }

    private func pitchName(_ pc: Int) -> String {
        let names = ["C", "C#", "D", "D#", "E", "F", "F#", "G", "G#", "A", "A#", "B"]
        return names[pc % names.count]
    }

    private func circularMeanHue01(_ hist: [Double]) -> Double {
        var sx = 0.0, sy = 0.0, wsum = 0.0
        for b in 0..<12 {
            let w = hist[b]
            let ang = Double(b) * 30.0 + 15.0
            let rad = ang * Double.pi / 180.0
            sx += w * cos(rad)
            sy += w * sin(rad)
            wsum += w
        }
        if wsum <= 1e-12 { return 0.0 }
        var deg = atan2(sy, sx) * 180.0 / Double.pi
        if deg < 0 { deg += 360.0 }
        return deg / 360.0
    }
}

private struct ChromaStackStrip: View {
    let valuesByChannel: [[Double]]   // [12][N]
    let duration: Double

    var body: some View {
        GeometryReader { _ in
            Canvas { ctx, size in
                guard valuesByChannel.count == 12, let count = valuesByChannel.first?.count, count > 0 else { return }
                let rowH = size.height / 12.0
                for i in 0..<count {
                    let t0 = Double(i)
                    let x0 = CGFloat(t0 / max(duration, 1)) * size.width
                    let w = CGFloat(1.0 / max(duration, 1)) * size.width

                    var sum = 0.0
                    for pc in 0..<12 { sum += max(0.0, valuesByChannel[pc][i]) }
                    if sum <= 1e-12 { continue }

                    for pc in 0..<12 {
                        let v = max(0.0, valuesByChannel[pc][i]) / sum
                        guard v > 1e-4 else { continue }
                        let y0 = CGFloat(pc) * rowH
                        let r = CGRect(x: x0, y: y0, width: w + 1, height: rowH)
                        ctx.fill(Path(r), with: .color(colorForPitchClass(pc).opacity(0.10 + 0.85 * v)))
                    }
                }

                // Row separators (subtle, keep grid stable)
                for pc in 1..<12 {
                    let y = CGFloat(pc) * (size.height / 12.0)
                    var p = Path()
                    p.move(to: CGPoint(x: 0, y: y))
                    p.addLine(to: CGPoint(x: size.width, y: y))
                    ctx.stroke(p, with: .color(.white.opacity(0.05)), lineWidth: 1)
                }
            }
        }
        .background(Color.black.opacity(0.10))
    }

    private func colorForPitchClass(_ pc: Int) -> Color {
        // Deterministic 12-color wheel: C=red-ish, then clockwise.
        let hue = Double(pc % 12) / 12.0
        return Color(hue: hue, saturation: 0.85, brightness: 0.95)
    }
}

private struct KeyTrack: View {
    let samples: [FilmAnalysis.Sample] // value: 0..23 (major/minor), -1 unknown
    let duration: Double
    @Binding var hoveredSecond: Double?

    var body: some View {
        GeometryReader { _ in
            ZStack(alignment: .topTrailing) {
                Canvas { ctx, size in
                    guard !samples.isEmpty else { return }
                    let seconds = Int(ceil(duration))
                    for sec in 0..<seconds {
                        let t0 = Double(sec)
                        let v = Int(round(value(at: t0)))
                        let x0 = CGFloat(t0 / max(duration, 1)) * size.width
                        let w = CGFloat(1.0 / max(duration, 1)) * size.width
                        let r = CGRect(x: x0, y: 0, width: w + 1, height: size.height)
                        ctx.fill(Path(r), with: .color(colorForKey(v)))
                    }
                }

                if let hoveredSecond {
                    let v = Int(round(value(at: hoveredSecond)))
                    Text(keyName(v))
                        .font(.caption2.monospacedDigit())
                        .foregroundStyle(.secondary)
                        .padding(.horizontal, 8)
                        .padding(.vertical, 4)
                        .background(.ultraThinMaterial, in: Capsule(style: .continuous))
                        .padding(6)
                }
            }
        }
        .background(Color.black.opacity(0.10))
        .clipShape(RoundedRectangle(cornerRadius: 6, style: .continuous))
    }

    private func value(at time: Double) -> Double {
        samples.last(where: { $0.time <= time })?.value ?? samples.first?.value ?? -1
    }

    private func keyName(_ v: Int) -> String {
        if v < 0 { return "Key: —" }
        let names = ["C", "C#", "D", "D#", "E", "F", "F#", "G", "G#", "A", "A#", "B"]
        let pc = v % 12
        let isMinor = v >= 12
        return "Key: \(names[pc])\(isMinor ? "m" : "")"
    }

    private func colorForKey(_ v: Int) -> Color {
        if v < 0 { return Color.white.opacity(0.06) }
        let pc = v % 12
        let hue = Double(pc) / 12.0
        let isMinor = v >= 12
        return Color(hue: hue, saturation: isMinor ? 0.55 : 0.80, brightness: isMinor ? 0.70 : 0.88).opacity(0.85)
    }
}

private struct QuantileBandLineTrack: View {
    let samples: [FilmAnalysis.Sample]
    let duration: Double
    let q: Double
    let bandColor: Color
    let strokeColor: Color
    @Binding var hoveredSecond: Double?

    var body: some View {
        let series = resampleBySecond(samples: samples, duration: duration)
        let thr = quantile(series, q: q)
        let mask = debouncedMask(series.map { $0 > thr })

        return ZStack {
            SegmentBackground(mask: mask, duration: duration, fill: bandColor)
            ChartLine(samples: samples, duration: duration, color: strokeColor, hoveredSecond: $hoveredSecond)
        }
        .background(Color.black.opacity(0.08))
        .clipShape(RoundedRectangle(cornerRadius: 8, style: .continuous))
    }

    private func resampleBySecond(samples: [FilmAnalysis.Sample], duration: Double) -> [Double] {
        let seconds = max(1, Int(ceil(duration)))
        var out = Array(repeating: 0.0, count: seconds)
        var idx = 0
        for sec in 0..<seconds {
            let t = Double(sec) + 0.5
            while idx + 1 < samples.count && samples[idx + 1].time <= t { idx += 1 }
            out[sec] = (idx < samples.count) ? samples[idx].value : 0.0
        }
        return out
    }

    private func quantile(_ v: [Double], q: Double) -> Double {
        guard !v.isEmpty else { return 0.0 }
        let qq = max(0.0, min(1.0, q))
        let s = v.sorted()
        let k = Int(round(qq * Double(max(0, s.count - 1))))
        return s[max(0, min(s.count - 1, k))]
    }

    private func debouncedMask(_ mask: [Bool]) -> [Bool] {
        guard mask.count >= 3 else { return mask }
        var out = mask
        for i in 1..<(mask.count - 1) {
            if out[i - 1] == out[i + 1], out[i] != out[i - 1] { out[i] = out[i - 1] }
        }
        return out
    }
}

private struct SegmentBackground: View {
    let mask: [Bool] // per-second
    let duration: Double
    let fill: Color

    var body: some View {
        GeometryReader { _ in
            Canvas { ctx, size in
                guard !mask.isEmpty else { return }
                let seconds = mask.count
                var inSeg = false
                var segStart = 0

                func flush(end: Int) {
                    let x0 = CGFloat(Double(segStart) / max(duration, 1)) * size.width
                    let x1 = CGFloat(Double(end) / max(duration, 1)) * size.width
                    let r = CGRect(x: x0, y: 0, width: max(1, x1 - x0), height: size.height)
                    ctx.fill(Path(r), with: .color(fill))
                }

                for sec in 0..<seconds {
                    if mask[sec] {
                        if !inSeg { inSeg = true; segStart = sec }
                    } else {
                        if inSeg { flush(end: sec); inSeg = false }
                    }
                }
                if inSeg { flush(end: seconds) }
            }
        }
    }
}

private struct ChartLine: View {
    let samples: [FilmAnalysis.Sample]
    let duration: Double
    let color: Color
    @Binding var hoveredSecond: Double?

    var body: some View {
        Chart {
            ForEach(Array(samples.enumerated()), id: \.offset) { _, s in
                LineMark(x: .value("Time", s.time),
                         y: .value("Value", s.value))
                .interpolationMethod(.catmullRom)
                .foregroundStyle(color)
                .lineStyle(StrokeStyle(lineWidth: 1.4))
            }
        }
        .chartXScale(domain: 0...duration)
        .chartYScale(domain: 0...1)
        .chartXAxis(.hidden)
        .chartYAxis(.hidden)
        .overlay(alignment: .topTrailing) {
            if let hoveredSecond,
               let i = samples.lastIndex(where: { $0.time <= hoveredSecond }) {
                Text(verbatim: TimelineFormatting.decimal(samples[i].value, maxFractionDigits: 2))
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

private struct ConfidenceBandTrack: View {
    let samples: [FilmAnalysis.Sample] // expected [0,1]
    let duration: Double
    @Binding var hoveredSecond: Double?

    var body: some View {
        let series = resampleBySecond(samples: samples, duration: duration)
        let thr = quantile(series, q: 0.25)
        return ZStack {
            // Low-confidence band (adaptive): conf < Q(0.25)
            SegmentBackground(mask: series.map { $0 < thr }, duration: duration, fill: Color.white.opacity(0.08))
            GeometryReader { _ in
                Canvas { ctx, size in
                    let seconds = series.count
                    for sec in 0..<seconds {
                        let t0 = Double(sec)
                        let x0 = CGFloat(t0 / max(duration, 1)) * size.width
                        let w = CGFloat(1.0 / max(duration, 1)) * size.width
                        let v = max(0.0, min(1.0, series[sec]))
                        let col = Color.white.opacity(0.12 + 0.70 * v)
                        let r = CGRect(x: x0, y: size.height * (1.0 - CGFloat(v)), width: w + 1, height: size.height * CGFloat(v))
                        ctx.fill(Path(r), with: .color(col))
                    }
                }
            }
        }
        .background(Color.black.opacity(0.10))
        .clipShape(RoundedRectangle(cornerRadius: 6, style: .continuous))
        .overlay(alignment: .topTrailing) {
            if let hoveredSecond {
                let i = max(0, min(series.count - 1, Int(floor(hoveredSecond))))
                Text(verbatim: L10n.format("overlay.conf_value", TimelineFormatting.decimal(series[i], maxFractionDigits: 2)))
                    .font(.caption2.monospacedDigit())
                    .foregroundStyle(.secondary)
                    .padding(.horizontal, 8)
                    .padding(.vertical, 4)
                    .background(.ultraThinMaterial, in: Capsule(style: .continuous))
                    .padding(6)
            }
        }
    }

    private func resampleBySecond(samples: [FilmAnalysis.Sample], duration: Double) -> [Double] {
        let seconds = max(1, Int(ceil(duration)))
        var out = Array(repeating: 0.0, count: seconds)
        var idx = 0
        for sec in 0..<seconds {
            let t = Double(sec) + 0.5
            while idx + 1 < samples.count && samples[idx + 1].time <= t { idx += 1 }
            out[sec] = (idx < samples.count) ? samples[idx].value : 0.0
        }
        return out
    }

    private func quantile(_ v: [Double], q: Double) -> Double {
        guard !v.isEmpty else { return 0.0 }
        let qq = max(0.0, min(1.0, q))
        let s = v.sorted()
        let k = Int(round(qq * Double(max(0, s.count - 1))))
        return s[max(0, min(s.count - 1, k))]
    }
}

private struct SpeechBandTrack: View {
    let samples: [FilmAnalysis.Sample] // expected [0,1]
    let duration: Double
    @Binding var hoveredSecond: Double?

    var body: some View {
        let series = resampleBySecond(samples: samples, duration: duration)
        let hi = quantile(series, q: 0.75)
        let mask = debouncedMask(series.map { $0 > hi })

        return ZStack {
            // Paragraph-style speech band (adaptive): prob > Q(0.75), debounced.
            SegmentBackground(mask: mask, duration: duration, fill: Color.cyan.opacity(0.14))
            GeometryReader { _ in
                Canvas { ctx, size in
                    let seconds = series.count
                    for sec in 0..<seconds {
                        let t0 = Double(sec)
                        let x0 = CGFloat(t0 / max(duration, 1)) * size.width
                        let w = CGFloat(1.0 / max(duration, 1)) * size.width
                        let v = max(0.0, min(1.0, series[sec]))
                        let y = size.height * (1.0 - CGFloat(v))
                        let r = CGRect(x: x0, y: y, width: w + 1, height: max(1, size.height * CGFloat(v)))
                        ctx.fill(Path(r), with: .color(Color.cyan.opacity(0.10 + 0.60 * v)))
                    }
                }
            }
        }
        .background(Color.black.opacity(0.10))
        .clipShape(RoundedRectangle(cornerRadius: 6, style: .continuous))
        .overlay(alignment: .topTrailing) {
            if let hoveredSecond {
                let i = max(0, min(series.count - 1, Int(floor(hoveredSecond))))
                let stateKey = mask[i] ? "overlay.state.on" : "overlay.state.off"
                let state = L10n.string(stateKey)
                Text(verbatim: L10n.format("overlay.state_value", state, TimelineFormatting.decimal(series[i], maxFractionDigits: 2)))
                    .font(.caption2.monospacedDigit())
                    .foregroundStyle(.secondary)
                    .padding(.horizontal, 8)
                    .padding(.vertical, 4)
                    .background(.ultraThinMaterial, in: Capsule(style: .continuous))
                    .padding(6)
            }
        }
    }

    private func resampleBySecond(samples: [FilmAnalysis.Sample], duration: Double) -> [Double] {
        let seconds = max(1, Int(ceil(duration)))
        var out = Array(repeating: 0.0, count: seconds)
        var idx = 0
        for sec in 0..<seconds {
            let t = Double(sec) + 0.5
            while idx + 1 < samples.count && samples[idx + 1].time <= t { idx += 1 }
            out[sec] = (idx < samples.count) ? samples[idx].value : 0.0
        }
        return out
    }

    private func quantile(_ v: [Double], q: Double) -> Double {
        guard !v.isEmpty else { return 0.0 }
        let qq = max(0.0, min(1.0, q))
        let s = v.sorted()
        let k = Int(round(qq * Double(max(0, s.count - 1))))
        return s[max(0, min(s.count - 1, k))]
    }

    private func debouncedMask(_ mask: [Bool]) -> [Bool] {
        guard mask.count >= 3 else { return mask }
        var out = mask
        for i in 1..<(mask.count - 1) {
            if out[i - 1] == out[i + 1], out[i] != out[i - 1] { out[i] = out[i - 1] }
        }
        return out
    }
}

private struct EventTickTrack: View {
    let samples: [FilmAnalysis.Sample] // value: 0/1
    let duration: Double

    var body: some View {
        GeometryReader { _ in
            Canvas { ctx, size in
                guard !samples.isEmpty else { return }
                let seconds = Int(ceil(duration))
                for sec in 0..<seconds {
                    let t0 = Double(sec)
                    let v = value(at: t0)
                    guard v > 0.5 else { continue }
                    let x = CGFloat(t0 / max(duration, 1)) * size.width
                    let r = CGRect(x: x, y: 0, width: 2.0, height: size.height)
                    ctx.fill(Path(r), with: .color(.white.opacity(0.65)))
                }
            }
        }
    }

    private func value(at time: Double) -> Double {
        samples.last(where: { $0.time <= time })?.value ?? 0
    }
}

private struct CounterpointBand: View {
    let samples: [FilmAnalysis.Sample] // value: 0/1
    let duration: Double

    var body: some View {
        GeometryReader { _ in
            Canvas { ctx, size in
                guard !samples.isEmpty else { return }
                let seconds = Int(ceil(duration))
                for sec in 0..<seconds {
                    let t0 = Double(sec)
                    let v = value(at: t0)
                    guard v > 0.5 else { continue }
                    let x0 = CGFloat(t0 / max(duration, 1)) * size.width
                    let w = CGFloat(1.0 / max(duration, 1)) * size.width
                    let r = CGRect(x: x0, y: 0, width: w + 1, height: size.height)
                    ctx.fill(Path(r), with: .color(.red.opacity(0.35)))
                }
            }
        }
        .background(Color.black.opacity(0.08))
        .clipShape(RoundedRectangle(cornerRadius: 6, style: .continuous))
    }

    private func value(at time: Double) -> Double {
        samples.last(where: { $0.time <= time })?.value ?? 0
    }
}

private struct AlignNormTrack: View {
    let samples: [FilmAnalysis.Sample]
    let counterpoint: [FilmAnalysis.Sample] // 0/1
    let duration: Double
    @Binding var hoveredSecond: Double?

    var body: some View {
        ZStack {
            CounterpointBackground(samples: counterpoint, duration: duration)
            NormalizedLineTrack(samples: samples, duration: duration, hoveredSecond: $hoveredSecond)
        }
        .clipShape(RoundedRectangle(cornerRadius: 8, style: .continuous))
    }
}

private struct CounterpointBackground: View {
    let samples: [FilmAnalysis.Sample] // value: 0/1
    let duration: Double

    var body: some View {
        GeometryReader { _ in
            Canvas { ctx, size in
                guard !samples.isEmpty else { return }
                let seconds = Int(ceil(duration))
                var inSeg = false
                var segStart: Int = 0

                func flush(at end: Int) {
                    let t0 = Double(segStart)
                    let t1 = Double(end)
                    let x0 = CGFloat(t0 / max(duration, 1)) * size.width
                    let x1 = CGFloat(t1 / max(duration, 1)) * size.width
                    let r = CGRect(x: x0, y: 0, width: max(1, x1 - x0), height: size.height)
                    ctx.fill(Path(r), with: .color(.red.opacity(0.16)))
                    ctx.stroke(Path(r.insetBy(dx: 0.5, dy: 0.5)), with: .color(.red.opacity(0.28)), lineWidth: 1)
                }

                for sec in 0..<seconds {
                    let t0 = Double(sec)
                    let v = value(at: t0)
                    if v > 0.5 {
                        if !inSeg { inSeg = true; segStart = sec }
                    } else {
                        if inSeg { flush(at: sec); inSeg = false }
                    }
                }
                if inSeg { flush(at: seconds) }
            }
        }
        .background(Color.black.opacity(0.08))
    }

    private func value(at time: Double) -> Double {
        samples.last(where: { $0.time <= time })?.value ?? 0
    }
}
