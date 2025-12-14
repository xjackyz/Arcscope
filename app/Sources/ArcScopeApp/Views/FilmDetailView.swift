import SwiftUI

struct FilmDetailView: View {
    let analysis: FilmAnalysis
    let referenceAnalysis: FilmAnalysis?

    @State private var hoveredSecond: Double?
    @State private var selectedRange: TimeRangeSelection?

    var body: some View {
        ScrollView {
            VStack(alignment: .leading, spacing: 16) {
                FilmHeaderView(analysis: analysis)
                CurveTimelineView(analysis: analysis,
                                  referenceAnalysis: referenceAnalysis,
                                  hoveredSecond: $hoveredSecond,
                                  selectedRange: $selectedRange)
                    .frame(height: 420)
                MetadataOverlaysView(overlays: analysis.overlays,
                                     duration: analysis.duration,
                                     hoveredSecond: $hoveredSecond,
                                     selectedRange: $selectedRange)
                DiagnosticsHeatmap(analysis: analysis, hoveredSecond: $hoveredSecond, selectedRange: $selectedRange)
                    .frame(height: 80)
                if !(analysis.shots.isEmpty && analysis.scenes.isEmpty && analysis.sequences.isEmpty && analysis.acts.isEmpty) {
                    StructureBandsView(analysis: analysis, hoveredSecond: $hoveredSecond, selectedRange: $selectedRange)
                }
                IssueListView(issues: analysis.issues)
                Spacer(minLength: 12)
            }
            .padding(24)
        }
    }
}

