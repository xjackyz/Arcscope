import SwiftUI

struct ExportCenterView: View {
    let analysis: FilmAnalysis

    @State private var mode: ExportMode = .poster

    var body: some View {
        VStack(spacing: 10) {
            Picker("export.mode", selection: $mode) {
                Text("export.mode.poster").tag(ExportMode.poster)
                Text("export.mode.report").tag(ExportMode.report)
            }
            .pickerStyle(.segmented)
            .frame(maxWidth: 520)

            switch mode {
            case .poster:
                PosterExportView(analysis: analysis)
            case .report:
                ReportExportView(analysis: analysis)
            }
        }
    }

    private enum ExportMode: Hashable {
        case poster
        case report
    }
}
