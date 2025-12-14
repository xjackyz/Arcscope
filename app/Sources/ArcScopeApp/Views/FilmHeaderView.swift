import SwiftUI

struct FilmHeaderView: View {
    let analysis: FilmAnalysis

    var body: some View {
        VStack(alignment: .leading, spacing: 10) {
            HStack(alignment: .center, spacing: 16) {
                VStack(alignment: .leading, spacing: 4) {
                    Text(analysis.title)
                        .font(.largeTitle.weight(.semibold))
                    Text(subtitle)
                        .font(.callout)
                        .foregroundStyle(.secondary)
                    if let path = analysis.filePath, !path.isEmpty {
                        Text(path)
                            .font(.caption.monospaced())
                            .foregroundStyle(.secondary)
                            .lineLimit(1)
                    }
                }
                Spacer()
            }

            if !analysis.curves.contains(where: { $0.kind == .faceAffect }) {
                Text("film.header.no_face")
                    .font(.caption)
                    .foregroundStyle(.secondary)
            }
        }
    }

    private var subtitle: String {
        let meta = TimelineFormatting.filmMetaLine(year: analysis.year,
                                                   director: analysis.director,
                                                   durationSeconds: analysis.duration)
        return L10n.format("film.header.summary", meta, analysis.curves.count, analysis.overlays.count)
    }
}
