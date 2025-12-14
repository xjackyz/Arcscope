import SwiftUI

struct ReferencePickerView: View {
    @EnvironmentObject private var viewModel: ArcScopeViewModel
    @Environment(\.dismiss) private var dismiss

    @State private var searchText: String = ""

    var body: some View {
        VStack(spacing: 14) {
            header
            Picker("reference.current", selection: Binding<String>(
                get: { viewModel.referenceFilmId ?? "" },
                set: { viewModel.selectReferenceFilm(id: $0.isEmpty ? nil : $0) }
            )) {
                Text("reference.none").tag("")
                ForEach(filteredReferenceLibrary) { film in
                    Text(film.title).tag(film.id)
                }
            }
            .pickerStyle(.menu)
            .frame(maxWidth: 520, alignment: .leading)

            List {
                Section("reference.section.library") {
                    ForEach(filteredFilms) { film in
                        ReferenceRow(film: film,
                                     isReference: isReference(film.id),
                                     isSelected: viewModel.referenceFilmId == film.id) { newValue in
                            viewModel.toggleReferenceLibraryFilm(id: film.id, isReference: newValue)
                        } select: {
                            viewModel.selectReferenceFilm(id: film.id)
                        }
                    }
                }
            }
            .listStyle(.inset)
        }
        .padding(18)
        .frame(minWidth: 720, minHeight: 520)
        .searchable(text: $searchText, placement: .toolbar)
    }

    private var header: some View {
        HStack(spacing: 12) {
            VStack(alignment: .leading, spacing: 4) {
                Text("reference.title")
                    .font(.title3.weight(.semibold))
                Text("reference.subtitle")
                    .font(.callout)
                    .foregroundStyle(.secondary)
            }
            Spacer()
            Button("reference.done") { dismiss() }
                .buttonStyle(.borderedProminent)
        }
    }

    private func isReference(_ id: String) -> Bool {
        viewModel.referenceLibrary.contains(where: { $0.id == id })
    }

    private var filteredReferenceLibrary: [FilmSummary] {
        let base = viewModel.referenceLibrary
        guard !searchText.isEmpty else { return base }
        return base.filter { $0.title.localizedCaseInsensitiveContains(searchText) || $0.id.localizedCaseInsensitiveContains(searchText) }
    }

    private var filteredFilms: [FilmSummary] {
        let base = viewModel.films
        guard !searchText.isEmpty else { return base }
        return base.filter { $0.title.localizedCaseInsensitiveContains(searchText) || $0.id.localizedCaseInsensitiveContains(searchText) }
    }
}

private struct ReferenceRow: View {
    let film: FilmSummary
    let isReference: Bool
    let isSelected: Bool
    let toggle: (Bool) -> Void
    let select: () -> Void

    var body: some View {
        HStack(spacing: 10) {
            Button {
                toggle(!isReference)
            } label: {
                Image(systemName: isReference ? "star.fill" : "star")
                    .foregroundStyle(isReference ? .yellow : .secondary)
            }
            .buttonStyle(.plain)

            VStack(alignment: .leading, spacing: 2) {
                Text(film.title)
                    .font(.callout.weight(.semibold))
                Text(verbatim: "\(film.year)  ·  \(TimelineFormatting.durationAbbreviated(film.duration))  ·  \(film.director)")
                    .font(.caption)
                    .foregroundStyle(.secondary)
                    .lineLimit(1)
            }
            Spacer()
            if isSelected {
                Text("reference.current")
                    .font(.caption2.weight(.semibold))
                    .padding(.horizontal, 8)
                    .padding(.vertical, 4)
                    .background(.ultraThinMaterial, in: Capsule(style: .continuous))
            }
            Button("reference.use") { select() }
                .buttonStyle(.bordered)
                .disabled(!isReference)
        }
        .contentShape(Rectangle())
    }
}
