import SwiftUI
import UniformTypeIdentifiers
#if canImport(AppKit)
import AppKit
#endif

struct LibraryView: View {
    @EnvironmentObject private var viewModel: ArcScopeViewModel

    @State private var showDatabaseImporter = false
    @State private var showVideoImporter = false
    @State private var showReferencePicker = false
    @State private var showExport = false
    @State private var didAppear = false

    var body: some View {
        ZStack {
            #if canImport(AppKit)
            Color(nsColor: .windowBackgroundColor).ignoresSafeArea()
            #else
            Color(.systemBackground).ignoresSafeArea()
            #endif

            NavigationSplitView {
                sidebar
            } detail: {
                detail
            }
            .navigationSplitViewStyle(.balanced)

            if let progress = viewModel.analysisProgress {
                ProgressHUD(progress: progress)
            }
        }
        .toolbar {
            ToolbarItemGroup(placement: .primaryAction) {
                Button {
                    showVideoImporter = true
                } label: {
                    Label("library.analyze_new_film", systemImage: "film")
                }

                Button {
                    showDatabaseImporter = true
                } label: {
                    Label("library.database", systemImage: "tray")
                }

                Button {
                    showReferencePicker = true
                } label: {
                    Label("library.references", systemImage: "star")
                }

                Button {
                    showExport = true
                } label: {
                    Label("library.export", systemImage: "square.and.arrow.up")
                }
                .disabled(viewModel.analysis == nil)
            }
        }
        .fileImporter(isPresented: $showDatabaseImporter,
                      allowedContentTypes: [UTType(filenameExtension: "db") ?? .data],
                      allowsMultipleSelection: false) { result in
            switch result {
            case .success(let urls):
                if let url = urls.first { viewModel.selectDatabase(url: url) }
            case .failure(let error):
                viewModel.errorMessage = error.localizedDescription
            }
        }
        .fileImporter(isPresented: $showVideoImporter,
                      allowedContentTypes: [.movie],
                      allowsMultipleSelection: false) { result in
            switch result {
            case .success(let urls):
                if let url = urls.first { viewModel.analyzeVideo(at: url) }
            case .failure(let error):
                viewModel.errorMessage = error.localizedDescription
            }
        }
        .sheet(isPresented: $showReferencePicker) {
            ReferencePickerView()
                .environmentObject(viewModel)
        }
        .sheet(isPresented: $showExport) {
            if let analysis = viewModel.analysis {
                ExportCenterView(analysis: analysis)
            } else {
                EmptyView()
            }
        }
        .focusedValue(\.libraryActions, LibraryActions(
            openDatabase: { showDatabaseImporter = true },
            analyzeNewFilm: { showVideoImporter = true },
            openReferences: { showReferencePicker = true },
            openExport: { showExport = true },
            canExport: viewModel.analysis != nil
        ))
        .onAppear {
            guard !didAppear else { return }
            didAppear = true
            viewModel.loadDefaultDatabase()
        }
    }

    private var sidebar: some View {
        VStack(alignment: .leading, spacing: 8) {
            HStack {
                Text("library.title")
                    .font(.title3.weight(.semibold))
                Spacer()
            }
            .padding(.horizontal, 12)
            .padding(.top, 10)

            if let db = viewModel.selectedDatabaseURL?.lastPathComponent {
                Text(verbatim: L10n.format("library.db_label", db))
                    .font(.caption)
                    .foregroundStyle(.secondary)
                    .padding(.horizontal, 12)
            }

            if let message = viewModel.errorMessage, viewModel.films.isEmpty {
                Text(message)
                    .font(.callout)
                    .foregroundStyle(.secondary)
                    .padding(.horizontal, 12)
            }

            List(selection: Binding<String?>(
                get: { viewModel.selectedFilmId },
                set: { newValue in
                    if let id = newValue { viewModel.selectFilm(id: id) }
                }
            )) {
                ForEach(viewModel.films) { film in
                    FilmCardView(film: film)
                        .tag(Optional(film.id))
                }
            }
            .listStyle(.sidebar)
        }
    }

    @ViewBuilder
    private var detail: some View {
        if let analysis = viewModel.analysis {
            FilmDetailView(analysis: analysis,
                           referenceAnalysis: viewModel.referenceAnalysis)
        } else if let message = viewModel.errorMessage {
            VStack(spacing: 12) {
                Text(message)
                Button("library.select_database") { showDatabaseImporter = true }
                    .buttonStyle(.borderedProminent)
            }
            .frame(maxWidth: .infinity, maxHeight: .infinity)
        } else if viewModel.selectedDatabaseURL == nil {
            VStack(spacing: 12) {
                Text("library.empty_state.title")
                    .font(.title3.weight(.semibold))
                Text("library.empty_state.body")
                    .foregroundStyle(.secondary)
                HStack(spacing: 10) {
                    Button("library.select_database") { showDatabaseImporter = true }
                        .buttonStyle(.borderedProminent)
                    Button("library.analyze_new_film") { showVideoImporter = true }
                        .buttonStyle(.bordered)
                }
            }
            .padding(24)
            .frame(maxWidth: .infinity, maxHeight: .infinity)
        } else {
            ProgressView("library.loading")
                .progressViewStyle(.circular)
                .frame(maxWidth: .infinity, maxHeight: .infinity)
        }
    }
}

private struct FilmCardView: View {
    let film: FilmSummary

    var body: some View {
        VStack(alignment: .leading, spacing: 4) {
            Text(film.title)
                .font(.callout.weight(.semibold))
                .lineLimit(1)
            Text(verbatim: "\(film.year) · \(TimelineFormatting.durationAbbreviated(film.duration)) · \(film.director)")
                .font(.caption)
                .foregroundStyle(.secondary)
                .lineLimit(1)
        }
        .padding(.vertical, 6)
    }
}

private struct ProgressHUD: View {
    let progress: Double
    var body: some View {
        VStack(spacing: 8) {
            ProgressView(value: progress)
                .progressViewStyle(.linear)
            Text(verbatim: L10n.format("progress.analyzing", TimelineFormatting.percent01(progress)))
                .font(.caption)
        }
        .padding(16)
        .background(.ultraThinMaterial, in: RoundedRectangle(cornerRadius: 16, style: .continuous))
        .shadow(radius: 10)
    }
}
