import Foundation

@MainActor
final class ArcScopeViewModel: ObservableObject {
    @Published var analysis: FilmAnalysis?
    @Published var referenceAnalysis: FilmAnalysis?
    @Published var films: [FilmSummary] = []
    @Published var referenceLibrary: [FilmSummary] = []
    @Published var isAnalyzing = false
    @Published var selectedFilmId: String?
    @Published var referenceFilmId: String?
    @Published var errorMessage: String?
    @Published var analysisProgress: Double?

    private let backend: ArcScopeBackend
    private let analysisService = VideoAnalysisService.shared
    private let libraryService = FilmLibraryService.shared
    private let referenceService = ReferenceLibraryService.shared
    private let bookmarkStore = SecurityScopedBookmarkStore.shared

    init(backend: ArcScopeBackend = BridgeBackend()) {
        self.backend = backend
    }

    private var defaultDatabaseURL: URL {
        AppPaths.defaultDatabaseURL()
    }

    @Published var selectedDatabaseURL: URL?

    func loadDefaultDatabase() {
        if let restored = bookmarkStore.restoreSelectedDatabase() {
            selectDatabase(url: restored)
            return
        }

        let defaultURL = defaultDatabaseURL
        if FileManager.default.fileExists(atPath: defaultURL.path) {
            selectDatabase(url: defaultURL)
            return
        }

        // First launch or no database yet: keep UI ready for user-driven import/analysis.
        selectedDatabaseURL = nil
        films = []
        analysis = nil
        referenceLibrary = []
        referenceAnalysis = nil
        errorMessage = nil
    }

    func selectDatabase(url: URL) {
        selectedDatabaseURL = url
        bookmarkStore.saveSelectedDatabase(url)
        Task {
            await loadLibrary(url: url, prefer: nil)
            await runAnalysis(url: url, filmId: selectedFilmId)
            await loadReferenceFilm(id: referenceFilmId, url: url)
        }
    }

    func selectFilm(id: String) {
        guard id != selectedFilmId else { return }
        selectedFilmId = id
        guard let dbURL = selectedDatabaseURL else { return }
        Task { await runAnalysis(url: dbURL, filmId: id) }
    }

    func selectReferenceFilm(id: String?) {
        referenceFilmId = id
        guard let dbURL = selectedDatabaseURL else { return }
        Task { await loadReferenceFilm(id: id, url: dbURL) }
    }

    func loadFilm(at url: URL?) {
        guard let dbURL = url ?? selectedDatabaseURL else { return }
        Task { await runAnalysis(url: dbURL, filmId: selectedFilmId) }
    }

    func analyzeVideo(at url: URL) {
        let destination = selectedDatabaseURL ?? defaultDatabaseURL
        selectedDatabaseURL = destination
        analysisProgress = 0.0
        Task { [weak self] in
            guard let self else { return }
            do {
                self.bookmarkStore.saveVideoFile(url)
                try await analysisService.analyzeVideo(inputURL: url,
                                                        databaseURL: destination,
                                                        title: self.deriveTitle(from: url),
                                                        year: 0,
                                                        director: "Unknown") { progress in
                    Task { @MainActor in
                        self.analysisProgress = progress
                    }
                }
                await self.loadLibrary(url: destination, prefer: nil)
                await self.runAnalysis(url: destination, filmId: self.selectedFilmId)
            } catch {
                self.analysis = nil
                self.errorMessage = error.localizedDescription
            }
            self.analysisProgress = nil
        }
    }

    private func loadLibrary(url: URL, prefer preferredFilmId: String?) async {
        do {
            let summaries = try await libraryService.loadFilms(databaseURL: url)
            films = summaries
            if let preferred = preferredFilmId, summaries.contains(where: { $0.id == preferred }) {
                selectedFilmId = preferred
            } else if let current = selectedFilmId, summaries.contains(where: { $0.id == current }) {
                // keep current selection
            } else {
                selectedFilmId = summaries.first?.id
            }
            referenceLibrary = (try? await referenceService.loadReferences(databaseURL: url)) ?? []

            if let ref = referenceFilmId,
               !referenceLibrary.contains(where: { $0.id == ref }) {
                referenceFilmId = nil
                referenceAnalysis = nil
            }
        } catch {
            errorMessage = error.localizedDescription
        }
    }

    func toggleReferenceLibraryFilm(id: String, isReference: Bool, note: String? = nil) {
        guard let dbURL = selectedDatabaseURL else { return }
        Task {
            do {
                try await referenceService.setReference(databaseURL: dbURL, filmId: id, isReference: isReference, note: note)
                referenceLibrary = (try? await referenceService.loadReferences(databaseURL: dbURL)) ?? []
                if !isReference, referenceFilmId == id {
                    referenceFilmId = nil
                    referenceAnalysis = nil
                }
            } catch {
                errorMessage = error.localizedDescription
            }
        }
    }

    private func runAnalysis(url: URL?, filmId: String?) async {
        guard let dbURL = url ?? selectedDatabaseURL else { return }
        isAnalyzing = true
        defer { isAnalyzing = false }
        do {
            let result = try await backend.loadFilm(databaseURL: dbURL, filmId: filmId ?? selectedFilmId)
            analysis = result
            errorMessage = nil
        } catch {
            analysis = nil
            errorMessage = error.localizedDescription
        }
    }

    private func loadReferenceFilm(id: String?, url: URL?) async {
        guard let id = id, let dbURL = url ?? selectedDatabaseURL else {
            referenceAnalysis = nil
            return
        }
        do {
            referenceAnalysis = try await backend.loadFilm(databaseURL: dbURL, filmId: id)
        } catch {
            referenceAnalysis = nil
        }
    }


    private func deriveTitle(from url: URL) -> String {
        url.deletingPathExtension().lastPathComponent
    }
}
