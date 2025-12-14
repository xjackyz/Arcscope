import CAArcScopeBridge
import Foundation

struct FilmSummary: Identifiable, Equatable {
    let id: String
    let title: String
    let duration: Double
    let year: Int
    let director: String
}

final class FilmLibraryService {
    static let shared = FilmLibraryService()

    func loadFilms(databaseURL: URL) async throws -> [FilmSummary] {
        try await SecurityScopedAccess.withAccessAsync(databaseURL) {
            try await Task.detached { () -> [FilmSummary] in
            var list = ArcScopeFilmList()
            let status = databaseURL.path.withCString { dbPath in
                arcscope_list_films(dbPath, &list)
            }
            guard status == 0 else {
                throw ArcScopeBackendError.decodeFailed
            }
            defer { arcscope_bridge_free_film_list(&list) }
            var results: [FilmSummary] = []
            if let items = list.items {
                for idx in 0..<list.count {
                    let item = items[idx]
                    let filmId = item.film_id.map { String(cString: $0) } ?? ""
                    let title = item.title.map { String(cString: $0) } ?? filmId
                    let director = item.director.map { String(cString: $0) } ?? ""
                    results.append(FilmSummary(id: filmId,
                                               title: title,
                                               duration: item.duration_seconds,
                                               year: Int(item.year),
                                               director: director))
                }
            }
            return results
            }.value
        }
    }
}
