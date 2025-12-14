import CAArcScopeBridge
import Foundation

final class ReferenceLibraryService {
    static let shared = ReferenceLibraryService()

    func loadReferences(databaseURL: URL) async throws -> [FilmSummary] {
        try await SecurityScopedAccess.withAccessAsync(databaseURL) {
            try await Task.detached { () -> [FilmSummary] in
            var list = ArcScopeFilmList()
            let status = databaseURL.path.withCString { dbPath in
                arcscope_list_reference_films(dbPath, &list)
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

    func setReference(databaseURL: URL, filmId: String, isReference: Bool, note: String? = nil) async throws {
        try await SecurityScopedAccess.withAccessAsync(databaseURL) {
            try await Task.detached {
                let status = databaseURL.path.withCString { dbPath in
                    filmId.withCString { filmCString in
                        note?.withCString { noteCString in
                            arcscope_set_reference_film(dbPath, filmCString, isReference ? 1 : 0, noteCString)
                        } ?? arcscope_set_reference_film(dbPath, filmCString, isReference ? 1 : 0, nil)
                    }
                }
                guard status == 0 else {
                    throw ArcScopeBackendError.decodeFailed
                }
            }.value
        }
    }
}
