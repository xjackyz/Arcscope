import Foundation

protocol ArcScopeBackend {
    /// 加载分析结果。databaseURL 为空时会尝试默认路径。
    func loadFilm(databaseURL: URL?, filmId: String?) async throws -> FilmAnalysis
}

enum ArcScopeBackendError: Error, LocalizedError {
    case fileMissing
    case decodeFailed

    var errorDescription: String? {
        switch self {
        case .fileMissing:
            return L10n.string("error.db_missing")
        case .decodeFailed:
            return L10n.string("error.decode_failed")
        }
    }
}
