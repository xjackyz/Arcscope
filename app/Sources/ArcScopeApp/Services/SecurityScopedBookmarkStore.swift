import Foundation

final class SecurityScopedBookmarkStore {
    static let shared = SecurityScopedBookmarkStore()

    private let selectedDatabaseKey = "bookmark.database.selected"
    private let videoBookmarksKey = "bookmark.video.byPath"

    private init() {}

    func saveSelectedDatabase(_ url: URL) {
        do {
            let data = try url.bookmarkData(options: [.withSecurityScope],
                                            includingResourceValuesForKeys: nil,
                                            relativeTo: nil)
            UserDefaults.standard.set(data, forKey: selectedDatabaseKey)
        } catch {
            // Best-effort; non-sandboxed apps may not require bookmarks.
        }
    }

    func restoreSelectedDatabase() -> URL? {
        guard let data = UserDefaults.standard.data(forKey: selectedDatabaseKey) else { return nil }
        var stale = false
        guard let url = try? URL(resolvingBookmarkData: data,
                                 options: [.withSecurityScope],
                                 relativeTo: nil,
                                 bookmarkDataIsStale: &stale) else {
            return nil
        }
        if stale {
            saveSelectedDatabase(url)
        }
        return url
    }

    func saveVideoFile(_ url: URL) {
        do {
            let data = try url.bookmarkData(options: [.withSecurityScope],
                                            includingResourceValuesForKeys: nil,
                                            relativeTo: nil)
            var dict = UserDefaults.standard.dictionary(forKey: videoBookmarksKey) as? [String: Data] ?? [:]
            dict[url.path] = data
            UserDefaults.standard.set(dict, forKey: videoBookmarksKey)
        } catch {
            // Best-effort.
        }
    }

    func restoreVideoFile(forPath path: String) -> URL? {
        guard let dict = UserDefaults.standard.dictionary(forKey: videoBookmarksKey) as? [String: Data],
              let data = dict[path] else {
            return nil
        }
        var stale = false
        guard let url = try? URL(resolvingBookmarkData: data,
                                 options: [.withSecurityScope],
                                 relativeTo: nil,
                                 bookmarkDataIsStale: &stale) else {
            return nil
        }
        if stale {
            saveVideoFile(url)
        }
        return url
    }
}

