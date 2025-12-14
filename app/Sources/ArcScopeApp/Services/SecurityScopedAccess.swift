import Foundation

enum SecurityScopedAccess {
    static func withAccess<T>(_ url: URL, _ body: () throws -> T) rethrows -> T {
        let didStart = url.startAccessingSecurityScopedResource()
        defer {
            if didStart {
                url.stopAccessingSecurityScopedResource()
            }
        }
        return try body()
    }

    static func withAccessAsync<T>(_ url: URL, _ body: () async throws -> T) async rethrows -> T {
        let didStart = url.startAccessingSecurityScopedResource()
        defer {
            if didStart {
                url.stopAccessingSecurityScopedResource()
            }
        }
        return try await body()
    }

    static func withAccessAsync<T>(_ urls: [URL], _ body: () async throws -> T) async rethrows -> T {
        let started = urls.map { $0.startAccessingSecurityScopedResource() }
        defer {
            for (url, didStart) in zip(urls, started) where didStart {
                url.stopAccessingSecurityScopedResource()
            }
        }
        return try await body()
    }
}

