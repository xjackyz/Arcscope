import Foundation

enum AppPaths {
    static func applicationSupportDirectory() -> URL {
        let base = FileManager.default.urls(for: .applicationSupportDirectory, in: .userDomainMask).first!
        let dir = base.appendingPathComponent("ArcScope", isDirectory: true)
        try? FileManager.default.createDirectory(at: dir, withIntermediateDirectories: true)
        return dir
    }

    static func defaultDatabaseURL() -> URL {
        applicationSupportDirectory().appendingPathComponent("ArcScope.sqlite", isDirectory: false)
    }
}

