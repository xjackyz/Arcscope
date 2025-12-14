import CAArcScopeBridge
import Foundation

enum VideoAnalysisError: LocalizedError {
    case invalidInput
    case analysisFailed(Int)

    var errorDescription: String? {
        switch self {
        case .invalidInput:
            return L10n.string("error.video_invalid")
        case .analysisFailed(let code):
            return L10n.format("error.analysis_failed", code)
        }
    }
}

final class VideoAnalysisService {
    static let shared = VideoAnalysisService()

    func analyzeVideo(inputURL: URL,
                      databaseURL: URL,
                      title: String,
                      year: Int,
                      director: String,
                      progress: @escaping (Double) -> Void) async throws {
        try await SecurityScopedAccess.withAccessAsync([inputURL, databaseURL]) {
            guard FileManager.default.fileExists(atPath: inputURL.path) else {
                throw VideoAnalysisError.invalidInput
            }
            try await Task.detached(priority: .userInitiated) {
                try self.performAnalysis(inputURL: inputURL,
                                         databaseURL: databaseURL,
                                         title: title,
                                         year: year,
                                         director: director,
                                         progress: progress)
            }.value
        }
    }

    private func performAnalysis(inputURL: URL,
                                 databaseURL: URL,
                                 title: String,
                                 year: Int,
                                 director: String,
                                 progress: @escaping (Double) -> Void) throws {
        let useEmotionProgressReparam = UserDefaults.standard.bool(forKey: "structureUseEmotionProgressReparam")
        let options: UInt32 = useEmotionProgressReparam ? UInt32(ARCSCOPE_ANALYZE_OPT_STRUCTURE_U_REPARAM) : 0

        let box = ProgressBox(handler: progress)
        let unmanaged = Unmanaged.passRetained(box)
        defer { unmanaged.release() }
        let status = inputURL.path.withCString { inputPath in
            databaseURL.path.withCString { dbPath in
                title.withCString { titlePtr in
                    director.withCString { directorPtr in
                        arcscope_analyze_film_with_options(inputPath,
                                                           dbPath,
                                                           nil,
                                                           titlePtr,
                                                           Int32(year),
                                                           directorPtr,
                                                           options,
                                                           { value, context in
                                                               guard let context else { return }
                                                               let box = Unmanaged<ProgressBox>.fromOpaque(context).takeUnretainedValue()
                                                               box.report(Double(value))
                                                           },
                                                           unmanaged.toOpaque())
                    }
                }
            }
        }
        guard status == 0 else {
            throw VideoAnalysisError.analysisFailed(Int(status))
        }
    }

    private final class ProgressBox {
        private let handler: (Double) -> Void
        init(handler: @escaping (Double) -> Void) {
            self.handler = handler
        }
        func report(_ value: Double) {
            handler(min(max(value, 0.0), 1.0))
        }
    }
}
