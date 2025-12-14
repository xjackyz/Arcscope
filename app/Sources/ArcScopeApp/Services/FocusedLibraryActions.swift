import SwiftUI

struct LibraryActions {
    let openDatabase: () -> Void
    let analyzeNewFilm: () -> Void
    let openReferences: () -> Void
    let openExport: () -> Void
    let canExport: Bool
}

private struct LibraryActionsKey: FocusedValueKey {
    typealias Value = LibraryActions
}

extension FocusedValues {
    var libraryActions: LibraryActions? {
        get { self[LibraryActionsKey.self] }
        set { self[LibraryActionsKey.self] = newValue }
    }
}

