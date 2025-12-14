import SwiftUI

struct ArcScopeCommands: Commands {
    @FocusedValue(\.libraryActions) private var actions

    var body: some Commands {
        CommandGroup(after: .newItem) {
            Button("library.select_database") { actions?.openDatabase() }
                .keyboardShortcut("o", modifiers: [.command])

            Button("library.analyze_new_film") { actions?.analyzeNewFilm() }
                .keyboardShortcut("n", modifiers: [.command])

            Divider()

            Button("library.references") { actions?.openReferences() }
                .keyboardShortcut("r", modifiers: [.command])

            Button("library.export") { actions?.openExport() }
                .keyboardShortcut("e", modifiers: [.command])
                .disabled(!(actions?.canExport ?? false))
        }
    }
}

