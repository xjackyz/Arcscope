import SwiftUI

@main
struct ArcScopeApp: App {
    @StateObject private var viewModel = ArcScopeViewModel()

    var body: some Scene {
        WindowGroup {
            LibraryView()
                .environmentObject(viewModel)
                .frame(minWidth: 1200, minHeight: 720)
        }
        .defaultSize(width: 1280, height: 800)
        .commands {
            ArcScopeCommands()
        }

        Settings {
            SettingsView()
        }
    }
}
