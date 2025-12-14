import SwiftUI

struct SettingsView: View {
    @Environment(\.dismiss) private var dismiss
    @AppStorage("structureUseEmotionProgressReparam") private var structureUseEmotionProgressReparam = false

    var body: some View {
        VStack(spacing: 12) {
            HStack {
                VStack(alignment: .leading, spacing: 4) {
                    Text("settings.title")
                        .font(.title3.weight(.semibold))
                    Text("settings.subtitle")
                        .font(.callout)
                        .foregroundStyle(.secondary)
                }
                Spacer()
                Button("settings.close") { dismiss() }
                    .buttonStyle(.bordered)
            }

            Form {
                Section("settings.section.structure") {
                    Toggle("settings.structure.reparam", isOn: $structureUseEmotionProgressReparam)
                }
            }
            .formStyle(.grouped)
        }
        .padding(18)
        .frame(minWidth: 640, minHeight: 360)
    }
}
