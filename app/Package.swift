// swift-tools-version: 5.9
import PackageDescription

let package = Package(
    name: "ArcScopeApp",
    defaultLocalization: "en",
    platforms: [
        .macOS(.v13)
    ],
    products: [
        .executable(name: "ArcScopeApp", targets: ["ArcScopeApp"])
    ],
    targets: [
        .systemLibrary(
            name: "CAArcScopeBridge",
            path: "Sources/CAArcScopeBridge"
        ),
        .executableTarget(
            name: "ArcScopeApp",
            dependencies: ["CAArcScopeBridge"],
            path: "Sources",
            resources: [
                .process("ArcScopeApp/Resources")
            ],
            linkerSettings: [
                .unsafeFlags([
                    "-L", "../bridge/build",
                    "-larcscope_bridge",
                    "-L", "../bridge/build/core_build",
                    "-larcscope_engine",
                    "-lc++",
                    "-lsqlite3",
                    "-framework", "AVFoundation",
                    "-framework", "CoreMedia",
                    "-framework", "CoreVideo",
                    "-framework", "Accelerate",
                    "-framework", "CoreGraphics",
                    "-framework", "Vision",
                    "-framework", "NaturalLanguage",
                    "-framework", "CoreML"
                ])
            ]
        )
    ]
)
