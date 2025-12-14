cmake_minimum_required(VERSION 3.16)
project(arcscope_bridge LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

add_subdirectory(../core ${CMAKE_BINARY_DIR}/core_build)

add_library(arcscope_bridge
    STATIC
    src/arcscope_bridge.mm
    VideoFeatureExtractor.mm
    AudioAnalysisEngine.mm
    FaceTrackingEngine.mm
    SubtitleNLP.mm)
set_source_files_properties(src/arcscope_bridge.mm VideoFeatureExtractor.mm AudioAnalysisEngine.mm FaceTrackingEngine.mm SubtitleNLP.mm PROPERTIES COMPILE_FLAGS "-fobjc-arc")
target_include_directories(arcscope_bridge
    PUBLIC
    ${CMAKE_CURRENT_SOURCE_DIR}/include
    ${CMAKE_CURRENT_SOURCE_DIR}
    ../core/include)
set(OBJC_FRAMEWORKS
    "-framework AVFoundation"
    "-framework CoreMedia"
    "-framework CoreVideo"
    "-framework Foundation"
    "-framework Accelerate"
    "-framework CoreGraphics"
    "-framework Vision"
    "-framework CoreML"
    "-framework NaturalLanguage"
)
target_link_libraries(arcscope_bridge PUBLIC arcscope_engine ${OBJC_FRAMEWORKS})
