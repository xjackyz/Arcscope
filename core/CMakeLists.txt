cmake_minimum_required(VERSION 3.16)
project(arcscope_engine LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_POSITION_INDEPENDENT_CODE ON)

	add_library(arcscope_engine
	    src/Normalization.cpp
	    src/CurveEngine.cpp
	    src/MetadataEngine.cpp
	    src/Diagnostics.cpp
	    src/SQLiteStore.cpp
	    src/FilmEngine.cpp
	    src/Utility.cpp
    src/features/SoundFeatures.cpp
    src/features/ColorFeatures.cpp
    src/features/FaceFeatures.cpp
    src/features/MotionFeatures.cpp
    src/features/InfoFeatures.cpp
    src/segmentation/StructureAnalyzer.cpp
)

target_include_directories(arcscope_engine PUBLIC include)
find_library(SQLITE3_LIBRARY sqlite3 REQUIRED)
target_link_libraries(arcscope_engine PUBLIC ${SQLITE3_LIBRARY})

