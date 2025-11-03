cmake_minimum_required(VERSION 3.20)
cmake_policy(SET CMP0144 NEW)
cmake_policy(SET CMP0167 NEW)

project(SKNanoAnalyzer VERSION 1.0.0)

set(BOOST_ROOT $ENV{CONDA_PREFIX})
set(Boost_NO_SYSTEM_PATHS ON)
set(Boost_USE_STATIC_LIBS OFF)
set(Boost_USE_MULTITHREADED ON)
set(Boost_USE_STATIC_RUNTIME OFF)

# find packages
find_package(ROOT REQUIRED)
find_package(Boost REQUIRED)
# find_package(Torch REQUIRED)  # Removed: Using Python PyTorch instead of C++ libtorch
find_package(correctionlib REQUIRED)

include(${ROOT_USE_FILE})
include_directories(${Boost_INCLUDE_DIRS})
link_directories(${Boost_LIBRARY_DIRS})

# Add subdirectories
add_subdirectory(external/KinematicFitter)
add_subdirectory(external/RoccoR)
add_subdirectory(DataFormats)
add_subdirectory(AnalyzerTools)
add_subdirectory(Analyzers)

# Fetch external dependencies
include(FetchContent)
FetchContent_Declare(
    yaml-cpp
    GIT_REPOSITORY https://github.com/jbeder/yaml-cpp.git
    GIT_PROGRESS TRUE
    GIT_SHALLOW TRUE
    GIT_TAG master
)
FetchContent_Declare(
    json
    GIT_REPOSITORY https://github.com/nlohmann/json.git
    GIT_PROGRESS TRUE
    GIT_SHALLOW TRUE
    GIT_TAG v3.11.3)

FetchContent_MakeAvailable(yaml-cpp json)

#set cxx standard
set(CMAKE_CXX_STANDARD ${ROOT_CXX_STANDARD})
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)


