cmake_minimum_required(VERSION 3.12)

# Find packages
find_package(ROOT REQUIRED)
# Find libraries with more robust method
find_library(LHAPDF_LIBRARY NAMES LHAPDF PATHS $ENV{LHAPDF_LIB_DIR}) 
find_library(CORRECTION_LIBRARY NAMES correctionlib PATHS $ENV{CORRECTION_LIB_DIR})
find_library(ONNXRUNTIME_LIBRARY NAMES onnxruntime PATHS $ENV{ONNXRUNTIME_LIB_DIR})

# MLHelper library
add_library(MLHelper SHARED 
    ${CMAKE_CURRENT_SOURCE_DIR}/src/MLHelper.cc 
    ${CMAKE_CURRENT_SOURCE_DIR}/include/MLHelper.h
)


target_include_directories(MLHelper PUBLIC 
    $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>
    $<BUILD_INTERFACE:$ENV{ONNXRUNTIME_INCLUDE_DIR}>
    #${TORCH_INCLUDE_DIRS}
    $<INSTALL_INTERFACE:include>
)

target_link_libraries(MLHelper PUBLIC 
    #${TORCH_LIBRARIES}
    ${ONNXRUNTIME_LIBRARY}
    ${ROOT_LIBRARIES}
)


#install(FILES include/MLHelper.h DESTINATION include)

ROOT_GENERATE_DICTIONARY(G__MLHelper
    include/MLHelper.h
    MODULE MLHelper
    LINKDEF include/MLHelperLinkDef.hpp
)

install(TARGETS MLHelper DESTINATION lib)
set(ROOTMAP "libMLHelper.rootmap")
set(PCM "libMLHelper_rdict.pcm")
install(FILES "${CMAKE_CURRENT_BINARY_DIR}/${ROOTMAP}" "${CMAKE_CURRENT_BINARY_DIR}/${PCM}" DESTINATION lib)
install(FILES include/MLHelper.h DESTINATION include/AnalyzerTools)
# Prepare headers for ROOT dictionary
file(GLOB HEADERS include/*.h)
list(FILTER DICTIONARY_HEADERS EXCLUDE REGEX "MLHelper.h")
foreach(HEADER ${HEADERS})
    file(RELATIVE_PATH REL_HEADER ${CMAKE_CURRENT_SOURCE_DIR} ${HEADER})
    list(APPEND REL_HEADERS ${REL_HEADER})
endforeach()
set(DICTIONARY_HEADERS ${REL_HEADERS})


# Prepare source files
file(GLOB SOURCES src/*.cc)
list(FILTER SOURCES EXCLUDE REGEX "MLHelper.cc")

# Main library target
add_library(AnalyzerTools SHARED
    ${SOURCES}
    #${CMAKE_CURRENT_BINARY_DIR}/G__AnalyzerTools.cxx
)

# Include directories
target_include_directories(AnalyzerTools PUBLIC
    $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>
    $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/../DataFormats/include>
    $<BUILD_INTERFACE:${ROOT_INCLUDE_DIRS}>
    $<BUILD_INTERFACE:$ENV{LHAPDF_INCLUDE_DIR}>
    $<BUILD_INTERFACE:$ENV{CORRECTION_INCLUDE_DIR}>
    $<INSTALL_INTERFACE:include>
)

# Link libraries
target_link_libraries(AnalyzerTools PUBLIC 
    ${ROOT_LIBRARIES}
    ${LHAPDF_LIBRARY}
    ${CORRECTION_LIBRARY}
    DataFormats
    yaml-cpp::yaml-cpp
    MLHelper
)

include_directories($ENV{LHAPDF_INCLUDE_DIR})
# Generate ROOT dictionary
ROOT_GENERATE_DICTIONARY(G__AnalyzerTools
    ${DICTIONARY_HEADERS}
    MODULE AnalyzerTools
    LINKDEF include/AnalyzerToolsLinkDef.hpp
)

# Installation
install(TARGETS AnalyzerTools DESTINATION lib)
install(DIRECTORY include/ DESTINATION include)

# Install ROOT dictionary files if they exist
set(ROOTMAP "libAnalyzerTools.rootmap")
set(PCM "libAnalyzerTools_rdict.pcm")

install(FILES "${CMAKE_CURRENT_BINARY_DIR}/${ROOTMAP}" "${CMAKE_CURRENT_BINARY_DIR}/${PCM}" DESTINATION lib)
