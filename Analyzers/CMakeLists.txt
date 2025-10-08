# find LHAPDF and CorrectionLib
set(CORRECTION_INCLUDE_DIR $ENV{CORRECTION_INCLUDE_DIR})
set(LHAPDF_INCLUDE_DIR $ENV{LHAPDF_INCLUDE_DIR})

# Define the library with source files and the generated dictionary
file(GLOB SOURCES src/*.cc)
add_library(Analyzers SHARED
    ${SOURCES}
    #${CMAKE_CURRENT_BINARY_DIR}/G__Analyzers.cxx
)

target_include_directories(Analyzers PRIVATE
    $<BUILD_INTERFACE:${LHAPDF_INCLUDE_DIR}>
    $<BUILD_INTERFACE:${CORRECTION_INCLUDE_DIR}>
    $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>
    $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/../AnalyzerTools/include>
    $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/../DataFormats/include>
    $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/../external/KinematicFitter/include>
    $<BUILD_INTERFACE:${ROOT_INCLUDE_DIRS}>
    $<BUILD_INTERFACE:${TORCH_INCLUDE_DIRS}>
    $<INSTALL_INTERFACE:include>
)

# Generate the ROOT dictionary before defining the library
file(GLOB HEADERS include/*.h)
foreach(HEADER ${HEADERS})
    file(RELATIVE_PATH REL_HEADER ${CMAKE_CURRENT_SOURCE_DIR} ${HEADER})
    list(APPEND REL_HEADERS ${REL_HEADER})
endforeach()
set(DICTIONARY_HEADERS ${REL_HEADERS})

ROOT_GENERATE_DICTIONARY(G__Analyzers 
    ${DICTIONARY_HEADERS}
    MODULE Analyzers
    LINKDEF include/AnalyzersLinkDef.hpp
)



# Link ROOT libraries
target_link_libraries(Analyzers PUBLIC ${ROOT_LIBRARIES})
target_link_libraries(Analyzers PUBLIC DataFormats)
target_link_libraries(Analyzers PUBLIC AnalyzerTools)
target_link_libraries(Analyzers PUBLIC KinFitter)
#target_link_libraries(Analyzers PUBLIC MLHelper)
target_link_libraries(Analyzers PUBLIC nlohmann_json::nlohmann_json)
target_link_libraries(Analyzers PUBLIC ${TORCH_LIBRARIES})

# Specify where to install the library
install(TARGETS Analyzers DESTINATION lib)


# Iinstall .rootmap and .pcm files to the lib
set(ROOTMAP "libAnalyzers.rootmap")
set(PCM "libAnalyzers_rdict.pcm")
message(STATUS "ROOTMAP: ${ROOTMAP}")
message(STATUS "PCM: ${PCM}")
install(FILES "${CMAKE_CURRENT_BINARY_DIR}/${ROOTMAP}" "${CMAKE_CURRENT_BINARY_DIR}/${PCM}" DESTINATION lib)
install(FILES ${HEADERS} DESTINATION include)
