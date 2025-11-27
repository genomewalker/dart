cmake_minimum_required(VERSION 3.18)
project(AncientGenePredictor VERSION 0.2.1 LANGUAGES CXX)

# Git-based version detection
find_package(Git QUIET)
if(GIT_FOUND)
    execute_process(
        COMMAND ${GIT_EXECUTABLE} describe --tags --always --dirty
        WORKING_DIRECTORY ${CMAKE_SOURCE_DIR}
        OUTPUT_VARIABLE GIT_VERSION
        OUTPUT_STRIP_TRAILING_WHITESPACE
        ERROR_QUIET
        RESULT_VARIABLE GIT_RESULT
    )
    execute_process(
        COMMAND ${GIT_EXECUTABLE} rev-parse --short HEAD
        WORKING_DIRECTORY ${CMAKE_SOURCE_DIR}
        OUTPUT_VARIABLE GIT_COMMIT
        OUTPUT_STRIP_TRAILING_WHITESPACE
        ERROR_QUIET
    )
    if(GIT_RESULT EQUAL 0 AND GIT_VERSION)
        # Extract version from tag (e.g., v0.2.1 -> 0.2.1)
        string(REGEX REPLACE "^v" "" AGP_VERSION "${GIT_VERSION}")
    else()
        set(AGP_VERSION "${PROJECT_VERSION}-${GIT_COMMIT}")
    endif()
else()
    set(AGP_VERSION "${PROJECT_VERSION}")
    set(GIT_COMMIT "unknown")
endif()

message(STATUS "AGP Version: ${AGP_VERSION}")

# Generate version header
configure_file(
    ${CMAKE_SOURCE_DIR}/cmake/version.h.in
    ${CMAKE_BINARY_DIR}/include/agp/version.h
    @ONLY
)
include_directories(${CMAKE_BINARY_DIR}/include)

# Set C++20 standard
set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

# Compiler flags for optimization
set(CMAKE_CXX_FLAGS_RELEASE "-O3 -march=native -mtune=native -DNDEBUG -ffast-math -funroll-loops")
set(CMAKE_CXX_FLAGS_DEBUG "-g -O0 -Wall -Wextra -Wpedantic")

# Enable AVX2/AVX512 if available
include(CheckCXXCompilerFlag)
check_cxx_compiler_flag("-mavx2" COMPILER_SUPPORTS_AVX2)
check_cxx_compiler_flag("-mavx512f" COMPILER_SUPPORTS_AVX512)

if(COMPILER_SUPPORTS_AVX512)
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -mavx512f -mavx512bw")
    add_compile_definitions(USE_AVX512)
elseif(COMPILER_SUPPORTS_AVX2)
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -mavx2")
    add_compile_definitions(USE_AVX2)
endif()

# Options
option(USE_OPENMP "Enable OpenMP parallelization" ON)
option(PGO_GENERATE "Build with profile generation" OFF)
option(PGO_USE "Build with profile-guided optimization" OFF)

# PGO flags
if(PGO_GENERATE)
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -fprofile-generate=${CMAKE_BINARY_DIR}/pgo")
    set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} -fprofile-generate=${CMAKE_BINARY_DIR}/pgo")
endif()

if(PGO_USE)
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -fprofile-use=${CMAKE_BINARY_DIR}/pgo -fprofile-correction")
    set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} -fprofile-use=${CMAKE_BINARY_DIR}/pgo")
endif()

# Find OpenMP
if(USE_OPENMP)
    find_package(OpenMP)
    if(OpenMP_CXX_FOUND)
        set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
    endif()
endif()

# Find zlib for compressed file support
find_package(ZLIB REQUIRED)

# Include directories
include_directories(
    ${PROJECT_SOURCE_DIR}/include
    ${ZLIB_INCLUDE_DIRS}
)

# Source files
set(AGP_SOURCES
    src/main.cpp
    src/types.cpp
    src/damage_model.cpp
    src/sequence_io.cpp
    src/codon_tables.cpp
    src/scoring.cpp
    src/damage_detection.cpp
    src/frame_selector.cpp
)

# Main executable
add_executable(agp ${AGP_SOURCES})
target_link_libraries(agp ${ZLIB_LIBRARIES})

if(OpenMP_CXX_FOUND)
    target_link_libraries(agp OpenMP::OpenMP_CXX)
endif()

# Unified validation tool
add_executable(agp-validate src/validate.cpp)
target_link_libraries(agp-validate ${ZLIB_LIBRARIES})

# Hexamer extraction tool
add_executable(extract-hexamers scripts/extract_refseq_hexamers.cpp)
target_link_libraries(extract-hexamers ${ZLIB_LIBRARIES})
if(OpenMP_CXX_FOUND)
    target_link_libraries(extract-hexamers OpenMP::OpenMP_CXX)
endif()

# Install
install(TARGETS agp agp-validate extract-hexamers DESTINATION bin)
