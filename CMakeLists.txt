cmake_minimum_required(VERSION 3.18)
project(AncientGenePredictor VERSION 1.0.0 LANGUAGES CXX)
include(CTest)

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
        # Extract version from tag (e.g., v0.3.0 -> 0.3.0)
        string(REGEX REPLACE "^v" "" DART_VERSION "${GIT_VERSION}")
    else()
        set(DART_VERSION "${PROJECT_VERSION}-${GIT_COMMIT}")
    endif()
else()
    set(DART_VERSION "${PROJECT_VERSION}")
    set(GIT_COMMIT "unknown")
endif()

message(STATUS "DART Version: ${DART_VERSION}")

# Generate version header
configure_file(
    ${CMAKE_SOURCE_DIR}/cmake/version.h.in
    ${CMAKE_BINARY_DIR}/include/dart/version.h
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

# Enable AVX2 for SIMD optimizations (x86_64 only)
if(CMAKE_SYSTEM_PROCESSOR MATCHES "x86_64|AMD64|x86")
    include(CheckCXXCompilerFlag)
    check_cxx_compiler_flag("-mavx2" COMPILER_SUPPORTS_AVX2)
    if(COMPILER_SUPPORTS_AVX2)
        set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -mavx2")
        add_compile_definitions(USE_AVX2)
    endif()
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

# Find zstd for columnar compression (required)
find_package(PkgConfig REQUIRED)
pkg_check_modules(ZSTD REQUIRED IMPORTED_TARGET libzstd)
add_compile_definitions(HAVE_ZSTD)
message(STATUS "ZSTD found: ${ZSTD_VERSION}")
if(FALSE)

    # TBB disabled - linking issues on this system
    # pkg_check_modules(TBB tbb)
    # if(TBB_FOUND)
    #     add_compile_definitions(HAVE_TBB)
    #     message(STATUS "TBB found: ${TBB_VERSION}")
    # endif()
endif()

# Include directories
include_directories(
    ${PROJECT_SOURCE_DIR}/include
    ${PROJECT_SOURCE_DIR}/src
    ${ZLIB_INCLUDE_DIRS}
)

# Source files
# Hexamer data is header-only (constexpr arrays in include/dart/*_hexamer_table.hpp)

# Core library sources (shared between subcommands)
set(DART_CORE_SOURCES
    src/cli/args.cpp
    src/cli/subcommand.cpp
    src/damage/profile.cpp
    src/damage/probability.cpp
    src/damage/analyzer.cpp
    src/types.cpp
    src/sequence_io.cpp
    src/codon_tables.cpp
    src/scoring.cpp
    src/sample_profile.cpp
    src/damage_probability.cpp
    src/damage_scoring.cpp
    src/frame_selector.cpp
    src/adaptive_damage.cpp
    src/damage_index_writer.cpp
    src/damage_index_reader.cpp
    src/reference_stats.cpp
    src/em_reassign.cpp
    src/em_index.cpp
    src/columnar_index.cpp
    src/coverage_em.cpp
)

# Subcommand sources
set(DART_SUBCOMMAND_SOURCES
    src/cli/dispatcher.cpp
    src/main.cpp                      # cmd_predict
    src/validate.cpp                  # cmd_validate
    src/cli/cmd_sample_damage.cpp     # cmd_sample_damage
    src/cli/cmd_damage_annotate.cpp   # cmd_damage_annotate
    src/cli/cmd_damage_profile.cpp    # cmd_damage_profile
    src/cli/cmd_hits2emi.cpp          # hits2emi (TSV to columnar)
)

# Main executable with subcommands
add_executable(dart ${DART_CORE_SOURCES} ${DART_SUBCOMMAND_SOURCES})
target_link_libraries(dart ${ZLIB_LIBRARIES} PkgConfig::ZSTD)
if(OpenMP_CXX_FOUND)
    target_link_libraries(dart OpenMP::OpenMP_CXX)
endif()

# Legacy standalone validation tool (for backwards compatibility)
add_executable(dart-validate-standalone src/validate.cpp ${DART_CORE_SOURCES})
target_link_libraries(dart-validate-standalone ${ZLIB_LIBRARIES} PkgConfig::ZSTD)
target_compile_definitions(dart-validate-standalone PRIVATE DART_STANDALONE_VALIDATE)

# Error analysis tool
add_executable(dart-analyze-errors src/cli/analyze_errors.cpp)
target_link_libraries(dart-analyze-errors ${ZLIB_LIBRARIES})

# Damage prediction benchmark
add_executable(dart-benchmark-damage src/cli/benchmark_damage.cpp)
target_link_libraries(dart-benchmark-damage ${ZLIB_LIBRARIES})

# Streaming EM comparison test
add_executable(test-streaming-em tests/test_streaming_em.cpp ${DART_CORE_SOURCES})
target_link_libraries(test-streaming-em ${ZLIB_LIBRARIES} PkgConfig::ZSTD)
if(OpenMP_CXX_FOUND)
    target_link_libraries(test-streaming-em OpenMP::OpenMP_CXX)
endif()

if(BUILD_TESTING)
    add_executable(
        dart-cli-args-tests
        tests/test_cli_args.cpp
        src/cli/args.cpp
    )
    add_test(
        NAME cli_args
        COMMAND dart-cli-args-tests
    )

    add_executable(
        dart-benchmark-fixture-test
        tests/test_benchmark_fixture.cpp
    )
    add_test(
        NAME benchmark_fixture
        COMMAND dart-benchmark-fixture-test ${CMAKE_SOURCE_DIR}/examples/tutorial/data/benchmark_summary.tsv
    )

    add_executable(
        dart-damage-annotate-regressions
        tests/test_damage_annotate_regressions.cpp
    )
    add_test(
        NAME damage_annotate_regressions
        COMMAND dart-damage-annotate-regressions $<TARGET_FILE:dart>
    )
endif()

# Install
install(TARGETS dart DESTINATION bin)
