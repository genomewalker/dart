cmake_minimum_required(VERSION 3.18)
project(AncientGenePredictor VERSION 1.0.0 LANGUAGES CXX)
include(CTest)
include(FetchContent)

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

# The pinned libtaph revision. Declared HERE, above configure_file, so it can be stamped into
# version.h -- the guard that enforces it lives further down, where libtaph is added. A DART artifact
# has to be able to say which damage estimator produced it; the guard below makes this SHA true (it
# hard-fails unless the checkout is at exactly this revision AND clean), so stamping it is not a
# guess.
set(LIBTAPH_PINNED_SHA "be1f750637cbc3e9f9909f4075a3a7e3cee9a90d")

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
option(DART_NATIVE_ARCH "Optimize for build machine CPU (-march=native, -ffast-math)" OFF)
if(DART_NATIVE_ARCH)
    set(CMAKE_CXX_FLAGS_RELEASE "-O3 -march=native -mtune=native -DNDEBUG -ffast-math -funroll-loops")
else()
    set(CMAKE_CXX_FLAGS_RELEASE "-O3 -DNDEBUG -funroll-loops")
endif()
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

# rapidgzip (via indexed_bzip2) — in-process parallel gzip decompression.
# Replaces the pigz subprocess; no PATH dependency, true multi-threaded inflate.
# WITH_ISAL ON: rapidgzip inflates via ISA-L igzip (hardware-accelerated), matching
# fqdup's decompression backend. pigz is no longer used at all (removed from
# sequence_io.cpp); the reader path is rapidgzip(ISA-L) -> zlib.
set(LIBRAPIDARCHIVE_BUILD_TESTS OFF CACHE BOOL "" FORCE)
set(LIBRAPIDARCHIVE_BUILD_TOOLS OFF CACHE BOOL "" FORCE)
set(LIBRAPIDARCHIVE_WITH_ISAL   ON  CACHE BOOL "" FORCE)
FetchContent_Declare(
    librapidarchive
    GIT_REPOSITORY https://github.com/mxmlnkn/indexed_bzip2.git
    GIT_TAG        master
    GIT_SHALLOW    TRUE
    UPDATE_DISCONNECTED TRUE
)
FetchContent_MakeAvailable(librapidarchive)

# libtaph: reference-free ancient DNA damage estimation and library-type classification
set(TAPH_BUILD_TESTS OFF CACHE BOOL "" FORCE)
# The damage estimates in every DART artifact come from libtaph, so the libtaph revision IS a build
# input. Pin it and prove the pin, or the numbers are not reproducible.
#
# 14b95f1 is verified green: 6/6 ctest and the oracle's route assertions all pass at this exact
# commit (SLURM 20140502, built from a worktree checked out at it). It is the first revision whose
# build inputs are fully reachable from a commit -- protein_damage.hpp / damage_summary.hpp were
# loose scratch files until 0b590e8.
#
# It is UNPUSHED (local main is ahead of origin), so the FetchContent path below CANNOT resolve it.
# That is exactly why the else() branch is fatal rather than a fallback: silently building against
# GitHub v1.2.1 would emit different numbers under the same version string. Fail loud, never drift.
set(LIBTAPH_LOCAL_PATH "$ENV{REPOS}/libtaph" CACHE PATH "Local libtaph checkout (overrides GitHub)")
if(EXISTS "${LIBTAPH_LOCAL_PATH}/CMakeLists.txt")
    execute_process(
        COMMAND git -C "${LIBTAPH_LOCAL_PATH}" rev-parse HEAD
        OUTPUT_VARIABLE LIBTAPH_HEAD_SHA OUTPUT_STRIP_TRAILING_WHITESPACE ERROR_QUIET)
    execute_process(
        COMMAND git -C "${LIBTAPH_LOCAL_PATH}" status --porcelain
        OUTPUT_VARIABLE LIBTAPH_DIRTY OUTPUT_STRIP_TRAILING_WHITESPACE ERROR_QUIET)
    if(NOT LIBTAPH_HEAD_SHA STREQUAL LIBTAPH_PINNED_SHA)
        message(FATAL_ERROR
            "libtaph at ${LIBTAPH_LOCAL_PATH} is not the pinned revision.\n"
            "  pinned: ${LIBTAPH_PINNED_SHA}\n"
            "  HEAD:   ${LIBTAPH_HEAD_SHA}\n"
            "Damage numbers depend on this. Check it out, or update LIBTAPH_PINNED_SHA deliberately "
            "and re-run the damage validation.")
    endif()
    if(LIBTAPH_DIRTY)
        message(FATAL_ERROR
            "libtaph at ${LIBTAPH_LOCAL_PATH} is at the pinned SHA but the tree is DIRTY.\n"
            "${LIBTAPH_DIRTY}\n"
            "A dirty tree means the built code is not the pinned code, and the SHA stamped into "
            "DART's artifacts would be a lie. Commit or stash first.")
    endif()
    message(STATUS "Using local libtaph, pinned + clean: ${LIBTAPH_PINNED_SHA}")
    add_subdirectory("${LIBTAPH_LOCAL_PATH}" "${CMAKE_BINARY_DIR}/libtaph-local" EXCLUDE_FROM_ALL)
else()
    message(FATAL_ERROR
        "No libtaph checkout at ${LIBTAPH_LOCAL_PATH}, and the pinned revision "
        "${LIBTAPH_PINNED_SHA} is not on any remote -- FetchContent cannot reach it.\n"
        "Falling back to a published tag here would silently build DIFFERENT damage estimates. "
        "Set LIBTAPH_LOCAL_PATH to a checkout at the pinned SHA.")
endif()

# Find zlib for compressed file support
find_package(ZLIB REQUIRED)

# ISA-L: igzip is DART's gzip codec for BOTH read (via rapidgzip) and write. No pigz
# subprocess anywhere -- SequenceWriter/GeneWriter/FastaWriter compress with isal_deflate
# in-process, matching fqdup's backend. Required: a mislabeled or missing isa-l is a
# silent regression to no compression, so fail loud.
find_path(ISAL_INCLUDE_DIR isa-l/igzip_lib.h)
find_library(ISAL_LIB NAMES isal)
if(NOT ISAL_INCLUDE_DIR OR NOT ISAL_LIB)
    message(FATAL_ERROR "ISA-L not found (need isa-l/igzip_lib.h and libisal). "
        "DART's gzip writers require it; there is no pigz fallback.")
endif()
message(STATUS "Found ISA-L: ${ISAL_LIB}")

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
    ${ISAL_INCLUDE_DIR}
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
    src/sequence_io_backend.cpp
    src/dart_columnar.cpp              # .dcx columnar interchange writer (zstd only)
)

# Subcommand sources
set(DART_SUBCOMMAND_SOURCES
    src/cli/dispatcher.cpp
    src/main.cpp                      # cmd_predict
    src/validate.cpp                  # cmd_validate
    src/cli/cmd_sample_damage.cpp     # cmd_sample_damage
    src/cli/cmd_damage_annotate.cpp   # cmd_damage_annotate
    src/cli/cmd_damage_info.cpp       # cmd_damage_info
    src/cli/cmd_damage_profile.cpp    # cmd_damage_profile
    src/cli/cmd_hits2emi.cpp          # hits2emi (TSV to columnar)
    src/cli/cmd_hits2emi_stream.cpp   # hits2emi-stream (mmseqs DB chunks -> one .emi)
    src/cli/cmd_emi_merge.cpp         # emi-merge (union resolved emis for joint EM)
    src/cli/cmd_emi_verify.cpp        # emi-verify (no dup read,target pairs; count integrity)
    src/cli/cmd_simulate.cpp          # simulate (synthetic damage reads)
    src/cli/cmd_eval.cpp              # eval (benchmark accuracy vs truth TSV)
    src/cli/cmd_fit_calib.cpp         # fit-calib (f_true + f_0 coherent-LR model)
)

# Main executable with subcommands
add_executable(dart ${DART_CORE_SOURCES} ${DART_SUBCOMMAND_SOURCES})
target_link_libraries(dart ${ZLIB_LIBRARIES} ${ISAL_LIB} PkgConfig::ZSTD rapidgzip::rapidgzip taph)
target_compile_definitions(dart PRIVATE HAVE_RAPIDGZIP)
if(OpenMP_CXX_FOUND)
    target_link_libraries(dart OpenMP::OpenMP_CXX)
endif()

# Legacy standalone validation tool (for backwards compatibility)
add_executable(dart-validate-standalone src/validate.cpp ${DART_CORE_SOURCES})
target_link_libraries(dart-validate-standalone ${ZLIB_LIBRARIES} ${ISAL_LIB} PkgConfig::ZSTD taph)
target_compile_definitions(dart-validate-standalone PRIVATE DART_STANDALONE_VALIDATE)

# Error analysis tool
add_executable(dart-analyze-errors src/cli/analyze_errors.cpp)
target_link_libraries(dart-analyze-errors ${ZLIB_LIBRARIES})

# Damage prediction benchmark
add_executable(dart-benchmark-damage src/cli/benchmark_damage.cpp)
target_link_libraries(dart-benchmark-damage ${ZLIB_LIBRARIES})

# Streaming EM comparison test
add_executable(test-streaming-em tests/test_streaming_em.cpp ${DART_CORE_SOURCES})
target_link_libraries(test-streaming-em ${ZLIB_LIBRARIES} ${ISAL_LIB} PkgConfig::ZSTD taph)
if(OpenMP_CXX_FOUND)
    target_link_libraries(test-streaming-em OpenMP::OpenMP_CXX)
endif()

# compact_em_map bitwise gate (fusion dual-b8247d43)
add_executable(gate-compact-em tests/gate_compact_em.cpp ${DART_CORE_SOURCES})
target_link_libraries(gate-compact-em ${ZLIB_LIBRARIES} ${ISAL_LIB} PkgConfig::ZSTD taph)
if(OpenMP_CXX_FOUND)
    target_link_libraries(gate-compact-em OpenMP::OpenMP_CXX)
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

    add_executable(
        dart-profile-regressions
        tests/test_profile_regressions.cpp
    )
    add_test(
        NAME profile_regressions
        COMMAND dart-profile-regressions $<TARGET_FILE:dart>
            ${CMAKE_SOURCE_DIR}/examples/tutorial/reads.fq.gz
    )

    add_executable(
        dart-formal-model-tests
        tests/test_formal_model.cpp
        ${DART_CORE_SOURCES}
    )
    target_link_libraries(dart-formal-model-tests ${ZLIB_LIBRARIES} ${ISAL_LIB} PkgConfig::ZSTD taph)
    if(OpenMP_CXX_FOUND)
        target_link_libraries(dart-formal-model-tests OpenMP::OpenMP_CXX)
    endif()
    add_test(
        NAME formal_model
        COMMAND dart-formal-model-tests
    )

    # .dcx writer round-trip (offset math, footer adler32, calib section, obs grain)
    add_executable(
        dart-dcx-roundtrip
        tests/test_dcx_roundtrip.cpp
        src/dart_columnar.cpp
    )
    target_link_libraries(dart-dcx-roundtrip PkgConfig::ZSTD)
    # The round-trip checks (incl. the write_dcx call itself) live inside assert(),
    # so keep asserts active even in Release (-DNDEBUG would strip the write).
    target_compile_options(dart-dcx-roundtrip PRIVATE -UNDEBUG)
    add_test(
        NAME dcx_roundtrip
        COMMAND dart-dcx-roundtrip
    )

    # ORF id uniqueness: the id grammar (legacy 3-field and current 4-field,
    # incl. a two-digit ordinal) plus an end-to-end check that `dart predict`
    # gives two same-frame ORFs of one read two distinct names.
    add_executable(
        dart-orf-id-uniqueness
        tests/test_orf_id_uniqueness.cpp
    )
    add_test(
        NAME orf_id_uniqueness
        COMMAND dart-orf-id-uniqueness $<TARGET_FILE:dart>
                ${CMAKE_CURRENT_SOURCE_DIR}/tests/data/orf_id_golden.tsv
    )

    # EMI v7 format round-trip: blocks, degree-CSR, hot/cold columns, var-length
    # qaln/taln, per-read random access, name-dict, and the two-Reader
    # cache-aliasing regression (thread_local caches keyed on a monotonic uid,
    # not `this`).
    add_executable(
        dart-emi-v7-roundtrip
        tests/emi_v7_roundtrip.cpp
    )
    target_link_libraries(dart-emi-v7-roundtrip PkgConfig::ZSTD)
    add_test(
        NAME emi_v7_roundtrip
        COMMAND dart-emi-v7-roundtrip
    )
endif()

# Install
install(TARGETS dart DESTINATION bin)
