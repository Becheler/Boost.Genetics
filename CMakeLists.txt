cmake_minimum_required(VERSION 3.14)
project(BoostGenetics VERSION 0.1.0 LANGUAGES CXX)

# Set C++ standard
set(CMAKE_CXX_STANDARD 11)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

# Header-only library
add_library(boost_genetics INTERFACE)
target_include_directories(boost_genetics INTERFACE
    $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>
    $<INSTALL_INTERFACE:include>
)

# Option to build examples
option(BUILD_EXAMPLES "Build examples" ON)

# Option to build tests
option(BUILD_TESTS "Build tests" ON)

# Option to build benchmarks
option(BUILD_BENCHMARKS "Build benchmarks" OFF)

# Examples
if(BUILD_EXAMPLES)
    add_executable(sequence_example examples/sequence_example.cpp)
    target_link_libraries(sequence_example PRIVATE boost_genetics)

    add_executable(variant_example examples/variant_example.cpp)
    target_link_libraries(variant_example PRIVATE boost_genetics)
endif()

# Tests
if(BUILD_TESTS)
    enable_testing()
    
    add_executable(test_sequence tests/test_sequence.cpp)
    target_link_libraries(test_sequence PRIVATE boost_genetics)
    target_include_directories(test_sequence PRIVATE ${CMAKE_CURRENT_SOURCE_DIR}/tests)
    add_test(NAME test_sequence COMMAND test_sequence)

    add_executable(test_variant tests/test_variant.cpp)
    target_link_libraries(test_variant PRIVATE boost_genetics)
    target_include_directories(test_variant PRIVATE ${CMAKE_CURRENT_SOURCE_DIR}/tests)
    add_test(NAME test_variant COMMAND test_variant)
endif()

# Benchmarks
if(BUILD_BENCHMARKS)
    # Google Benchmark - you may need to install it or use FetchContent
    find_package(benchmark QUIET)
    
    if(NOT benchmark_FOUND)
        # Try to fetch Google Benchmark
        include(FetchContent)
        FetchContent_Declare(
            benchmark
            GIT_REPOSITORY https://github.com/google/benchmark.git
            GIT_TAG v1.8.3
        )
        set(BENCHMARK_ENABLE_TESTING OFF CACHE BOOL "" FORCE)
        FetchContent_MakeAvailable(benchmark)
    endif()

    add_executable(benchmark_sequence benchmarks/benchmark_sequence.cpp)
    target_link_libraries(benchmark_sequence PRIVATE boost_genetics benchmark::benchmark)

    add_executable(benchmark_variant benchmarks/benchmark_variant.cpp)
    target_link_libraries(benchmark_variant PRIVATE boost_genetics benchmark::benchmark)
endif()

# Installation
install(DIRECTORY include/ DESTINATION include)
install(TARGETS boost_genetics EXPORT BoostGeneticsTargets)
install(EXPORT BoostGeneticsTargets
    FILE BoostGeneticsTargets.cmake
    NAMESPACE BoostGenetics::
    DESTINATION lib/cmake/BoostGenetics
)

# Package configuration
include(CMakePackageConfigHelpers)
write_basic_package_version_file(
    "${CMAKE_CURRENT_BINARY_DIR}/BoostGeneticsConfigVersion.cmake"
    VERSION ${PROJECT_VERSION}
    COMPATIBILITY AnyNewerVersion
)

configure_package_config_file(
    "${CMAKE_CURRENT_SOURCE_DIR}/cmake/BoostGeneticsConfig.cmake.in"
    "${CMAKE_CURRENT_BINARY_DIR}/BoostGeneticsConfig.cmake"
    INSTALL_DESTINATION lib/cmake/BoostGenetics
    NO_CHECK_REQUIRED_COMPONENTS_MACRO
)

install(FILES
    "${CMAKE_CURRENT_BINARY_DIR}/BoostGeneticsConfig.cmake"
    "${CMAKE_CURRENT_BINARY_DIR}/BoostGeneticsConfigVersion.cmake"
    DESTINATION lib/cmake/BoostGenetics
)
