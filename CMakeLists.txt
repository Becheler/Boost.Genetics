cmake_minimum_required(VERSION 3.14)
project(BoostGenetics VERSION 0.1.0 LANGUAGES CXX)

# Set C++ standard
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

# Add profiling and optimization flags
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -O3 -g -fno-omit-frame-pointer")

# Optional gperftools for profiling
find_library(PROFILER_LIB profiler HINTS /opt/homebrew/opt/gperftools/lib)
find_path(PROFILER_INCLUDE gperftools/profiler.h HINTS /opt/homebrew/opt/gperftools/include)
if(PROFILER_LIB AND PROFILER_INCLUDE)
    message(STATUS "gperftools found - profiling enabled")
    set(PROFILER_AVAILABLE TRUE)
else()
    message(STATUS "gperftools not found - profiling disabled")
    set(PROFILER_AVAILABLE FALSE)
endif()

# Find dependencies
find_package(Boost COMPONENTS iostreams)
find_package(ZLIB)

# Header-only library
add_library(boost_genetics INTERFACE)
target_include_directories(boost_genetics INTERFACE
    $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>
    $<INSTALL_INTERFACE:include>
)

# Optional BCF support with Boost.IOStreams and zlib
if(Boost_FOUND AND ZLIB_FOUND)
    message(STATUS "Boost.IOStreams and zlib found - enabling BCF support")
    target_compile_definitions(boost_genetics INTERFACE BOOST_GENETICS_HAS_IOSTREAMS)
    target_link_libraries(boost_genetics INTERFACE Boost::iostreams ZLIB::ZLIB)
else()
    message(STATUS "Boost.IOStreams or zlib not found - BCF support disabled (text VCF only)")
endif()

# Option to build examples
option(BUILD_EXAMPLES "Build examples" ON)

# Option to build tests
option(BUILD_TESTS "Build tests" ON)

# Examples
if(BUILD_EXAMPLES)
    add_executable(vcf_example examples/vcf_example.cpp)
    target_link_libraries(vcf_example PRIVATE boost_genetics)
endif()

# Tests
if(BUILD_TESTS)
    enable_testing()
    
    add_executable(test_vcf tests/test_vcf.cpp)
    target_link_libraries(test_vcf PRIVATE boost_genetics)
    target_include_directories(test_vcf PRIVATE ${CMAKE_CURRENT_SOURCE_DIR}/tests)
    add_test(NAME test_vcf COMMAND test_vcf)
    
    add_executable(test_validation tests/test_validation.cpp)
    target_link_libraries(test_validation PRIVATE boost_genetics)
    target_include_directories(test_validation PRIVATE ${CMAKE_CURRENT_SOURCE_DIR}/tests)
    add_test(NAME test_validation COMMAND test_validation WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR})
endif()

# Option to build benchmarks
option(BUILD_BENCHMARKS "Build benchmarks" OFF)

# Benchmarks
if(BUILD_BENCHMARKS)
    add_executable(benchmark_vcf_parsing benchmarks/benchmark_vcf_parsing.cpp)
    target_link_libraries(benchmark_vcf_parsing PRIVATE boost_genetics)
    
    # Add gperftools if available
    if(PROFILER_AVAILABLE)
        target_compile_definitions(benchmark_vcf_parsing PRIVATE USE_GPERFTOOLS)
        target_include_directories(benchmark_vcf_parsing PRIVATE ${PROFILER_INCLUDE})
        target_link_libraries(benchmark_vcf_parsing PRIVATE ${PROFILER_LIB})
    endif()
    
    # Parallel VCF parsing benchmark
    add_executable(benchmark_vcf_parallel benchmarks/benchmark_vcf_parallel.cpp)
    target_link_libraries(benchmark_vcf_parallel PRIVATE boost_genetics)
    find_package(Threads REQUIRED)
    target_link_libraries(benchmark_vcf_parallel PRIVATE Threads::Threads)
    
    # Memory profiling benchmark
    add_executable(benchmark_memory benchmarks/benchmark_memory.cpp)
    target_link_libraries(benchmark_memory PRIVATE boost_genetics)
    
    # Writer performance benchmark
    add_executable(benchmark_writer benchmarks/benchmark_writer.cpp)
    target_link_libraries(benchmark_writer PRIVATE boost_genetics)
    
    # Scalability benchmark
    add_executable(benchmark_scalability benchmarks/benchmark_scalability.cpp)
    target_link_libraries(benchmark_scalability PRIVATE boost_genetics)
    target_link_libraries(benchmark_scalability PRIVATE Threads::Threads)
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
