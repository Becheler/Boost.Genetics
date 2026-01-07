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
