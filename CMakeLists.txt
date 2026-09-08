cmake_minimum_required(VERSION 3.16)
project(hampr VERSION 0.1.0 LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

if(NOT CMAKE_BUILD_TYPE)
    set(CMAKE_BUILD_TYPE Release CACHE STRING "Build type" FORCE)
endif()

set(HAMPR_ENABLE_GPU "OFF" CACHE STRING "GPU acceleration backend: OFF, SYCL, or OPENCL")
set_property(CACHE HAMPR_ENABLE_GPU PROPERTY STRINGS OFF SYCL OPENCL)

if(HAMPR_ENABLE_GPU)
    message(STATUS "GPU acceleration backend: ${HAMPR_ENABLE_GPU}")
endif()

option(HAMPR_BUILD_TESTS "Build unit tests" ON)
option(HAMPR_BUILD_EXAMPLES "Build example applications" ON)
option(HAMPR_BUILD_PYTHON "Build Python bindings" OFF)

find_package(PkgConfig REQUIRED)
pkg_check_modules(FFTW3 IMPORTED_TARGET fftw3)
find_package(Eigen3 REQUIRED)

set(HAMPR_SOURCES
    src/dsp/channel_isolator.cpp
    src/dsp/beamformer.cpp
    src/dsp/clutter_canceler.cpp
    src/dsp/detector.cpp
    src/dsp/doa_estimator.cpp
    src/dsp/metric_extractor.cpp
    src/dsp/array_geometry.cpp
    src/dsp/synchronizer.cpp
    src/dsp/fusion.cpp
    src/dsp/localization.cpp
    src/dsp/multi_site_pipeline.cpp
    src/dsp/multi_static_correlator.cpp
    src/dsp/multi_static_fusion.cpp
    src/dsp/pipeline.cpp
    src/io/text_data_source.cpp
    src/io/console_sink.cpp
    src/io/file_sink.cpp
    src/io/multi_site_data_source.cpp
    src/io/streaming_source.cpp
    src/pipeline/graph.cpp
    src/pipeline/multi_site_graph.cpp
    src/accel/fftw_backend.cpp
    src/accel/eigen_backend.cpp
)

if(HAMPR_ENABLE_GPU STREQUAL "SYCL")
    find_package(SYCL REQUIRED)
    list(APPEND HAMPR_SOURCES src/accel/sycl_fft_backend.cpp src/accel/sycl_linalg_backend.cpp)
elseif(HAMPR_ENABLE_GPU STREQUAL "OPENCL")
    find_package(OpenCL REQUIRED)
    list(APPEND HAMPR_SOURCES src/accel/opencl_fft_backend.cpp)
endif()

add_library(hampr ${HAMPR_SOURCES})

target_include_directories(hampr
    PUBLIC
        ${CMAKE_CURRENT_SOURCE_DIR}/include
)

target_link_libraries(hampr
    PUBLIC
        PkgConfig::FFTW3
        Eigen3::Eigen
)

if(HAMPR_ENABLE_GPU STREQUAL "SYCL")
    target_compile_definitions(hampr PRIVATE HAMPR_GPU_BACKEND=SYCL)
    target_link_libraries(hampr PRIVATE SYCL::SYCL)
elseif(HAMPR_ENABLE_GPU STREQUAL "OPENCL")
    target_compile_definitions(hampr PRIVATE HAMPR_GPU_BACKEND=OPENCL)
    target_link_libraries(hampr PRIVATE OpenCL::OpenCL)
endif()

target_compile_options(hampr PRIVATE -Wall -Wextra -Wpedantic)

if(HAMPR_BUILD_EXAMPLES)
    add_executable(benchmark_app examples/benchmark_app.cpp)
    target_link_libraries(benchmark_app PRIVATE hampr)
    set_target_properties(benchmark_app PROPERTIES
        RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/bin
    )
endif()

if(HAMPR_BUILD_TESTS)
    enable_testing()
    add_executable(test_core tests/test_core.cpp)
    target_link_libraries(test_core PRIVATE hampr)
    add_test(NAME test_core COMMAND test_core)
    add_executable(test_dsp tests/test_dsp.cpp)
    target_link_libraries(test_dsp PRIVATE hampr)
    add_test(NAME test_dsp COMMAND test_dsp)
    add_executable(test_dsp_extended tests/test_dsp_extended.cpp)
    target_link_libraries(test_dsp_extended PRIVATE hampr)
    add_test(NAME test_dsp_extended COMMAND test_dsp_extended)
    add_executable(test_regression tests/test_regression.cpp)
    target_link_libraries(test_regression PRIVATE hampr)
    add_test(NAME test_regression COMMAND test_regression)
    add_executable(test_accel tests/test_accel.cpp)
    target_link_libraries(test_accel PRIVATE hampr)
    add_test(NAME test_accel COMMAND test_accel)
    add_executable(test_streaming tests/test_streaming.cpp)
    target_link_libraries(test_streaming PRIVATE hampr)
    add_test(NAME test_streaming COMMAND test_streaming)
    add_executable(test_multi_site tests/test_multi_site.cpp)
    target_link_libraries(test_multi_site PRIVATE hampr)
    add_test(NAME test_multi_site COMMAND test_multi_site)
    add_executable(test_phase8 tests/test_phase8.cpp)
    target_link_libraries(test_phase8 PRIVATE hampr)
    add_test(NAME test_phase8 COMMAND test_phase8)
endif()

install(TARGETS hampr
    LIBRARY DESTINATION lib
    ARCHIVE DESTINATION lib
)
install(DIRECTORY include/hampr
    DESTINATION include
)
