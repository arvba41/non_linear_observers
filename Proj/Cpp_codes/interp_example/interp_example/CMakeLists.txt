cmake_minimum_required(VERSION 3.20.0)
project(spline_test VERSION 0.1.0 LANGUAGES C CXX)

find_package(nlohmann_json REQUIRED)
find_package(LAPACKE REQUIRED
    HINTS /opt/homebrew/opt/lapack/)  # Hint to homebrew installation on Mac

add_executable(spline_test test_interp1.cc interp.cc)
target_compile_features(spline_test PRIVATE cxx_std_17)
target_compile_options(spline_test PRIVATE -Wall -Wextra -Wpedantic)
target_link_libraries(spline_test PRIVATE nlohmann_json::nlohmann_json)
if (UNIX AND NOT APPLE)  # I don't know why this is needed for Ubuntu
    add_library(lpk INTERFACE)
    target_link_libraries(lpk INTERFACE lapacke gfortran)
    target_link_libraries(spline_test PRIVATE lpk)
else()
    target_link_libraries(spline_test PRIVATE lapacke)
endif()
