cmake_minimum_required(VERSION 3.20.0)
project(test_EKF VERSION 0.1.0 LANGUAGES C CXX)

find_package(nlohmann_json REQUIRED)
find_package(LAPACKE REQUIRED
    HINTS /opt/homebrew/opt/lapack/)  # Hint to homebrew installation on Mac

add_executable(test_EKF test_EKF.cc interp.cc)
target_compile_features(test_EKF PRIVATE cxx_std_17)
target_compile_options(test_EKF PRIVATE -Wall -Wextra -Wpedantic)
target_link_libraries(test_EKF PRIVATE nlohmann_json::nlohmann_json)
if (UNIX AND NOT APPLE)  # I don't know why this is needed for Ubuntu
    add_library(lpk INTERFACE)
    target_link_libraries(lpk INTERFACE lapacke gfortran)
    target_link_libraries(test_EKF PRIVATE lpk)
else()
    target_link_libraries(test_EKF PRIVATE lapacke)
endif()