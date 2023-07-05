cmake_minimum_required(VERSION 3.5)

set(CMAKE_SKIP_BUILD_RPATH FALSE)
set(CMAKE_BUILD_WITH_INSTALL_RPATH TRUE)
set(CMAKE_INSTALL_RPATH "$ORIGIN/../lib")
set(CMAKE_INSTALL_RPATH_USE_LINK_PATH TURE)

project(D2C)

set(CMAKE_BUILD_TYPE Release)

add_subdirectory(d2c)
message(STATUS "Compile D2C")

