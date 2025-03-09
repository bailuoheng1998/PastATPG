cmake_minimum_required(VERSION 3.16)


# Compile as C++17
set(CMAKE_CXX_STANDARD 17)
# Require C++11 support
set(CMAKE_CXX_STANDARD_REQUIRED ON)
# We expect external libraries to be linked statically
set(BUILD_SHARED_LIBS OFF)

#set(CMAKE_CXX_COMPILER "/usr/bin/g++-9")
#set(CMAKE_C_COMPILER  "/usr/bin/gcc-9")

project(ICTest
        VERSION 1.1.0
        DESCRIPTION "CASTest DFT Tool"
        LANGUAGES C CXX   # Support C & CPP language
        )
find_package(OpenMP)
find_package(ZLIB QUIET)
#set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -O0 -g")
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -O3")
# Include guards
if (PROJECT_BINARY_DIR STREQUAL PROJECT_SOURCE_DIR)
    message(FATAL_ERROR "The binary directory of CMake cannot be the same as source directory!"
            "Please consider use 'cmake -B build' to specify a different binary directory."
            "Otherwise ZENO may fail to build, or make Git hard to exclude binary files.")
endif ()

# Default is build release
# Debug is build by cmake -B build -DCMAKE_BUILD_TYPE=Debug
if (NOT CMAKE_BUILD_TYPE)
    set(CMAKE_BUILD_TYPE Release)
endif ()
# print ${CMAKE_BUILD_TYPE} info
message(STATUS "Build directory: ${PROJECT_BINARY_DIR}")
message(STATUS "ICTest Build Type: ${CMAKE_BUILD_TYPE}")


# CMake modules
add_subdirectory(src)
add_subdirectory(test)
add_subdirectory(thirdparty)

# Enable unit test
enable_testing()