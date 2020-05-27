cmake_minimum_required(VERSION 3.10)
project(tiebrush)

set (CMAKE_STATIC_LINKER_FLAGS "-Wl,--as-needed.-lcurl")

include(ExternalProject)

set(CMAKE_CXX_STANDARD 11)

include("cmake/zlib.cmake")
include("cmake/htslib.cmake")

include(FetchContent)

FetchContent_Declare(
        gclib
        GIT_REPOSITORY https://github.com/gpertea/gclib.git
)

FetchContent_MakeAvailable(gclib)

message(STATUS "gclib source dir: ${gclib_SOURCE_DIR}")
message(STATUS "gclib binary dir: ${gclib_BINARY_DIR}")

include_directories(${gclib_SOURCE_DIR})

add_executable(tiebrush src/tiebrush.cpp
                        src/GSam.cpp
                        src/GSam.h
                        src/tmerge.cpp
                        src/tmerge.h
                        ${gclib_SOURCE_DIR}/GStr.h
                        ${gclib_SOURCE_DIR}/GStr.cpp
                        ${gclib_SOURCE_DIR}/GArgs.h
                        ${gclib_SOURCE_DIR}/GArgs.cpp
                        ${gclib_SOURCE_DIR}/GBase.h
                        ${gclib_SOURCE_DIR}/GBase.cpp)

add_dependencies(tiebrush htslib zlib)
target_link_libraries(tiebrush ${htslib_LIB} ${zlib_LIB} -lbz2 -llzma -pthread -lcrypto -lcurl)