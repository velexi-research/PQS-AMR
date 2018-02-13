# ----------------------------------------------------------------------
# CMake initialization
# ----------------------------------------------------------------------

cmake_minimum_required(VERSION 3.5)
project(snap-amr
        LANGUAGES C CXX Fortran
        VERSION 0.1.0)

# ----------------------------------------------------------------------
# Configuration options
# ----------------------------------------------------------------------

# Build shared libraries
#option(BUILD_SHARED_LIBS "Build shared libraries" ON)

# Floating-point precision
#option(SINGLE_PRECISION "Build single-precision libraries" OFF)
#if (SINGLE_PRECISION)
#    message("-- Setting floating-point precision to 'single'")
#else (SINGLE_PRECISION)
#    message("-- Setting floating-point precision to 'double'")
#endif (SINGLE_PRECISION)

# Boost
set(BOOST "" CACHE FILEPATH "Location of Boost library")

# HDF5
set(HDF5 "" CACHE FILEPATH "Location of HDF5 library")

# MPI
set(MPI "" CACHE FILEPATH "Location of MPI library")

# SAMRAI
set(SAMRAI "" CACHE FILEPATH "Location of SAMRAI library")

# CMake build type
set(CMAKE_BUILD_TYPE Release CACHE STRING "CMake build type")

# ----------------------------------------------------------------------
# Build parameters
# ----------------------------------------------------------------------

# Path variables
set(EXECUTABLE_OUTPUT_PATH ${CMAKE_CURRENT_BINARY_DIR}/bin)
set(LIBRARY_OUTPUT_PATH ${CMAKE_CURRENT_BINARY_DIR}/lib)
set(EXT_PROJECTS_DIR ${PROJECT_SOURCE_DIR}/ext)
set(BINARY_INCLUDE_DIR ${CMAKE_BINARY_DIR}/include)

# Compiler flags
set(CMAKE_C_FLAGS "${CFLAGS} ${CMAKE_C_FLAGS}")
set(CMAKE_CXX_FLAGS ${CMAKE_CXX_FLAGS} -std=c++11)
set(CMAKE_Fortran_FLAGS "${FFLAGS} ${CMAKE_Fortran_FLAGS}")

# Search path for header files
include_directories(${BINARY_INCLUDE_DIR})

# CMake options
set(CMAKE_INCLUDE_CURRENT_DIR ON)

# ----------------------------------------------------------------------
# Build requirements
# ----------------------------------------------------------------------

# --- Software with pre-defined CMake package modules

# Boost
if (BOOST)
    set(BOOST_ROOT ${BOOST})
endif (BOOST)
find_package(Boost REQUIRED)

# Git
find_package(Git)

# HDF5
if (HDF5)

    set(HDF5_ROOT ${HDF5})
    find_package(HDF5 QUIET)
    if (HDF5_FOUND)
        message("-- Found HDF5: " ${HDF5_LIBRARIES} "(found version \"${HDF5_VERSION}\"")
    endif (HDF5_FOUND)
endif (HDF5)

if (NOT HDF5_FOUND)
    unset(HDF5_ROOT)

    # Deal with bugs(?) in HDF5 CMake package
    # by unsetting problematic HDF5 package variables
    unset(HDF5_C_INCLUDE_DIR)
    unset(HDF5_C_INCLUDE_DIRS)

    find_package(HDF5 REQUIRED)
endif (NOT HDF5_FOUND)

# MPI
if (MPI)
    set(CMAKE_PREFIX_PATH_SAVE ${CMAKE_PREFIX_PATH})
    set(CMAKE_PREFIX_PATH ${MPI})
endif (MPI)
find_package(MPI REQUIRED)
if (MPI)
    set(CMAKE_PREFIX_PATH ${CMAKE_PREFIX_PATH_SAVE})
    unset(CMAKE_PREFIX_PATH_SAVE)
endif (MPI)

# Threads
find_package(Threads REQUIRED)

# --- Software with customized checks

# include-what-you-use
find_program(IWYU NAMES include-what-you-use iwyu)

if (IWYU)
    message("-- Found include-what-you-use: ${IWYU}")
    set(CMAKE_CXX_INCLUDE_WHAT_YOU_USE ${IWYU} ${CMAKE_CXX_FLAGS})
else (IWYU)
    message("-- Could not find include-what-you-use.  "
            "Analysis of #includes will be skipped.")
endif (IWYU)

# SAMRAI
add_subdirectory(${EXT_PROJECTS_DIR}/samrai)

# ----------------------------------------------------------------------
# Process template files
# ----------------------------------------------------------------------

## Generate configuration header files
#configure_file(config/LSMLIB_config.h.in
#               ${BINARY_INCLUDE_DIR}/LSMLIB_config.h)

# ----------------------------------------------------------------------
# Build targets
# ----------------------------------------------------------------------

# Testing: SAMRAI LinAdv
add_subdirectory(${PROJECT_SOURCE_DIR}/src/LinAdv)

# examples
#add_subdirectory(examples/toolbox EXCLUDE_FROM_ALL)

# ----------------------------------------------------------------------
# Testing
# ----------------------------------------------------------------------

# --- Build and configure googletest

configure_file(ext/googletest/CMakeLists.txt.in ext/googletest/CMakeLists.txt)
execute_process(
    COMMAND "${CMAKE_COMMAND}" -G "${CMAKE_GENERATOR}" .
        WORKING_DIRECTORY "${CMAKE_BINARY_DIR}/ext/googletest")
execute_process(
    COMMAND "${CMAKE_COMMAND}" --build .
        WORKING_DIRECTORY "${CMAKE_BINARY_DIR}/ext/googletest")

# Prevent GoogleTest from overriding our compiler/linker options
# when building with Visual Studio
set(gtest_force_shared_crt ON CACHE BOOL "" FORCE)

# Add googletest directly to our build. This adds the following targets:
# gtest, gtest_main, gmock and gmock_main
add_subdirectory("${CMAKE_BINARY_DIR}/ext/googletest/googletest-src")
#add_subdirectory("${CMAKE_BINARY_DIR}/ext/googletest/googletest-src"
#                 "${CMAKE_BINARY_DIR}/ext/googletest/googletest-build")

# The gtest/gmock targets carry header search path dependencies automatically
# when using CMake 2.8.11 or later. Otherwise we have to add them here
# ourselves.
if(CMAKE_VERSION VERSION_LESS 2.8.11)
    include_directories("${gtest_SOURCE_DIR}/include"
                        "${gmock_SOURCE_DIR}/include")
endif()

# --- Configure tests

# Enable CMake testing
enable_testing()

# Emulate 'make check'
add_custom_target(check COMMAND ${CMAKE_CTEST_COMMAND})

# Add test directory to build
add_subdirectory(${PROJECT_SOURCE_DIR}/test)

# ----------------------------------------------------------------------
# Installation
# ----------------------------------------------------------------------

# Header files
#file(GLOB HEADER_FILES ${SRCDIR}/*.h)
#install(FILES ${HEADER_FILES} DESTINATION include)

# Library files
#install(TARGETS lsm_parallel DESTINATION lib)

# Executable files
# TODO

# --- Add `uninstall` target

# Generate cmake_uninstall.cmake file
configure_file(
    "${CMAKE_CURRENT_SOURCE_DIR}/cmake_uninstall.cmake.in"
    "${CMAKE_CURRENT_BINARY_DIR}/cmake_uninstall.cmake"
    IMMEDIATE @ONLY
)

# Define `uninstall` target
add_custom_target(uninstall
    COMMAND ${CMAKE_COMMAND}
            -P ${CMAKE_CURRENT_BINARY_DIR}/cmake_uninstall.cmake
)

# ----------------------------------------------------------------------
# Source code distribution
# ----------------------------------------------------------------------

if (GIT_FOUND)
    # Add target to build source code distribution
    set(ARCHIVE_NAME ${CMAKE_PROJECT_NAME}-${PROJECT_VERSION})
    set(ARCHIVE_TARBALL ${ARCHIVE_NAME}.tar.gz)
    add_custom_target(dist
        COMMAND git archive --prefix=${ARCHIVE_NAME}/ HEAD |
                gzip > ${CMAKE_BINARY_DIR}/${ARCHIVE_TARBALL}
        WORKING_DIRECTORY ${PROJECT_SOURCE_DIR}
    )

    # Add target to verify source code distribution
    if (UNIX)
        file(GENERATE OUTPUT check-uninstall.sh
             CONTENT "
remaining=\"`find \${PWD}/${ARCHIVE_NAME}/build/_inst -type f | wc -l`\"
if test \"\${remaining}\" -ne 0; then
    echo \"** Stage directory not empty: \${remaining} file (s) remaining!\"
    exit 1
fi"
        )

        add_custom_target(distcheck
            COMMAND rm -rf ${ARCHIVE_NAME}
            COMMAND gzip -cd ${ARCHIVE_TARBALL} | tar xfz -
            COMMAND cd ${ARCHIVE_NAME} && cmake -H. -Bbuild
            COMMAND cd ${ARCHIVE_NAME}/build && make
            COMMAND cd ${ARCHIVE_NAME}/build && make test
            COMMAND cd ${ARCHIVE_NAME}/build && make DESTDIR=$${PWD}/_inst install
            COMMAND cd ${ARCHIVE_NAME}/build && make DESTDIR=$${PWD}/_inst uninstall
            COMMAND sh check-uninstall.sh
            COMMAND cd ${ARCHIVE_NAME}/build && make clean
            COMMAND rm -rf ${ARCHIVE_NAME}
            COMMAND echo "*** Package ${ARCHIVE_TARBALL} ready for distribution."
        )
        add_dependencies(distcheck dist)

    endif (UNIX)
endif (GIT_FOUND)

# ----------------------------------------------------------------------
# Build preparations
# ----------------------------------------------------------------------

# Create build output directories
file(MAKE_DIRECTORY ${EXECUTABLE_OUTPUT_PATH})
file(MAKE_DIRECTORY ${LIBRARY_OUTPUT_PATH})
