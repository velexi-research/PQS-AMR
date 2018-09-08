# ----------------------------------------------------------------------
# CMake initialization
# ----------------------------------------------------------------------

cmake_minimum_required(VERSION 3.5)
project(pqs-amr
        LANGUAGES C CXX Fortran
        VERSION 0.1.0)

# ----------------------------------------------------------------------
# Configuration options
# ----------------------------------------------------------------------

# Build shared libraries
#option(BUILD_SHARED_LIBS "Build shared libraries" ON)

# Floating-point precision
option(SINGLE_PRECISION "Build single-precision libraries" OFF)
if (SINGLE_PRECISION)
    message("-- Setting floating-point precision to 'single'")
else (SINGLE_PRECISION)
    message("-- Setting floating-point precision to 'double'")
endif (SINGLE_PRECISION)

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

# --- Path variables

# Source directory
set(EXTERNAL_PACKAGES_SOURCE_DIR ${PROJECT_SOURCE_DIR}/ext)

# Build directory
set(BUILD_EXECUTABLE_DIR ${CMAKE_CURRENT_BINARY_DIR}/bin)
set(BUILD_LIBRARY_DIR ${CMAKE_CURRENT_BINARY_DIR}/lib)
set(LIBRARY_OUTPUT_PATH ${BUILD_LIBRARY_DIR})
set(BUILD_INCLUDE_DIR ${CMAKE_BINARY_DIR}/include)
set(BUILD_INCLUDE_DIR_PQS ${BUILD_INCLUDE_DIR}/PQS)

# Install directory
set(INSTALL_LIBRARY_DIR lib)
set(INSTALL_INCLUDE_DIR include)
set(INSTALL_INCLUDE_DIR_PQS ${INSTALL_INCLUDE_DIR}/PQS)

# Search paths for header files
include_directories(${BUILD_INCLUDE_DIR})

# --- Compiler flags

set(CMAKE_C_FLAGS "${CFLAGS} ${CMAKE_C_FLAGS}")
set(CMAKE_CXX_FLAGS ${CMAKE_CXX_FLAGS} -std=c++11)
set(CMAKE_Fortran_FLAGS "${FFLAGS} ${CMAKE_Fortran_FLAGS}")

# --- CMake options

set(CMAKE_INCLUDE_CURRENT_DIR ON)

# ----------------------------------------------------------------------
# Build requirements
# ----------------------------------------------------------------------

# --- Software with pre-defined CMake package modules

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
add_subdirectory(${EXTERNAL_PACKAGES_SOURCE_DIR}/samrai)

# ----------------------------------------------------------------------
# Process template files
# ----------------------------------------------------------------------

# Floating-point precision
if (SINGLE_PRECISION)

    set(REAL "float")
    set(ZERO_TOL 1.e-5)
    set(REAL_MAX "FLT_MAX")
    set(REAL_MIN "FLT_MIN")
    set(REAL_EPSILON "FLT_EPSILON")

else (SINGLE_PRECISION)

    set(PQS_DOUBLE_PRECISION ON)

    set(REAL "double")
    set(ZERO_TOL 1.e-11)
    set(REAL_MAX "DBL_MAX")
    set(REAL_MIN "DBL_MIN")
    set(REAL_EPSILON "DBL_EPSILON")

endif (SINGLE_PRECISION)

# Generate configuration header file
configure_file(config/PQS_config.h.in ${BUILD_INCLUDE_DIR_PQS}/PQS_config.h)

# ----------------------------------------------------------------------
# Build targets
# ----------------------------------------------------------------------

# PQS
#add_subdirectory(${PROJECT_SOURCE_DIR}/src/)

# pqs_utilities
add_subdirectory(${PROJECT_SOURCE_DIR}/src/utilities)

# pqs_math
add_subdirectory(${PROJECT_SOURCE_DIR}/src/math)

# pqs
add_subdirectory(${PROJECT_SOURCE_DIR}/src/pqs)

# PQS examples
add_subdirectory(${PROJECT_SOURCE_DIR}/examples)

# Sample input files
add_subdirectory(${PROJECT_SOURCE_DIR}/src/sample-input-files)

# ----------------------------------------------------------------------
# Prepare build directory
# ----------------------------------------------------------------------

# Create build output directories
file(MAKE_DIRECTORY ${BUILD_EXECUTABLE_DIR})
file(MAKE_DIRECTORY ${BUILD_LIBRARY_DIR})

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
add_subdirectory(
    "${CMAKE_BINARY_DIR}/ext/googletest/src/googletest"
    "${CMAKE_BINARY_DIR}/ext/googletest/build"
    EXCLUDE_FROM_ALL)

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

# Header files - handled by subdirectory CMakeList.txt files
install(FILES ${BUILD_INCLUDE_DIR_PQS}/PQS_config.h
        DESTINATION ${INSTALL_INCLUDE_DIR_PQS})

# Library files - handled by subdirectory CMakeList.txt files

# Executable files - handled by subdirectory CMakeList.txt files

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
        file(GENERATE OUTPUT ${CMAKE_BINARY_DIR}/check-uninstall.sh
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
