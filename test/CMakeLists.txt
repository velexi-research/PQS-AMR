# ----------------------------------------------------------------------
# Testing customizations
# ----------------------------------------------------------------------

# CMake initialization: No customizations
# Build requirements: No customizations

# ----------------------------------------------------------------------
# Build parameters
# ----------------------------------------------------------------------

# Set output directory
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/bin)

# ----------------------------------------------------------------------
# Build targets
# ----------------------------------------------------------------------

set(PQS_AMR_TESTS ${CMAKE_PROJECT_NAME}-tests)

# Source files containing tests
file(GLOB TEST_SOURCE_FILES
     ${PROJECT_SOURCE_DIR}/test/*/*.cc)

add_executable(${PQS_AMR_TESTS} EXCLUDE_FROM_ALL ${TEST_SOURCE_FILES})
target_link_libraries(${PQS_AMR_TESTS}
                      pqs
                      ${SAMRAI_LIBS}
                      ${HDF5_C_LIBRARY_hdf5}
                      ${MPI_LIBRARY}
                      gtest
                      gtest_main
                      ${CMAKE_THREAD_LIBS_INIT})

# ----------------------------------------------------------------------
# Test setup
# ----------------------------------------------------------------------
add_test(${PQS_AMR_TESTS}
    ${CMAKE_RUNTIME_OUTPUT_DIRECTORY}/${PQS_AMR_TESTS})
add_dependencies(check ${PQS_AMR_TESTS})
