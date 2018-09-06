# ----------------------------------------------------------------------
# Testing customizations
# ----------------------------------------------------------------------

# CMake initialization: No customizations
# Build requirements: No customizations

# ----------------------------------------------------------------------
# Build parameters
# ----------------------------------------------------------------------

# Set output directory
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/../bin)

# --------------------------------------------------------------
# Generate Fortran files
# --------------------------------------------------------------

# Floating-point precision
if (SINGLE_PRECISION)

    set(REAL "real*4")
    set(TINY_NONZERO_NUMBER 1.d-35)

else (SINGLE_PRECISION)

    set(REAL "real*8")
    set(TINY_NONZERO_NUMBER 1.d-99)

endif (SINGLE_PRECISION)

# Generate fortran files
configure_file(kernels/kernels_2d.f.in kernels/kernels_2d.f)
configure_file(kernels/kernels_3d.f.in kernels/kernels_3d.f)

# ----------------------------------------------------------------------
# Build targets
# ----------------------------------------------------------------------

set(PQS_AMR_TESTS ${CMAKE_PROJECT_NAME}-lsm-tests)

# Source files containing tests
file(GLOB TEST_SOURCE_FILES
     ${CMAKE_CURRENT_SOURCE_DIR}/*.h
     ${CMAKE_CURRENT_SOURCE_DIR}/*.cc
     ${CMAKE_CURRENT_SOURCE_DIR}/kernels/*.h
     ${CMAKE_CURRENT_BINARY_DIR}/kernels/*.f)

add_executable(${PQS_AMR_TESTS} EXCLUDE_FROM_ALL ${TEST_SOURCE_FILES})
target_link_libraries(${PQS_AMR_TESTS}
                      pqs
                      pqs_kernels
                      lsm
                      lsm_kernels
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
