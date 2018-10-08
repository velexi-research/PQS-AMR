# --------------------------------------------------------------
# Customizations for example application
# --------------------------------------------------------------

# CMake initialization: No customizations
# Build requirements: No customizations
# Installation: No customizations

# --------------------------------------------------------------
# Build parameters
# --------------------------------------------------------------

# Set output directory for executables
set (CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CURRENT_BINARY_DIR})

# Set output directory for sample input files
set (SAMPLE_INPUT_FILE_DIRECTORY
     ${CMAKE_BINARY_DIR}/sample-input-files)

# --------------------------------------------------------------
# Process template files
# --------------------------------------------------------------

# Process Fortran files for numerical kernels
add_subdirectory(kernels)

# --------------------------------------------------------------
# Build targets
# --------------------------------------------------------------

# Source files
file(GLOB SOURCE_FILES
     ${CMAKE_CURRENT_SOURCE_DIR}/*.h
     ${CMAKE_CURRENT_SOURCE_DIR}/*.cc
     ${CMAKE_CURRENT_SOURCE_DIR}/kernels/*.h
     ${CMAKE_CURRENT_BINARY_DIR}/kernels/*.f)

# Add executable and dependencies
add_executable(pqs-amr-2d-no-solid-phase ${SOURCE_FILES})
add_dependencies(pqs-amr-2d-no-solid-phase
                 pqs_utilities
                 pqs pqs_kernels
                 pqs_math pqs_math_kernels
                 samrai)

# Set search path for header files
target_include_directories(pqs-amr-2d-no-solid-phase
    PUBLIC ${MPI_CXX_INCLUDE_PATH}
    PUBLIC ${SAMRAI_INCLUDE_DIR})

# Set link libraries
target_link_libraries(pqs-amr-2d-no-solid-phase
    pqs_utilities
    pqs pqs_kernels
    pqs_math pqs_math_kernels
    ${SAMRAI_LIBS}
    ${HDF5_C_LIBRARIES}
    ${MPI_LIBRARY})

# Add files to copy directly to build directory
configure_file(sample.conf sample.conf COPYONLY)
