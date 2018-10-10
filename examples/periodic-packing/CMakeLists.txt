
# CMake initialization: No customizations
# Build requirements: No customizations
# Installation: No customizations

# --------------------------------------------------------------
# Build parameters
# --------------------------------------------------------------

# Set output directory for executables
set (CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CURRENT_BINARY_DIR})

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
add_executable(pqs-amr-periodic-packing ${SOURCE_FILES})
add_dependencies(pqs-amr-periodic-packing
                 pqs_utilities
                 pqs pqs_kernels
                 pqs_math pqs_math_kernels
                 samrai)

# Set search path for header files
target_include_directories(pqs-amr-periodic-packing
    PUBLIC ${MPI_CXX_INCLUDE_PATH}
    PUBLIC ${SAMRAI_INCLUDE_DIR})

# Set link libraries
target_link_libraries(pqs-amr-periodic-packing
    pqs_utilities
    pqs pqs_kernels
    pqs_math pqs_math_kernels
    ${SAMRAI_LIBS}
    ${HDF5_C_LIBRARIES}
    ${MPI_LIBRARY})

# Add files to copy directly to build directory
configure_file(sample-2D.conf sample-2D.conf COPYONLY)
#configure_file(sample-3D.conf sample-3D.conf COPYONLY)
