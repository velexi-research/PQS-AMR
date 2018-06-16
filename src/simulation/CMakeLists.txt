# --------------------------------------------------------------
# Customizations for PQS simulation program
# --------------------------------------------------------------

# CMake initialization: No customizations
# Build requirements: No customizations

# --------------------------------------------------------------
# Build parameters
# --------------------------------------------------------------

# Set output directory for executables
set (CMAKE_RUNTIME_OUTPUT_DIRECTORY ${BUILD_EXECUTABLE_DIR})

# Set output directory for sample input files
set (SAMPLE_INPUT_FILE_DIRECTORY
     ${CMAKE_BINARY_DIR}/sample-input-files)

# --------------------------------------------------------------
# Build targets
# --------------------------------------------------------------

# Source files
file(GLOB SOURCE_FILES
     ${CMAKE_CURRENT_SOURCE_DIR}/*.h
     ${CMAKE_CURRENT_SOURCE_DIR}/*.cc)

# Add executable and dependencies
add_executable(pqs ${SOURCE_FILES})
add_dependencies(pqs pqs_kernels samrai)

# Set search path for header files
target_include_directories(pqs
    PUBLIC ${Boost_INCLUDE_DIR}
    PUBLIC ${MPI_CXX_INCLUDE_PATH}
    PUBLIC ${SAMRAI_INCLUDE_DIR})

# Set link libraries
target_link_libraries(pqs
    ${SAMRAI_LIBS}
    ${HDF5_C_LIBRARIES}
    ${MPI_LIBRARY})

# --- Add files to copy directly to build directory

# Create directory for example input files
#file(MAKE_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/example_inputs)

# Copy sample input files to build directory
configure_file(sample-input-files/curvature-2d-circle.input
               ${SAMPLE_INPUT_FILE_DIRECTORY}/curvature-2d-circle.input
               COPYONLY)
configure_file(sample-input-files/curvature-3d-sphere.input
               ${SAMPLE_INPUT_FILE_DIRECTORY}/curvature-3d-sphere.input
               COPYONLY)
