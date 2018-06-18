# --------------------------------------------------------------
# Customizations for sample input files
# --------------------------------------------------------------

# CMake initialization: No customizations
# Build requirements: No customizations

# --------------------------------------------------------------
# Build parameters
# --------------------------------------------------------------

# Set output directory for sample input files
set (SAMPLE_INPUT_FILE_DIRECTORY ${CMAKE_BINARY_DIR}/sample-input-files)

# --------------------------------------------------------------
# Build targets
# --------------------------------------------------------------

# Copy sample input files to build directory
configure_file(curvature-2d-circle.input
               ${SAMPLE_INPUT_FILE_DIRECTORY}/curvature-2d-circle.input
               COPYONLY)
configure_file(curvature-3d-sphere.input
               ${SAMPLE_INPUT_FILE_DIRECTORY}/curvature-3d-sphere.input
               COPYONLY)
