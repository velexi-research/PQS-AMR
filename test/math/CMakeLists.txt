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

# lsm
add_subdirectory(${CMAKE_CURRENT_SOURCE_DIR}/lsm)
