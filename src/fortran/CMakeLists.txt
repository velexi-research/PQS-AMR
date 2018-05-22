# --------------------------------------------------------------
# Customizations for pqs_fortran library
# --------------------------------------------------------------

# CMake initialization: No customizations
# Build parameters: No customizations
# Build requirements: No customizations

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
configure_file(spatial_derivatives2d.f.in spatial_derivatives2d.f)
configure_file(spatial_derivatives3d.f.in spatial_derivatives3d.f)

# --------------------------------------------------------------
# Build targets
# --------------------------------------------------------------

# Library
file(GLOB SOURCE_FILES
     ${CMAKE_CURRENT_SOURCE_DIR}/*.h
     ${CMAKE_CURRENT_SOURCE_DIR}/*.c
     ${CMAKE_CURRENT_SOURCE_DIR}/*.f
     ${CMAKE_CURRENT_BINARY_DIR}/*.f)
add_library(pqs_fortran ${SOURCE_FILES})

# Copy header files to build directory
file(GLOB HEADER_FILES ${CMAKE_CURRENT_SOURCE_DIR}/*.h)
file(COPY ${HEADER_FILES} DESTINATION ${CMAKE_BINARY_DIR}/include/pqs)

# --------------------------------------------------------------
# Prepare build directory
# --------------------------------------------------------------

# Create build output directories
file(MAKE_DIRECTORY ${BUILD_LIBRARY_DIR})

# ----------------------------------------------------------------------
# Installation
# ----------------------------------------------------------------------

# Header files
install(DIRECTORY . DESTINATION ${INSTALL_INCLUDE_DIR}
        FILES_MATCHING PATTERN "*.h")

# Library
install(TARGETS pqs_fortran DESTINATION ${INSTALL_LIBRARY_DIR})
