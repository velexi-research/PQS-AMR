# --------------------------------------------------------------
# Customizations for pqs_math_kernels library
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
configure_file(spatial_derivatives_2d.f.in spatial_derivatives_2d.f)
configure_file(spatial_derivatives_3d.f.in spatial_derivatives_3d.f)
configure_file(time_integration_2d.f.in time_integration_2d.f)
configure_file(time_integration_3d.f.in time_integration_3d.f)

# --------------------------------------------------------------
# Build targets
# --------------------------------------------------------------

# Library
file(GLOB SOURCE_FILES
     ${CMAKE_CURRENT_SOURCE_DIR}/*.h
     ${CMAKE_CURRENT_SOURCE_DIR}/*.c
     ${CMAKE_CURRENT_SOURCE_DIR}/*.f
     ${CMAKE_CURRENT_BINARY_DIR}/*.f)
add_library(pqs_math_kernels ${SOURCE_FILES})

# Copy header files to build directory (whenever header files change)
file(GLOB HEADER_FILES ${CMAKE_CURRENT_SOURCE_DIR}/*.h)
foreach(HEADER_FILE ${HEADER_FILES})
    get_filename_component(FILENAME ${HEADER_FILE} NAME)
    configure_file(${HEADER_FILE}
                   ${BUILD_INCLUDE_DIR_PQS}/math/${FILENAME}
                   COPYONLY)
endforeach(HEADER_FILE)

# ----------------------------------------------------------------------
# Installation
# ----------------------------------------------------------------------

# Header files
install(DIRECTORY . DESTINATION ${INSTALL_INCLUDE_DIR_PQS}/math
        FILES_MATCHING PATTERN "*.h")

# Library
install(TARGETS pqs_math_kernels DESTINATION ${INSTALL_LIBRARY_DIR})
