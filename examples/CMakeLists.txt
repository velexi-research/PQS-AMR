# --------------------------------------------------------------
# Customizations for pqs-amr program
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
add_executable(pqs-amr ${SOURCE_FILES})
add_dependencies(pqs-amr
                 pqs_utilities
                 pqs pqs_kernels
                 pqs_math pqs_math_kernels
                 samrai)

# Set search path for header files
target_include_directories(pqs-amr
    PUBLIC ${Boost_INCLUDE_DIR}
    PUBLIC ${MPI_CXX_INCLUDE_PATH}
    PUBLIC ${SAMRAI_INCLUDE_DIR})

# Set link libraries
target_link_libraries(pqs-amr
    pqs_utilities
    pqs pqs_kernels
    pqs_math pqs_math_kernels
    ${SAMRAI_LIBS}
    ${HDF5_C_LIBRARIES}
    ${MPI_LIBRARY})

# ----------------------------------------------------------------------
# Installation
# ----------------------------------------------------------------------

# Executables - TODO
#install(TARGETS pqs DESTINATION ${INSTALL_LIBRARY_DIR})
