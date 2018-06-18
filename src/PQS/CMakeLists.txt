# --------------------------------------------------------------
# Customizations for PQS Solver library
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

 # Library
file(GLOB SOURCE_FILES
     ${CMAKE_CURRENT_SOURCE_DIR}/*.h
     ${CMAKE_CURRENT_SOURCE_DIR}/*.cc)
add_library(pqs ${SOURCE_FILES})

# Set search path for header files
target_include_directories(pqs
    PUBLIC ${Boost_INCLUDE_DIR}
    PUBLIC ${MPI_CXX_INCLUDE_PATH}
    PUBLIC ${SAMRAI_INCLUDE_DIR})

# Copy header files to build directory
file(GLOB HEADER_FILES ${CMAKE_CURRENT_SOURCE_DIR}/*.h)
file(COPY ${HEADER_FILES} DESTINATION ${CMAKE_BINARY_DIR}/include/pqs)

# ----------------------------------------------------------------------
# Installation
# ----------------------------------------------------------------------

# Header files
install(DIRECTORY . DESTINATION ${INSTALL_INCLUDE_DIR}
        FILES_MATCHING PATTERN "*.h")

# Library
install(TARGETS pqs DESTINATION ${INSTALL_LIBRARY_DIR})
