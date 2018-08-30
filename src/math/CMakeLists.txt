# --------------------------------------------------------------
# Customizations for pqs_math library
# --------------------------------------------------------------

# CMake initialization: No customizations
# Build parameters: No customizations
# Build requirements: No customizations

# --------------------------------------------------------------
# Build targets
# --------------------------------------------------------------

# Library
file(GLOB SOURCE_FILES
     ${CMAKE_CURRENT_SOURCE_DIR}/*.h
     ${CMAKE_CURRENT_SOURCE_DIR}/*.cc)
add_library(pqs_math ${SOURCE_FILES})

# Set search path for header files
target_include_directories(pqs_math
    PUBLIC ${MPI_CXX_INCLUDE_PATH}
    PUBLIC ${SAMRAI_INCLUDE_DIR})

# Copy header files to build directory (whenever header files change)
file(GLOB HEADER_FILES ${CMAKE_CURRENT_SOURCE_DIR}/*.h)
foreach(HEADER_FILE ${HEADER_FILES})
    get_filename_component(FILENAME ${HEADER_FILE} NAME)
    configure_file(${HEADER_FILE}
                   ${BUILD_INCLUDE_DIR_PQS}/math/${FILENAME}
                   COPYONLY)
endforeach(HEADER_FILE)

# --- Sub-targets

# Add library for numerical kernels
add_subdirectory(kernels)

# ----------------------------------------------------------------------
# Installation
# ----------------------------------------------------------------------

# Header files
install(DIRECTORY . DESTINATION ${INSTALL_INCLUDE_DIR_PQS}/math
        FILES_MATCHING PATTERN "*.h")

# Library
install(TARGETS pqs_math DESTINATION ${INSTALL_LIBRARY_DIR})
