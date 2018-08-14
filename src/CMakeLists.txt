# --------------------------------------------------------------
# Customizations for PQS package
# --------------------------------------------------------------

# CMake initialization: No customizations
# Build parameters: No customizations
# Build requirements: No customizations

# --------------------------------------------------------------
# Build targets
# --------------------------------------------------------------

# Copy header files to build directory (whenever header files change)
file(GLOB HEADER_FILES ${CMAKE_CURRENT_SOURCE_DIR}/*.h)
foreach(HEADER_FILE ${HEADER_FILES})
    get_filename_component(FILENAME ${HEADER_FILE} NAME)
    configure_file(${HEADER_FILE}
                   ${BUILD_INCLUDE_DIR_PQS}/${FILENAME}
                   COPYONLY)
endforeach(HEADER_FILE)

# ----------------------------------------------------------------------
# Installation
# ----------------------------------------------------------------------

# Header files
install(DIRECTORY . DESTINATION ${INSTALL_INCLUDE_DIR_PQS}
        FILES_MATCHING PATTERN "*.h")
