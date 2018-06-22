# --------------------------------------------------------------
# Customizations for pqs library
# --------------------------------------------------------------

# CMake initialization: No customizations
# Build requirements: No customizations

# --------------------------------------------------------------
# Build parameters
# --------------------------------------------------------------

# Set output directory for executables
set (CMAKE_RUNTIME_OUTPUT_DIRECTORY ${BUILD_EXECUTABLE_DIR})

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

# Copy header files to build directory (whenever header files change)
file(GLOB HEADER_FILES ${CMAKE_CURRENT_SOURCE_DIR}/*.h)
foreach(HEADER_FILE ${HEADER_FILES})
    get_filename_component(FILENAME ${HEADER_FILE} NAME)
    configure_file(${HEADER_FILE}
                   ${BUILD_INCLUDE_DIR_PQS}/pqs/${FILENAME}
                   COPYONLY)
endforeach(HEADER_FILE)

# ----------------------------------------------------------------------
# Installation
# ----------------------------------------------------------------------

# Header files
install(DIRECTORY . DESTINATION ${INSTALL_INCLUDE_DIR_PQS}/pqs
        FILES_MATCHING PATTERN "*.h")

# Library
install(TARGETS pqs DESTINATION ${INSTALL_LIBRARY_DIR})
