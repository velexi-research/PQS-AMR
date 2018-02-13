# ----------------------------------------------------------------------
# Testing customizations
# ----------------------------------------------------------------------

# CMake initialization: No customizations
# Build requirements: No customizations

# ----------------------------------------------------------------------
# Build parameters
# ----------------------------------------------------------------------

# Include search path
include_directories(${GTEST_INCLUDE_DIR})

# Set output directory
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR})

# ----------------------------------------------------------------------
# Build targets
# ----------------------------------------------------------------------

set(AMR_TESTS ${CMAKE_PROJECT_NAME}-tests)

# Source files containing tests
file(GLOB TEST_SOURCE_FILES
     ${PROJECT_SOURCE_DIR}/test/*.cc)

add_executable(${AMR_TESTS} EXCLUDE_FROM_ALL ${TEST_SOURCE_FILES})
include_directories(${gtest_SOURCE_DIR}/include)
target_link_libraries(${AMR_TESTS}
                      gtest
                      gtest_main
                      ${CMAKE_THREAD_LIBS_INIT})

# ----------------------------------------------------------------------
# Test setup
# ----------------------------------------------------------------------
add_test(${AMR_TESTS}
         ${CMAKE_CURRENT_BINARY_DIR}/${AMR_TESTS})
add_dependencies(check ${AMR_TESTS})
