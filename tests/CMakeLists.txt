cmake_minimum_required(VERSION 3.1)
project(filib_tests)

### Conditionally compile certain modules depending on libraries found on the system
list(APPEND CMAKE_MODULE_PATH ${PROJECT_SOURCE_DIR}/../cmake)

### Download data
# igl_download_tutorial_data()

### Output directories
if(MSVC)
  set(CMAKE_RUNTIME_OUTPUT_DIRECTORY_DEBUG ${CMAKE_BINARY_DIR})
  set(CMAKE_RUNTIME_OUTPUT_DIRECTORY_RELEASE ${CMAKE_BINARY_DIR})
else()
  set(CMAKE_RUNTIME_OUTPUT_DIRECTORY "../")
endif()

add_subdirectory("bisection")
add_subdirectory("comp_exp")
add_subdirectory("fi_test")
add_subdirectory("hornerc")
add_subdirectory("hornercpp")
add_subdirectory("xinewton")
add_subdirectory("comp_sin")



