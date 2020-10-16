cmake_minimum_required(VERSION 3.8)
project(filib)

################################################################################
set(CMAKE_CXX_STANDARD 11)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

### Configuration
set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} ${CMAKE_CURRENT_SOURCE_DIR}/cmake)

include(Warnings)
include(UseColors)
include(${PROJECT_NAME}Utils)


# inclusion FILIB
add_subdirectory(filib)
prepend_current_path(FILIB_SOURCES)
filib_copy_headers(${FILIB_SOURCES})
filib_set_source_group(${FILIB_SOURCES})
# message(WARNING "src files: ${FILIB_SOURCES}")
add_library(${PROJECT_NAME} STATIC ${FILIB_SOURCES})
target_include_directories(${PROJECT_NAME} PUBLIC ${PROJECT_BINARY_DIR}/include)
target_link_libraries(${PROJECT_NAME} PRIVATE warnings::all)

if(CMAKE_SOURCE_DIR STREQUAL CMAKE_CURRENT_SOURCE_DIR)
	set(FILIB_TOPLEVEL_PROJECT ON)
else()
	set(FILIB_TOPLEVEL_PROJECT OFF)
endif()

if(FILIB_TOPLEVEL_PROJECT)
	# Unit tests
	include(CTest)
	enable_testing()
	add_subdirectory(tests)
endif()
