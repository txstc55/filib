get_filename_component(PROJECT_NAME ${CMAKE_CURRENT_SOURCE_DIR} NAME)
project(${PROJECT_NAME})

add_executable(${PROJECT_NAME}_bin comp_exp.c)
target_link_libraries(${PROJECT_NAME}_bin filib)
