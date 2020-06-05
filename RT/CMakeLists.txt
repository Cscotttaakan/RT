
set(RT_SOURCES ${CMAKE_CURRENT_LIST_DIR}/src/stb_write_impl.cpp)
set(RT_INCLUDES ${CMAKE_CURRENT_LIST_DIR}/include/stb_image_write.h)

add_executable(RT RT.cpp)
target_include_directories(RT PRIVATE ${CMAKE_SOURCE_DIR}/eigen ${RT_INCLUDES})
target_sources(RT PRIVATE ${RT_SOURCE})